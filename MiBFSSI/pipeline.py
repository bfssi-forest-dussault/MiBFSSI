#!/usr/bin/env python3

__version__ = "0.0.1"
__author__ = "Forest Dussault"
__email__ = "forest.dussault@canada.ca"

import os
import click
import shutil
import logging
import pandas as pd
import multiprocessing

from pathlib import Path
from dataclasses import dataclass
from subprocess import Popen, PIPE

script = os.path.basename(__file__)
logger = logging.getLogger()
logging.basicConfig(
    format=f'\033[92m \033[1m {script} %(levelname)-2s:\033[0m %(message)s ',
    level=logging.INFO)


def print_version(ctx, param, value):
    if not value or ctx.resilient_parsing:
        return
    logging.info(f"Version: {__version__}")
    logging.info(f"Author: {__author__}")
    logging.info(f"Email: {__email__}")
    quit()


@dataclass
class Sample(object):
    def __lt__(self, other):
        """Allows for Sample objects to be sorted"""
        return self.sample_id < other.sample_id

    # Mandatory attributes
    sample_id: str
    r1: Path
    r2: Path
    outdir: Path

    # Optional attributes
    reference_genome: Path = None
    taxid: str = None
    taxname: str = None
    bamfile: Path = None
    mapping_stats: Path = None
    mapping_stats_df: pd.DataFrame = None
    assembly: Path = None


@click.command()
@click.option('-i', '--inputdir',
              type=click.Path(exists=True),
              required=True,
              default=None,
              help='Root directory to store all output files')
@click.option('-o', '--outdir',
              type=click.Path(exists=False),
              required=True,
              default=None,
              help='Root directory to store all output files')
@click.option('-f', '--forward-id',
              type=click.STRING,
              required=False,
              default='_R1',
              help='Pattern to recognize forward reads in input directory. Defaults to "R1".')
@click.option('-r', '--reverse-id',
              type=click.STRING,
              required=False,
              default='_R2',
              help='Pattern to recognize forward reads in input directory. Defaults to "R2".')
@click.option('-r', '--reference',
              type=click.Path(exists=True),
              required=True,
              default=None,
              help='Path to a reference .FASTA to use instead of the automatically acquired '
                   'references from sendsketch.sh')
@click.option('-t', '--threads',
              type=click.INT,
              required=False,
              default=multiprocessing.cpu_count() - 1,
              help='Number of threads to dedicate to parallelizable steps of the pipeline.'
                   'Will take all available threads - 1 by default.')
@click.option('--snippy',
              help='Specify this flag to run Snippy against each sample. Useful if you just want to call variants.',
              is_flag=True,
              default=False)
@click.option('--bbmap',
              help='Specify this flag to run BBMap.sh against each sample to generate sorted BAM files and coverage '
                   'statistics.',
              is_flag=True,
              default=False)
@click.option('--nullarbor',
              help='Specify this flag to run Nullarbor against each sample. This will run the full pipeline on each '
                   'sample. https://github.com/tseemann/nullarbor',
              is_flag=True,
              default=False)
@click.option('--version',
              help='Specify this flag to print the version and exit.',
              is_flag=True,
              is_eager=True,
              callback=print_version,
              expose_value=False)
def pipeline(inputdir, outdir, forward_id, reverse_id, reference, threads, snippy, bbmap, nullarbor):
    logging.info("Starting MiBFSSI Pipeline")

    # Path setup
    inputdir = Path(inputdir)
    outdir = Path(outdir)
    reference = Path(reference)

    # Reference validation
    if reference is not None and reference.suffix == ".gz":
        logging.error("ERROR: Please provide an uncompressed reference FASTA file.")
        quit()

    # Pipeline flags
    pipeline_flags = {"snippy": snippy, "bbmap": bbmap, "nullarbor": nullarbor}
    active_flags = [y for x, y in pipeline_flags.items() if y is True]

    # Alignment pipeline validation
    if len(active_flags) > 1:
        logging.error(f"ERROR: Please only specify one of the following: {pipeline_flags.keys()}")
        quit()
    elif len(active_flags) == 0:
        logging.error(f"ERROR: Please specify one of the following {pipeline_flags.keys()}")
        quit()

    # Output directory validation
    try:
        os.makedirs(str(outdir), exist_ok=False)
    except FileExistsError:
        logging.error("ERROR: Output directory already exists.")
        quit()

    # Pair up FASTQ files and prepare a dictionary of samples
    fastq_dict = get_sample_dictionary(inputdir, forward_id, reverse_id)

    # Put everything into Sample objects
    sample_list = list()
    for sample_id, reads in fastq_dict.items():
        try:
            sample = Sample(sample_id=sample_id,
                            r1=reads[0],
                            r2=reads[1],
                            outdir=outdir / sample_id)
            sample_list.append(sample)
            if reference is not None:
                sample.reference_genome = reference
        except KeyError:
            pass
    sample_list = sorted(sample_list)

    if reference is None:
        # Populate taxid, taxname, and reference genome for each Sample object
        taxid_reference_retrieval(sample_list=sample_list, outdir=outdir)

    # Run analyses
    pipeline_flags = {"snippy": snippy, "bbmap": bbmap, "nullarbor": nullarbor}
    logging.info(f"Running {[x for x, y in pipeline_flags.items() if y is True]}")
    for sample in sample_list:
        try:
            # Generate BAM + sorted BAM with bbmap
            if bbmap:
                (sample.bamfile, sample.mapping_stats) = call_bbmap(fwd_reads=sample.r1,
                                                                    rev_reads=sample.r2,
                                                                    reference=sample.reference_genome,
                                                                    outdir=sample.outdir,
                                                                    sample_id=sample.sample_id,
                                                                    threads=threads)

                # Parse bbmap output and drop it into a dataframe
                sample.mapping_stats_df = parse_genome_results(sample.mapping_stats, sample_id=sample.sample_id)

                # Create final mapping stats report for every sample
                logging.info("Collating mapping statistics")
                mapping_stats_list = [sample.mapping_stats_df for sample in sample_list]
                mapping_stats_report = combine_dataframes(mapping_stats_list)
                mapping_stats_report.to_csv(outdir / "project_mapping_stats.csv", sep=",", index=None)

                # Qualimap
                sample.qualimap_dir = call_qualimap(bamfile=sample.bamfile, outdir=sample.outdir, threads=threads)
            # Call snippy against each sample
            elif snippy:
                sample.snippy_dir = call_snippy(fwd_reads=sample.r1,
                                                rev_reads=sample.r2,
                                                reference=sample.reference_genome,
                                                outdir=sample.outdir,
                                                threads=threads,
                                                prefix=sample.sample_id)
                (sample.snippy_summary, sample.snippy_vcf,
                 sample.snippy_consensus, sample.bamfile) = parse_snippy(snippy_dir=sample.snippy_dir)

                # Qualimap
                sample.qualimap_dir = call_qualimap(bamfile=sample.bamfile, outdir=sample.outdir, threads=threads)

        except KeyError as e:
            logging.error(f"ERROR: Could not call bbmap.sh on {sample.sample_id}. Traceback:")
            logging.error(e)

    # Run nullarbor on all samples simultaneously
    if nullarbor:
        sample_file = prepare_nullarbor_sample_file(samples=sample_list, outdir=outdir)
        call_nullarbor(project_name=outdir.name, reference=reference,
                       outdir=outdir / 'nullarbor', samples=sample_file, threads=threads)


def taxid_reference_retrieval(sample_list, outdir):
    # Grab taxIDs for each sample
    logging.info("Running sendsketch.sh to find closest reference genome for each sample")
    for sample in sample_list:
        sample.taxid, sample.taxname = call_sendsketch(fwd_reads=sample.r1, rev_reads=sample.r2)
        logging.info(f"{sample.sample_id}: {sample.taxname} (taxid {sample.taxid})")

    # Prepare list of unique taxids to request for ncbi-genome-download
    taxid_download_list = list(set([x.taxid for x in sample_list]))

    # Download the genomes, return dict with key: taxid and value: path to reference
    logging.info("Calling ncbi-genome-download to retrieve reference genome(s)")
    reference_genome_dict = taxid_reference_dict(taxids=taxid_download_list, outdir=outdir)

    # Update sample objects with corresponding reference genomes
    for sample in sample_list:
        try:
            sample.reference_genome = reference_genome_dict[sample.taxid]
        except KeyError:
            pass


def call_sendsketch(fwd_reads: Path, rev_reads: Path) -> tuple:
    tmpreads = fwd_reads.parent / 'tmpfile.fastq.gz'
    tmptxt = fwd_reads.parent / 'tmpfile.txt'

    # Merge reads
    cmd = f"bbmerge.sh in1={fwd_reads} in2={rev_reads} out={tmpreads} overwrite=true"
    run_subprocess(cmd)

    # Call sendsketch
    cmd = f"sendsketch.sh in={tmpreads} out={tmptxt} " \
          f"overwrite=true address=refseq"
    run_subprocess(cmd)

    # Delete merged reads
    tmpreads.unlink()

    # Retrieve taxid
    df = pd.read_csv(tmptxt, delimiter="\t", skiprows=2)
    taxid = df['TaxID'][0]
    taxname = df['taxName'][0]

    # Delete tmptxt
    tmptxt.unlink()

    return str(taxid), str(taxname)


def call_fuse(reference: Path) -> Path:
    fused_reference = reference.with_suffix(".fused.fna.gz")  # TODO: Make this less brittle
    cmd = f"fuse.sh in={reference} out={fused_reference} overwrite=true"
    run_subprocess(cmd)
    return fused_reference


def call_bbmap(fwd_reads: Path, rev_reads: Path, reference: Path, outdir: Path, sample_id: str, threads: int) -> tuple:
    logging.debug(f"Running bbmap on {sample_id} against the reference {reference.name}")
    # Create sample folder
    os.makedirs(str(outdir), exist_ok=True)

    # Setup output files for bbmap
    outstats = outdir / (sample_id + ".mapping_stats.txt")
    outbam = outdir / (sample_id + ".bam")
    outsortedbam = outdir / (sample_id + "_sorted.bam")

    # Call bbmap
    cmd = f"bbmap.sh in={fwd_reads} in2={rev_reads} ref={reference} outm={outbam} covstats={outstats} " \
          f"rebuild=t overwrite=true threads={threads} bamscript=bs.sh; sh bs.sh"
    run_subprocess(cmd)

    # Remove unsorted bam
    outbam.unlink()

    return outsortedbam, outstats


def call_qualimap(bamfile: Path, outdir: Path, threads: int) -> Path:
    outdir = outdir / "qualimap"
    logging.debug(f"Running Qualimap on {bamfile.name}")
    cmd = f"qualimap bamqc -bam {bamfile} -outdir {outdir} -nt {threads} -nw 3000"
    run_subprocess(cmd)
    return outdir


def call_snippy(fwd_reads: Path, rev_reads: Path, reference: Path, outdir: Path, prefix: str, threads: int) -> Path:
    """
    Note that this is still using version 3.x because it was conda installed. Should be updated to 4.x soon...
    https://github.com/tseemann/snippy
    """
    outdir = outdir / "snippy"
    cmd = f"snippy --cpus {threads} --outdir {outdir} --ref {reference} --prefix {prefix} " \
          f"--R1 {fwd_reads} --R2 {rev_reads}"
    run_subprocess(cmd)
    return outdir


def parse_snippy(snippy_dir: Path) -> tuple:
    snippy_summary = list(snippy_dir.glob('*.tab'))[0]
    snippy_vcf = list(snippy_dir.glob('*.vcf'))[0]
    snippy_consensus = list(snippy_dir.glob('*.aligned.fa'))[0]
    snippy_bam = list(snippy_dir.glob('*.bam'))[0]
    return snippy_summary, snippy_vcf, snippy_consensus, snippy_bam


def taxid_reference_dict(taxids: list, outdir: Path):
    reference_genome_dict = dict()
    for taxid in taxids:
        # Download the closest matching genome from RefSeq
        reference_genome = retrieve_reference_genome(taxid, outdir)

        # Use fuse.sh to fuse contigs (arbitrary padding of 300 Ns for any gaps)
        reference_genome = call_fuse(reference_genome)

        # Update reference genome dictionary
        reference_genome_dict[taxid] = reference_genome

    # Cleanup junk download folder from ncbi-genome-download
    shutil.rmtree(outdir / "refseq")
    return reference_genome_dict


def taxid_list_to_text(taxid_list: list, outdir: Path) -> Path:
    taxfile = Path(outdir) / "taxids.txt"
    with open(taxfile, mode="w") as outfile:
        for s in taxid_list:
            outfile.write(f"{s}\n")
    return taxfile


def retrieve_reference_genome(taxid_string: str, outdir: Path):
    # Retrieve reference genomes
    cmd = f"ncbi-genome-download " \
          f"-o {outdir} --taxid {taxid_string} -F fasta bacteria"
    run_subprocess(cmd)
    logging.info(cmd)
    fasta = list(outdir.glob("*/**/*.fna.gz"))[0]  # Recursively glob for FASTA files, grab first (and only) fasta
    fasta_path = Path(outdir / fasta.name)
    fasta.rename(fasta_path)
    return fasta_path


def get_sample_dictionary(directory: Path, forward_id: str, reverse_id: str) -> dict:
    """
    Chains several functions together to create a sample dictionary with unique/valid sample IDs as keys
    and paths to forward and reverse reads as values
    :param directory: Path to a directory containing .fastq.gz files
    :param forward_id: ID indicating forward read in filename (e.g. _R1)
    :param reverse_id: ID indicating reverse read in filename (e.g. _R2)
    :return: Validated sample dictionary with sample_ID:R1,R2 structure
    """
    fastq_file_list = retrieve_fastqgz(directory)
    sample_id_list = retrieve_unique_sampleids(fastq_file_list)
    sample_dictionary = populate_sample_dictionary(sample_id_list, fastq_file_list, forward_id, reverse_id)
    logging.info(f"Successfully paired {len(sample_dictionary)} samples")
    return sample_dictionary


def retrieve_fastqgz(directory: Path) -> [Path]:
    """
    :param directory: Path to folder containing output from MiSeq run
    :return: LIST of all .fastq.gz files in directory
    """
    fastq_file_list = list(directory.glob("*.f*q*"))
    return fastq_file_list


def retrieve_unique_sampleids(fastq_file_list: [Path]) -> list:
    """
    :param fastq_file_list: List of fastq.gz filepaths generated by retrieve_fastqgz()
    :return: List of valid OLC Sample IDs
    """
    # Iterate through all of the fastq files and grab the sampleID, append to list
    sample_id_list = list()
    for f in fastq_file_list:
        sample_id = f.name.split('_')[0]
        sample_id_list.append(sample_id)

    # Get unique sample IDs
    sample_id_list = list(set(sample_id_list))
    return sample_id_list


def prepare_nullarbor_sample_file(samples: [Sample], outdir: Path) -> Path:
    outfile = outdir / 'Nullarbor_SampleSheet.tab'

    # Create dataframe with relevant sample object metadata
    d = []
    for sample in samples:
        d.append({'SampleID': sample.sample_id, 'R1': sample.r1, 'R2': sample.r2})
    df = pd.DataFrame(d)

    # Reorder columns
    df = df[['SampleID', 'R1', 'R2']]

    # Export to .tab file
    df.to_csv(outfile, sep="\t", index=None, header=False)
    return outfile


def call_nullarbor(project_name: str, reference: Path, samples: Path, outdir: Path, threads: int):
    """
    Let nullarbor+SKESA do all of the hard work. Currently not available for install via conda, but should be soon.

    https://github.com/tseemann/nullarbor
    """
    cmd = f"nullarbor.pl --name {project_name} --ref {reference} --input {samples} --outdir {outdir} --cpus {threads}"
    run_subprocess(cmd)


def populate_sample_dictionary(sample_id_list: list, fastq_file_list: [Path], forward_id: str,
                               reverse_id: str) -> dict:
    """
    :param sample_id_list: List of unique Sample IDs generated by retrieve_unique_sampleids()
    :param fastq_file_list: List of fastq.gz file paths generated by retrieve_fastqgz()
    :param forward_id: ID indicating forward read in filename (e.g. _R1)
    :param reverse_id: ID indicating reverse read in filename (e.g. _R2)
    :return: dictionary with each Sample ID as a key and the read pairs as values
    """
    # Find file pairs for each unique sample ID
    sample_dictionary = {}
    for sample_id in sample_id_list:
        read_pair = get_readpair(sample_id, fastq_file_list, forward_id, reverse_id)
        if read_pair is not None:
            sample_dictionary[sample_id] = read_pair
        else:
            pass
    return sample_dictionary


def index_bam(bam: Path):
    cmd = f"samtools index {bam}"
    run_subprocess(cmd)


def mpileup(bam: Path, reference: Path) -> Path:
    logging.info("Generating mpileup file")

    outmpileup = bam.with_suffix(".mpileup")
    cmd = f"samtools mpileup -f {reference} -s {bam} > {outmpileup}"
    run_subprocess(cmd)
    return outmpileup


def consensus_sequence(vcfgz: Path, reference: Path) -> Path:
    outconsensus = vcfgz.with_suffix(".consensus.fasta")
    logging.info("Generating consensus sequence with vcfutils.pl")
    cmd = f"cat {reference} | bcftools " \
          f"consensus {vcfgz} > {outconsensus}"
    run_subprocess(cmd)
    return outconsensus


def run_subprocess(cmd):
    p = Popen(cmd, shell=True)
    # p = Popen(cmd, shell=True, stdout=PIPE, stderr=PIPE)
    p.wait()


def parse_genome_results(mapping_stats: Path, sample_id=None) -> pd.DataFrame:
    # Read in stats .txt file
    df = pd.read_csv(mapping_stats, delimiter="\t", index_col=None)

    # Rename columns from bbmap output
    column_names = [
        "Reference_genome", "Avg_fold", "Length(bp)", "Ref_GC",
        "Covered_percent", "Covered_bases", "Plus_reads", "Minus_reads",
        "Read_GC", "Median_fold", "Std_dev"
    ]
    df.columns = column_names

    # Add SampleID
    if sample_id is None:
        df["SampleID"] = mapping_stats.name.split(".")[0]
    else:
        df["SampleID"] = sample_id

    # Rearrange column headers
    column_names_sorted = [
        "SampleID", "Reference_genome", "Avg_fold", "Median_fold",
        "Std_dev", "Length(bp)", "Ref_GC", "Covered_percent",
        "Covered_bases", "Plus_reads", "Minus_reads", "Read_GC"
    ]
    df = df[column_names_sorted]
    return df


def combine_dataframes(dfs: list) -> pd.DataFrame:
    df = pd.concat(dfs)
    return df


def get_readpair(sample_id: str, fastq_file_list: [Path], forward_id: str, reverse_id: str) -> list:
    """
    :param sample_id: String of sample ID
    :param fastq_file_list: List of fastq.gz file paths generated by retrieve_fastqgz()
    :param forward_id: ID indicating forward read in filename (e.g. _R1)
    :param reverse_id: ID indicating reverse read in filename (e.g. _R2)
    :return: the absolute filepaths of R1 and R2 for a given sample ID
    """

    r1, r2 = None, None
    for f in fastq_file_list:
        if sample_id in f.name:
            if forward_id in f.name:
                r1 = f
            elif reverse_id in f.name:
                r2 = f
    if r1 is not None:
        return [r1, r2]
    else:
        logging.debug('Could not pair {}'.format(sample_id))


if __name__ == "__main__":
    pipeline()
