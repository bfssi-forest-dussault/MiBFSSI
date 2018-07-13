import os
import shutil
import logging
import pandas as pd

from pathlib import Path
from subprocess import Popen
from bin.accessories import run_subprocess, run_subprocess_stdout


def call_blastn(database: Path, query_fasta: Path, outdir: Path) -> Path:
    # Index database
    cmd = f"makeblastdb -in {database} -parse_seqids -dbtype nucl"
    run_subprocess(cmd)

    outfile = outdir / query_fasta.with_suffix(".BLASTn").name
    cmd = f"blastn -query {query_fasta} -db {database} -out {outfile} " \
          f"-outfmt '6 qseqid stitle slen length qstart qend pident score'"
    run_subprocess(cmd)
    return outfile


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


def call_fuse(reference: Path) -> Path:
    fused_reference = reference.with_suffix(".fused.fna.gz")  # TODO: Make this less brittle
    cmd = f"fuse.sh in={reference} out={fused_reference} overwrite=true"
    run_subprocess(cmd)
    return fused_reference


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
    https://github.com/tseemann/snippy
    """
    outdir = outdir / "snippy"
    cmd = f"snippy --cpus {threads} --outdir {outdir} --ref {reference} --prefix {prefix} " \
          f"--R1 {fwd_reads} --R2 {rev_reads}"
    run_subprocess(cmd)
    return outdir


def call_nullarbor(project_name: str, reference: Path, samples: Path, outdir: Path, threads: int):
    """
    Let nullarbor+SKESA do all of the hard work. Currently not available for install via conda, but should be soon.
    Requires a minimum of 4 samples to run and uses a single reference genome for the whole batch for variant calling.

    https://github.com/tseemann/nullarbor
    """
    cmd = f"nullarbor.pl --name {project_name} --ref {reference} --input {samples} --outdir {outdir} " \
          f"--cpus {threads}"  # TODO: Add --trim parameter when it gets fixed in a new version of nullarbor
    out, err = run_subprocess_stdout(cmd)

    # Parse out nullarbor's instruction to run the 'nice make' command + run it
    for line in err.split("\n"):
        if "nice make" in line:
            nullarbor_cmd = line.split("] ")[1]
            print(nullarbor_cmd)
            p = Popen(['/bin/bash', '-c', nullarbor_cmd], cwd=str(Path.home()))
            p.wait()


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
