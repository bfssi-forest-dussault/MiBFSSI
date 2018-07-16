#!/usr/bin/env python3

__version__ = "0.0.1"
__author__ = "Forest Dussault"
__email__ = "forest.dussault@canada.ca"

import os
import click
import logging
import multiprocessing

from pathlib import Path

from bin.tool_wrappers import call_qualimap, \
    call_snippy, \
    call_nullarbor, \
    call_bbmap, \
    taxid_reference_retrieval
from bin.accessories import get_sample_dictionary, \
    combine_dataframes, \
    parse_genome_results, \
    prepare_nullarbor_sample_file, \
    parse_snippy, \
    Sample

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
              required=False,
              default=None,
              help='Path to a reference .FASTA to use instead of the automatically acquired '
                   'references from sendsketch.sh. Must be specified if calling nullarbor.')
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
              help='Specify this flag to run Nullarbor against each sample. https://github.com/tseemann/nullarbor',
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
    if reference is None and nullarbor:
        logging.error("ERROR: Please provide a single reference with --reference in order to use nullarbor.")
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
        taxid_reference_retrieval(sample_list=sample_list,
                                  outdir=outdir)

    # Run analyses
    pipeline_flags = {"snippy": snippy, "bbmap": bbmap, "nullarbor": nullarbor}
    logging.info(f"Running {[x for x, y in pipeline_flags.items() if y is True]}")
    for sample in sample_list:
        try:
            # Generate BAM + sorted BAM with bbmap
            if bbmap:
                sample.bamfile, sample.mapping_stats = call_bbmap(fwd_reads=sample.r1,
                                                                  rev_reads=sample.r2,
                                                                  reference=sample.reference_genome,
                                                                  outdir=sample.outdir,
                                                                  sample_id=sample.sample_id,
                                                                  threads=threads)

                # Parse bbmap output and drop it into a dataframe
                sample.mapping_stats_df = parse_genome_results(mapping_stats=sample.mapping_stats,
                                                               sample_id=sample.sample_id)

                # Create final mapping stats report for every sample
                logging.info("Collating mapping statistics")
                mapping_stats_list = [sample.mapping_stats_df for sample in sample_list]
                mapping_stats_report = combine_dataframes(mapping_stats_list)
                mapping_stats_report.to_csv(outdir / "project_mapping_stats.csv", sep=",", index=None)

                # Qualimap
                sample.qualimap_dir = call_qualimap(bamfile=sample.bamfile,
                                                    outdir=sample.outdir,
                                                    threads=threads)

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
                sample.qualimap_dir = call_qualimap(bamfile=sample.bamfile,
                                                    outdir=sample.outdir,
                                                    threads=threads)

        except KeyError as e:
            logging.error(f"ERROR: Could not process {sample.sample_id}. Traceback:")
            logging.error(e)

    # Run nullarbor on all samples simultaneously
    if nullarbor:
        sample_file = prepare_nullarbor_sample_file(samples=sample_list,
                                                    outdir=outdir)
        call_nullarbor(project_name=outdir.name,
                       reference=reference,
                       outdir=outdir / 'nullarbor',
                       samples=sample_file,
                       threads=threads)


if __name__ == "__main__":
    pipeline()
