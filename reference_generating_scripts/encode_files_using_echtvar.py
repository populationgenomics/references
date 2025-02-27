"""
all processes related to echtvar - generating a VCF annotation resource from raw inputs, or applying those annotations
this is not confined to a specific cohort/project

This script is intended to be run using analysis-runner, generating one or more jobs.

The default behaviour is to run for gnomAD, writing to cpg-common
Alternatively, any number of input VCFs can be provided, and the script will run echtvar on them, creating a single file

Custom configs can be provided by adding a config to the echtvar_config section of the config file
See https://github.com/brentp/echtvar?tab=readme-ov-file#configuration-file-for-encode

The default config in the echtvar config is:
[
    {"field": "AC", "alias": "gnomad_AC", "missing_value": -2147483648},
    {"field": "AN", "alias": "gnomad_AN", "missing_value": -2147483648},
    {"field": "AF", "alias": "gnomad_AF", "missing_value": -0, "multiplier": 2000000},
    {"field": "nhomalt", "alias": "gnomad_HomAlt", "missing_value": -2147483648},
]

Which would be encoded in toml as:
[[echtvar_config]]
field = "AC"
alias = "gnomad_AC"
missing_value = -2147483648

[[echtvar_config]]
field = "AN"
alias = "gnomad_AN"
missing_value = -2147483648

[[echtvar_config]]
field = "AF"
alias = "gnomad_AF"
missing_value = 0
multiplier = 2000000

[[echtvar_config]]
field = "nhomalt"
alias = "gnomad_HomAlt"
missing_value = -2147483648

It may be easier to generate the config file you want, then run it through toml.dumps() to get the toml string.
"""

import logging
import random
from argparse import ArgumentParser

from os.path import join

from cpg_utils import Path, to_path
from cpg_utils.config import config_retrieve, image_path, ConfigError
from cpg_utils.hail_batch import get_batch

CANONICAL_CHROMOSOMES = [f'chr{x}' for x in list(range(1, 23)) + ['X', 'Y']] + ['whole_genome']


# pull from config, defaulting to the image at time of writing
try:
    echtvar_image = image_path('echtvar')
except ConfigError:
    echtvar_image = 'australia-southeast1-docker.pkg.dev/cpg-common/images/echtvar:v0.2.1'

word_file = "/usr/share/dict/words"
WORDS = open(word_file).read().splitlines()


def storage_with_buffer(file_path: str, buffer: int = 10) -> int:
    """
    determine the storage requirement for a file, adding a buffer
    Args:
        file_path ():
        buffer (int): number of GiB to add to the storage requirement
    """
    # determine exact storage requirement, add a buffer for safety and outputs
    return (to_path(file_path).stat().st_size // 1024**3) + buffer


def encode_gnomad() -> StageOutput | None:
    """
    run echtvar encode on all gnomadV4 contigs, separately and combined
    we need to do this once ever, estimated cost $5

    This job only needs to run once, ever
    """

    common_folder = join(config_retrieve(['storage', 'common']), 'gnomad', 'echtvar')
    output_template = join(common_folder, 'gnomad_4.1_{chrom}.zip')

    contig_files = []
    storage_running_total = 0
    for contig in CANONICAL_CHROMOSOMES:

        contig_output = output_template.format(chrom=contig)
        if to_path(contig_output).exists():
            logging.info(f'Skipping echtvar on {contig}, output already exists')
            continue

        # don't do this for the whole genome output
        if contig == 'whole_genome':
            continue

        # localise this one file
        file_path = config_retrieve(['references', 'gnomad_4.1_vcfs', contig])
        contig_localised = get_batch().read_input(file_path)
        # add to the list of inputs for the whole genome job
        contig_files.append(contig_localised)

        # create and resource a job
        contig_job = get_batch().new_job(f'Run echtvar on gnomad v4.1, {contig}')
        contig_job.image(echtvar_image)
        job_storage = storage_with_buffer(contig_job.output)
        contig_job.storage(f'{job_storage}Gi')
        contig_job.cpu(4)
        contig_job.memory('highmem')

        # run the echtvar encode command
        contig_job.command(f'echtvar encode {contig_job.output} $ECHTVAR_CONFIG {contig_localised}')
        get_batch().write_output(contig_job.output, contig_output)

        # add to the total storage required for the whole genome job
        storage_running_total += job_storage

    whole_genome_output = output_template.format(chrom='whole_genome')
    if not to_path(whole_genome_output).exists():
        logging.info('Running echtvar on whole genome')
        job = get_batch().new_job('Run echtvar on gnomad v4.1, whole genome')
        job.image(echtvar_image)
        job.storage(f'{storage_running_total}Gi')
        job.cpu(4)
        job.memory('highmem')
        # the input files were all localised individually
        job.command(f'echtvar encode {job.output} $ECHTVAR_CONFIG {" ".join(contig_files)}')
        get_batch().write_output(job.output, whole_genome_output)

    get_batch().run(wait=False)


def encode_anything(input_list: list[str], output: str):
    """
    Takes a list of VCFs, and runs echtvar on them to encode the results into a minimised representation
    Echtvar requires a config, so we'll pull one from the run config and log the contents

    Writes the output to a single zip file (e.g. for per-chromsome encoding, submit multiple jobs)

    Args:
        input_list ():
        output ():
    """

    # pull the echtvar config from the run config
    # if successful write it to a temp file, then read into the batch
    # the intention here is to make the config used more flexible, as it has to be a localised JSON in the container
    # if nothing was supplied, we'll use the default config (gnomAD)
    try:
        echtvar_config = config_retrieve(['echtvar_config'])
        random_string = '_'.join(random.choices(WORDS, k=3)) + '.json'
        temp_config = join(config_retrieve(['storage', 'default', 'tmp']), random_string)
        with to_path(temp_config).open('w') as handle:
            json.dump(echtvar_config, handle, indent=2)
        config = get_batch().read_input(temp_config)

    except ConfigError:
        logging.info('No echtvar config entry found, using default (set up for gnomAD)')
        config = '$ECHTVAR_CONFIG'

    total_storage = 0
    localised_inputs = []
    for input in input_list:
        localised_inputs.append(get_batch().read_input(input))
        total_storage += storage_with_buffer(input)

    # create a job to run echtvar on all the input VCFs
    job = get_batch().new_job('Run echtvar on all input VCFs')
    job.image(echtvar_image)
    job.storage(f'{total_storage}Gi')
    job.cpu(4)
    job.memory('highmem')
    job.command(f'echtvar encode {job.output} {config} {" ".join(localised_inputs)}')
    get_batch().write_output(job.output, output)
    get_batch().run(wait=False)


if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO)
    parser = ArgumentParser()
    parser.add_argument('--input', help='Path to input data. If not supplied will default to gnomad_4.1 vcfs', nargs='+', default=[],)
    parser.add_argument('--output', help='Path to write the result - all arguments supplied with --input will be processed together')
    args = parser.parse_args()

    if len(args.input) == 0:
        encode_gnomad()
    else:
        encode_anything(input_list=args.input, output=args.output)
