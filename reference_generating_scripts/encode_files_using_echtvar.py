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
import json
import random
from argparse import ArgumentParser
from pathlib import Path
from string import ascii_lowercase

from os.path import join

from cpg_utils import to_path
from cpg_utils.config import config_retrieve, image_path, ConfigError
from cpg_utils.hail_batch import get_batch

CANONICAL_CHROMOSOMES = [f'chr{x}' for x in list(range(1, 23)) + ['X', 'Y']] + [
    'whole_genome'
]


# pull images from config, defaulting to the images at time of writing
try:
    echtvar_image = image_path('echtvar')
    bcftools_image = image_path('bcftools_120')
except ConfigError:
    echtvar_image = 'australia-southeast1-docker.pkg.dev/cpg-common/images/echtvar:v0.2.1'
    bcftools_image = 'australia-southeast1-docker.pkg.dev/cpg-common/images/bcftools_120:1.20'


def storage_with_buffer(file_path: str, buffer: int = 10) -> int:
    """
    determine the storage requirement for a file, adding a buffer
    Args:
        file_path ():
        buffer (int): number of GiB to add to the storage requirement
    """
    # determine exact storage requirement, add a buffer for safety and outputs
    return (to_path(file_path).stat().st_size // 1024**3) + buffer


def encode_gnomad(region: str | None = None) -> None:
    """
    run echtvar encode on all gnomadV4 contigs, separately and combined
    we need to do this once ever, estimated cost $5

    This includes the potential for region filtering, which will reduce the size of the output files and processing time
    Useful in Talos (which is gene-centric), but not relevant for the core pipeline which requires whole-genome AFs

    Args:
        region (str): optional, path to a BED file containing a subset of regions to encode

    Returns:
        None, writes to storage directly
    """

    common_folder = join(
        config_retrieve(['storage', 'common', 'default']),
        'gnomad',
        'echtvar',
    )

    # set the filename template to use for this run
    if region is not None:
        logging.info(f'Running echtvar on gnomad v4.1, region subset: {region}')

        # include the name (no extensions) of the region file in the output name
        region_name = Path(region).name.split('.')[0]

        # double layered templating - 'chrom' will be inserted later
        output_template = join(common_folder, f'gnomad_4.1_region_{region_name}_{{chrom}}.zip')

    else:
        output_template = join(common_folder, 'gnomad_4.1_{chrom}.zip')

    # get the name of the final output - if that already exists, this becomes much less work
    wg_output = output_template.format(chrom='whole_genome')
    wg_exists: bool = to_path(wg_output).exists()

    contig_files = []
    storage_running_total = 0
    for contig in CANONICAL_CHROMOSOMES:

        # don't do this for the whole genome output
        if contig == 'whole_genome':
            continue

        # output path for this contig
        contig_output = output_template.format(chrom=contig)
        contig_output_exists = to_path(contig_output).exists()

        if contig_output_exists and wg_exists:
            logging.info(f'Skipping echtvar on {contig}, output already exists. No need to retain for genome-wide')
            continue

        # select the file, read in, and determine size
        file_path = config_retrieve(['references', 'gnomad_4.1_vcfs', contig])

        contig_vcf = get_batch().read_input_group(vcf=file_path, index=f'{file_path}.tbi')['vcf']

        job_storage = storage_with_buffer(file_path)

        localised_region = get_batch().read_input(region)

        # if we only want to run on a subset of the genome, read in the BED file
        if region is not None:
            trim_job = get_batch().new_bash_job(f'Trim {contig} to specified region')
            trim_job.image(bcftools_image)
            trim_job.storage(f'{job_storage}Gi')
            trim_job.cpu(4)
            trim_job.command(
                f'bcftools view '
                f'-R {localised_region} '
                f'--regions-overlap 2 '
                f'{contig_vcf} '
                f'-Oz -o {trim_job.output}',
            )
            # update this value to the trimmed file
            contig_vcf = trim_job.output

        # if we need this for the whole-genome encoding, add to the list of inputs
        if not wg_exists:
            contig_files.append(contig_vcf)

        # if we don't need to generate this contig's file, continue
        if contig_output_exists:
            continue

        # create and resource a job
        contig_job = get_batch().new_job(f'Run echtvar on gnomad v4.1, {contig}, Region: {region or "unrestricted"}')
        contig_job.image(echtvar_image)
        contig_job.storage(f'{job_storage}Gi')
        contig_job.cpu(4)
        contig_job.memory('highmem')

        # run the echtvar encode command
        contig_job.command(f'echtvar encode {contig_job.output} $ECHTVAR_CONFIG {contig_vcf}')
        get_batch().write_output(contig_job.output, contig_output)

        # add to the total storage required for the whole genome job. Plan for worst case (region == whole genome)
        storage_running_total += job_storage

    # finally, take all the contig files (region filtered, or not), and run echtvar (unless genome-wide already exists)
    # this job becomes implicitly dependent on any previous region-filtering jobs from use of prior output as input
    if not wg_exists:
        logging.info(f'Running echtvar on whole genome, Region: {region or "unrestricted"}')
        job = get_batch().new_job(f'Run echtvar on gnomad v4.1, whole genome, Region: {region or "unrestricted"}')
        job.image(echtvar_image)
        job.storage(f'{storage_running_total}Gi')
        job.cpu(4)
        job.memory('highmem')
        # the input files were all localised individually
        job.command(
            f'echtvar encode {job.output} $ECHTVAR_CONFIG {" ".join(contig_files)}'
        )
        get_batch().write_output(job.output, wg_output)

    get_batch().run(wait=False)


def encode_anything(input_list: list[str], output: str, region: str | None = None):
    """
    Takes a list of VCFs, and runs echtvar on them to encode the results into a minimised representation
    Echtvar requires a config, so we'll pull one from the run config and log the contents

    Writes the output to a single zip file (e.g. for per-chromsome encoding, submit multiple jobs)

    Args:
        input_list ():
        output ():
        region (str): optional, path to a BED file containing a subset of regions to encode
    """

    # pull the echtvar config from the run config
    # if successful write it to a temp file, then read into the batch
    # the intention here is to make the config used more flexible, as it has to be a localised JSON in the container
    # if nothing was supplied, we'll use the default config (gnomAD)
    try:
        echtvar_config = config_retrieve(['echtvar_config'])
        random_string = '_'.join(random.choices(ascii_lowercase, k=20)) + '.json'
        temp_config = join(
            config_retrieve(['storage', 'default', 'tmp']), random_string
        )
        with to_path(temp_config).open('w') as handle:
            json.dump(echtvar_config, handle, indent=2)
        config = get_batch().read_input(temp_config)

    except ConfigError:
        logging.info('No echtvar config entry found, using default (set up for gnomAD)')
        config = '$ECHTVAR_CONFIG'

    total_storage = 0
    localised_inputs = []
    for input_file in input_list:
        single_localised_input = get_batch().read_input(input_file)
        input_size = storage_with_buffer(input_file)
        if region is not None:
            localised_region = get_batch().read_input(region)
            # if we only want to run on a subset of the genome, read in the BED file
            trim_job = get_batch().new_bash_job(f'Trim {input_file} to specified region')
            trim_job.image(bcftools_image)
            trim_job.storage(f'{input_size}Gi')
            trim_job.cpu(4)
            trim_job.command(
                f'bcftools view '
                f'-R {localised_region} '
                f'--regions-overlap 2 '
                f'{single_localised_input} '
                f'-Oz -o {trim_job.output}',
            )
            localised_inputs.append(trim_job.output)
        else:
            localised_inputs.append(single_localised_input)
        total_storage += input_size

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
    parser.add_argument(
        '--input',
        help='Path to input data. If not supplied will default to gnomad_4.1 vcfs',
        nargs='+',
        default=[],
    )
    parser.add_argument(
        '--output',
        help='Path to write the result - all arguments supplied with --input will be processed together',
    )
    parser.add_argument(
        '--region',
        help='Path to a BED file, used to limit encoded data to a specific region',
        required=False,
        default=None,
    )
    args = parser.parse_args()

    if len(args.input) == 0:
        encode_gnomad(region=args.region)
    else:
        encode_anything(input_list=args.input, output=args.output, region=args.region)
