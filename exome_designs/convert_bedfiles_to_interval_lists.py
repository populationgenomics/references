#!/usr/bin/env python3

"""
This script converts a folder of bedfiles into picard interval_lists

analysis-runner \
    --dataset common \
    --access-level full \
    --description  'ref beds to interval_lists' \
    convert_bedfiles_to_interval_lists.py



NOTES:
1. Run in MAIN!
2. Assumes bedfiles are in TEST.
3. Defaults to references source exome_probesets and hg38: references/exome-probesets/hg38
4. Copies BED and INTERVAL_LIST to cpg-common-main
5. Takes a glob of the bed_source path: if re-run and previously converted BED files are still in test, they will be converted anew.

POSSIBLE TO DO:
arg for overwrite + check if interval_list already exists on cpg-common-main? see NOTES 5.

It is heavily indebted to a script which converts a single bed file:
https://github.com/populationgenomics/references/blob/main/reference_generating_scripts/convert_bed_to_interval_list.py
"""

import os
import sys
from collections.abc import Generator, Iterator
from pathlib import Path

import click
from cloudpathlib import AnyPath, CloudPath
from cpg_utils.config import ConfigError, config_retrieve, cpg_test_dataset_path, reference_path
from cpg_utils.hail_batch import get_batch, image_path

SOURCE = 'exome_probesets'


def get_bedfile_paths(exome_path: str) -> Iterator[Path] | Generator[CloudPath, None, None]:
    """
    Returns a generator of CloudPath for bed files found at "exome_path" in
    cpg-common-test/references.
    """
    
    bed_path = AnyPath(cpg_test_dataset_path(exome_path))
    bedfile_paths = bed_path.rglob('*.bed')
    return (bedfile_paths)


def make_interval_lists(bedfile_paths: Iterator[Path] | Generator[CloudPath, None, None], sd_ref: str , exome_ref_dict: dict[str, str]) -> None:
    """
    Converts a list of BED files in cpg-common-test to .interval_list files 
    using Picard BedToIntervalList.
    
    Saves the .interval_list file to the cpg-common-main reference directory.
    Copies the .bed file to the cpg-common-main reference directory.
    """
    sd_file = reference_path(sd_ref)
    b = get_batch()

    for bed_path in bedfile_paths:
        bed_file = bed_path.parts[-1]
        bed_stem = bed_path.stem

        picard_job = b.new_job(name=f'Convert {bed_file} to .interval_list file')
        picard_job.image(image_path('picard'))

        # lookup the bed and interval_list outpaths by matching their 'resource name' in references. will FAIL
        # if not in references.py
        bed_out_path = reference_path('exome_probesets/' + exome_ref_dict[bed_file])
        interval_list_out_path = reference_path('exome_probesets/' + exome_ref_dict[bed_stem + '.interval_list'])

        picard_job.command(
            f"""
            gcloud storage cp {str(bed_path)} {str(bed_out_path)}
            
            picard BedToIntervalList \
            -I {b.read_input(str(bed_path))} \
            -O {picard_job.ofile} \
            -SD {b.read_input(str(sd_file))}
            """,
        )

        b.write_output(picard_job.ofile, interval_list_out_path)

    b.run(wait=False)  


@click.command()
@click.option(
    '--exome-path',
    required=True,
    help='String identifying PATH (in cpg-common-test) of the bed files to be converted',
    default='references/exome-probesets/hg38',
)
@click.option(
    '--sd-ref',
    required=True,
    help='String identifier for the sequence dictionary file from the references.',
    default='broad/genome_calling_interval_lists',
)
def main(exome_path, sd_ref):
    """
    Converts a directory of BED files to .interval_list files 
    using Picard BedToIntervalList.
    
    Both the input BEDFILES and output INTERVAL_LIST files 
    must have been added to references.py.
    """

    bedfile_paths = get_bedfile_paths(exome_path)

    # make reverse map of actual_file_name -> reference
    try:
        ref_exome_dict = config_retrieve(['references', SOURCE])
    except ConfigError:
        sys.exit(f'No reference entry for provided source: {SOURCE}')
        
    exome_ref_dict = {os.path.basename(v): k for k, v in ref_exome_dict.items()}

    make_interval_lists(bedfile_paths, sd_ref, exome_ref_dict)


if __name__ == '__main__':
    main()
