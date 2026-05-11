#!/usr/bin/env python3

"""
This script lifts over hg19 exome BED files to hg38 with picard
LiftOverIntervalList, and copies the resulting BED and interval_list files
to cpg-common-main references.

analysis-runner \
    --dataset common \
    --access-level full \
    --description  'liftover hg19 exome beds to hg38' \
    --output-dir cpg-common-main/references \
    liftover_and_convert_hg19_bedfiles.py


NOTES:
1. Run in MAIN!
2. Assumes hg19 bedfiles are in TEST, in a dedicated hg19/ subfolder
   alongside the existing hg38/ layout:
   gs://cpg-common-test/references/exome-probesets/hg19/*_hg19.bed
3. Requires hg19_dict to be present in references (run generate_hg19_dict.py
   once before this script).
4. Per-file pipeline:
     picard BedToIntervalList    (hg19 BED -> hg19 interval_list, SD=hg19_dict)
     picard LiftOverIntervalList (hg19 -> hg38 interval_list, target SD=broad)
     picard IntervalListToBed    (hg38 interval_list -> hg38 BED)
5. Both the BED and interval_list output filenames (with _hg19 substituted for
   _hg38) must exist in references.py under the exome_probesets Source.
6. Logs the count of rejected intervals from picard LiftOverIntervalList so
   the operator can sanity-check liftover loss before merging into ICA.
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


def get_hg19_bedfile_paths(exome_path: str) -> Iterator[Path] | Generator[CloudPath, None, None]:
    """
    Returns a generator of CloudPath for hg19 bed files (``*_hg19.bed``)
    found at ``exome_path`` in cpg-common-test/references.
    """

    bed_path = AnyPath(cpg_test_dataset_path(exome_path))
    return bed_path.rglob('*_hg19.bed')


def liftover_hg19_bedfiles(
    bedfile_paths: Iterator[Path] | Generator[CloudPath, None, None],
    hg38_sd_ref: str,
    chain_ref: str,
    hg19_dict_ref: str,
    exome_ref_dict: dict[str, str],
) -> None:
    """
    For each ``*_hg19.bed`` in cpg-common-test, run a single picard job that
    chains BedToIntervalList -> LiftOverIntervalList -> IntervalListToBed and
    writes the resulting hg38 BED and interval_list to cpg-common-main at the
    paths declared by references.py.
    """

    b = get_batch()
    hg38_sd_input = b.read_input(reference_path(hg38_sd_ref))
    hg19_dict_input = b.read_input(reference_path(hg19_dict_ref))
    chain_input = b.read_input(reference_path(chain_ref))

    for bed_path in bedfile_paths:
        hg19_bed_file = bed_path.parts[-1]
        hg38_bed_file = hg19_bed_file.replace('_hg19.bed', '_hg38.bed')
        hg38_interval_list_file = hg19_bed_file.replace('_hg19.bed', '_hg38.interval_list')

        try:
            bed_out_key = exome_ref_dict[hg38_bed_file]
            interval_out_key = exome_ref_dict[hg38_interval_list_file]
        except KeyError as missing:
            sys.exit(
                f'No reference entry for expected hg38 output {missing.args[0]} '
                f'(derived from {hg19_bed_file}). Check references.py.',
            )

        bed_out_path = reference_path(f'{SOURCE}/{bed_out_key}')
        interval_list_out_path = reference_path(f'{SOURCE}/{interval_out_key}')

        picard_job = b.new_job(name=f'Liftover {hg19_bed_file} hg19 -> hg38')
        picard_job.image(image_path('picard'))
        picard_job.declare_resource_group(
            hg38_outputs={
                'bed': '{root}.bed',
                'interval_list': '{root}.interval_list',
            },
        )

        picard_job.command(
            f"""
            set -euo pipefail

            picard BedToIntervalList \
                I={b.read_input(str(bed_path))} \
                O=$BATCH_TMPDIR/hg19.interval_list \
                SD={hg19_dict_input}

            picard LiftOverIntervalList \
                I=$BATCH_TMPDIR/hg19.interval_list \
                O=$BATCH_TMPDIR/hg38.interval_list \
                SD={hg38_sd_input} \
                CHAIN={chain_input} \
                REJECT=$BATCH_TMPDIR/rejected.interval_list

            rejected=$(grep -cv '^@' $BATCH_TMPDIR/rejected.interval_list || true)
            echo "Liftover rejected intervals for {hg19_bed_file}: $rejected"

            picard IntervalListToBed \
                I=$BATCH_TMPDIR/hg38.interval_list \
                O=$BATCH_TMPDIR/hg38.bed

            cp $BATCH_TMPDIR/hg38.interval_list {picard_job.hg38_outputs.interval_list}
            cp $BATCH_TMPDIR/hg38.bed {picard_job.hg38_outputs.bed}
            """,
        )

        b.write_output(picard_job.hg38_outputs.bed, bed_out_path)
        b.write_output(picard_job.hg38_outputs.interval_list, interval_list_out_path)

    b.run(wait=False)


@click.command()
@click.option(
    '--exome-path',
    required=True,
    help='String identifying PATH (in cpg-common-test) of the hg19 bed files.',
    default='references/exome-probesets/hg19',
)
@click.option(
    '--hg38-sd-ref',
    required=True,
    help='Reference key for the hg38 sequence dictionary used as liftover target SD.',
    default='broad/genome_calling_interval_lists',
)
@click.option(
    '--chain-ref',
    required=True,
    help='Reference key for the hg19 -> hg38 liftover chain.',
    default='liftover_37_to_38',
)
@click.option(
    '--hg19-dict-ref',
    required=True,
    help='Reference key for the hg19 picard sequence dictionary.',
    default='hg19_dict',
)
def main(exome_path, hg38_sd_ref, chain_ref, hg19_dict_ref):
    """
    Lift over hg19 exome BED files in cpg-common-test to hg38 with picard,
    and write hg38 BED + interval_list outputs to the paths declared in
    references.py under exome_probesets.
    """

    bedfile_paths = get_hg19_bedfile_paths(exome_path)

    # reverse map of actual_file_name -> reference subkey
    try:
        ref_exome_dict = config_retrieve(['references', SOURCE])
    except ConfigError:
        sys.exit(f'No reference entry for provided source: {SOURCE}')

    exome_ref_dict = {os.path.basename(v): k for k, v in ref_exome_dict.items()}

    liftover_hg19_bedfiles(
        bedfile_paths, hg38_sd_ref, chain_ref, hg19_dict_ref, exome_ref_dict,
    )


if __name__ == '__main__':
    main()
