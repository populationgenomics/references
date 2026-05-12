#!/usr/bin/env python3

"""
One-shot bootstrap that generates hg19.dict alongside the hg19 fasta in
gs://cpg-common-main/references/hg19/v0/. The dict is consumed by
LiftOverIntervalList inside liftover_and_convert_hg19_bedfiles.py.

analysis-runner \
    --dataset common \
    --access-level full \
    --description  'Bootstrap hg19 picard sequence dictionary' \
    --output-dir cpg-common-main/references \
    generate_hg19_dict.py

NOTES:
1. Run in MAIN (hg19.dict lives in cpg-common-main).
2. Requires hg19.fa.gz to have been transferred by CI via the `hg19_fasta`
   Source in references.py.
3. Idempotent: exits without submitting a batch if hg19.dict is already
   present at reference_path('hg19_dict').
"""

import sys

from cloudpathlib import AnyPath
from cpg_utils.config import reference_path
from cpg_utils.hail_batch import get_batch, image_path


def main() -> None:
    """
    Submit a single picard CreateSequenceDictionary job to produce hg19.dict
    next to hg19.fa.gz in references.
    """

    fasta_path = reference_path('hg19_fasta')
    dict_path = reference_path('hg19_dict')

    if AnyPath(dict_path).exists():
        print(f'hg19 dict already present at {dict_path}; nothing to do.')
        return

    if not AnyPath(fasta_path).exists():
        sys.exit(
            f'hg19 fasta not found at {fasta_path}. Run CI on the references '
            f'repo to transfer hg19.fa.gz before generating the dict.',
        )

    b = get_batch()
    fasta_input = b.read_input(fasta_path)

    picard_job = b.new_job(name='Generate hg19 sequence dictionary')
    picard_job.image(image_path('picard'))
    picard_job.storage('20Gi')

    picard_job.command(
        f"""
        gunzip -c {fasta_input} > $BATCH_TMPDIR/hg19.fa
        picard CreateSequenceDictionary \
            R=$BATCH_TMPDIR/hg19.fa \
            O={picard_job.ofile}
        """,
    )

    b.write_output(picard_job.ofile, dict_path)
    b.run(wait=False)


if __name__ == '__main__':
    main()
