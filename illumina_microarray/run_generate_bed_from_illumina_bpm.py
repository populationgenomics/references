#!/usr/bin/env python3

"""
Driver that submits a Hail Batch job to generate
GDA-8v1-0_D2-biallelic-snps.bed from the Illumina GDA BPM manifest, polarised
against the GRCh38 fasta. Both inputs already live under
gs://cpg-common-main/references/illumina_microarray/.

analysis-runner \
    --dataset common \
    --access-level full \
    --description 'Generate GDA biallelic SNV BED from BPM manifest' \
    --output-dir references/illumina_microarray \
    illumina_microarray/run_generate_bed_from_illumina_bpm.py

NOTES:
1. The output filename is fixed at GDA-8v1-0_D2-biallelic-snps.bed, matching
   the `GDA_8v1_0_D2_biallelic_snps_bed` entry under the `illumina_microarray`
   Source in references.py.
2. Idempotent: exits without submitting a batch if the BED is already present
   at the canonical references path.
3. Run at --access-level full --dataset common to land outputs at the
   canonical references location (gs://cpg-common-main/references/...).
   --access-level test writes to gs://cpg-common-test/references/... for
   staging/validation.
4. The worker (reference_generating_scripts/generate_bed_from_illumina_bpm.py)
   is shipped into the Batch job via a base64 heredoc, and its only deps
   (numpy, pyfaidx, IlluminaBeadArrayFiles) are pip-installed at job start in
   a plain python:3.11-slim image — no custom image needed.
"""

import base64
from pathlib import Path

from cloudpathlib import AnyPath
from cpg_utils.config import output_path, reference_path
from cpg_utils.hail_batch import get_batch, image_path


OUTPUT_FILENAME = 'GDA-8v1-0_D2-biallelic-snps.bed'
REPO_ROOT = Path(__file__).resolve().parent.parent
WORKER_SCRIPT = REPO_ROOT / 'reference_generating_scripts' / 'generate_bed_from_illumina_bpm.py'


def main() -> None:
    bpm_path = reference_path('illumina_microarray/GDA_8v1_0_D2_bpm')
    fasta_path = reference_path(
        'illumina_microarray/GCA_000001405_15_GRCh38_no_alt_analysis_set_fna',
    )
    fai_path = reference_path(
        'illumina_microarray/GCA_000001405_15_GRCh38_no_alt_analysis_set_fna_fai',
    )
    bed_out_path = output_path(OUTPUT_FILENAME)

    if AnyPath(bed_out_path).exists():
        print(f'{OUTPUT_FILENAME} already present at {bed_out_path}; nothing to do.')
        return

    worker_b64 = base64.b64encode(WORKER_SCRIPT.read_bytes()).decode()

    b = get_batch()
    bpm_input = b.read_input(bpm_path)
    # Dotted keys are supported: Batch localises these as {root}.fna and
    # {root}.fna.fai, which is what pyfaidx's default index lookup expects.
    fasta_group = b.read_input_group(**{'fna': fasta_path, 'fna.fai': fai_path})

    job = b.new_job(name='Generate GDA biallelic SNV BED from BPM')
    job.image(image_path('python'))
    job.storage('20Gi')
    job.cpu(2)
    job.memory('8Gi')

    job.command(
        f"""
        set -euo pipefail
        pip install --quiet --no-cache-dir \\
            numpy pyfaidx \\
            git+https://github.com/Illumina/BeadArrayFiles.git

        echo {worker_b64} | base64 -d > $BATCH_TMPDIR/generate_bed.py
        python3 $BATCH_TMPDIR/generate_bed.py \\
            --bpm {bpm_input} \\
            --fasta {fasta_group.fna} \\
            --output {job.ofile}
        """,
    )

    b.write_output(job.ofile, bed_out_path)
    b.run(wait=False)


if __name__ == '__main__':
    main()
