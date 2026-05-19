#!/usr/bin/env python3

"""
Driver that generates GDA-8v1-0_D2-biallelic-snps.bed from the Illumina GDA BPM
manifest, polarised against the GRCh38 fasta. Both inputs already live under
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
2. Idempotent: exits without doing any work if the BED is already present at
   the canonical references path.
3. Run at --access-level full --dataset common to land outputs at the
   canonical references location (gs://cpg-common-main/references/...).
   --access-level test writes to gs://cpg-common-test/references/... for
   staging/validation.
4. Runs entirely in the analysis-runner driver container — no Batch job is
   submitted. pyfaidx and IlluminaBeadArrayFiles are pip-installed at startup;
   inputs are copied in with gcloud storage; the worker is imported from the
   same repo and called in-process; output is copied back out with gcloud.
"""

import subprocess
import sys
import tempfile
from pathlib import Path

from cloudpathlib import AnyPath
from cpg_utils.config import output_path, reference_path


OUTPUT_FILENAME = 'GDA-8v1-0_D2-biallelic-snps.bed'
BEADARRAYFILES_TARBALL = (
    'https://github.com/Illumina/BeadArrayFiles/archive/refs/tags/1.3.4.tar.gz'
)
REPO_ROOT = Path(__file__).resolve().parent.parent


def main() -> None:
    bed_out_path = output_path(OUTPUT_FILENAME)
    if AnyPath(bed_out_path).exists():
        print(f'{OUTPUT_FILENAME} already present at {bed_out_path}; nothing to do.')
        return

    bpm_src = reference_path('illumina_microarray/GDA_8v1_0_D2_bpm')
    fasta_src = reference_path(
        'illumina_microarray/GCA_000001405_15_GRCh38_no_alt_analysis_set_fna',
    )
    fai_src = reference_path(
        'illumina_microarray/GCA_000001405_15_GRCh38_no_alt_analysis_set_fna_fai',
    )

    subprocess.check_call([
        sys.executable, '-m', 'pip', 'install', '--quiet', '--no-cache-dir',
        'numpy', 'pyfaidx', BEADARRAYFILES_TARBALL,
    ])

    with tempfile.TemporaryDirectory() as tmp:
        tmpdir = Path(tmp)
        bpm = tmpdir / 'manifest.bpm'
        fasta = tmpdir / 'ref.fna'
        fai = tmpdir / 'ref.fna.fai'
        bed = tmpdir / OUTPUT_FILENAME

        for src, dst in ((bpm_src, bpm), (fasta_src, fasta), (fai_src, fai)):
            subprocess.check_call(['gcloud', 'storage', 'cp', src, str(dst)])

        sys.path.insert(0, str(REPO_ROOT))
        from reference_generating_scripts.generate_bed_from_illumina_bpm import (
            main as generate_bed,
        )
        generate_bed(str(bpm), str(fasta), str(bed))

        subprocess.check_call(['gcloud', 'storage', 'cp', str(bed), bed_out_path])


if __name__ == '__main__':
    main()
