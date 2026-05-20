#!/usr/bin/env python3

"""
Promote the EGT-derived sample-less BCF (+ CSI index) from
gs://cpg-common-test/references/illumina_microarray/, where
extract_egt_info_bcf.sh --upload stages it, to the canonical references
location under gs://cpg-common-main/references/illumina_microarray/.

analysis-runner \\
    --dataset common \\
    --access-level full \\
    --description 'Promote EGT-derived BCF to cpg-common-main' \\
    --output-dir references/illumina_microarray \\
    illumina_microarray/copy_egt_info_bcf_to_main.py

NOTES:
1. Filenames are pinned to the entries declared under the `illumina_microarray`
   Source in references.py (GDA_8v1_0_D1_ClusterFile_egt_info_bcf{,_index}).
2. Idempotent: per-file existence check; objects already at the destination
   are skipped.
3. Source bucket is fixed at cpg-common-test (writable by analysts via the
   build-side extract_egt_info_bcf.sh --upload). Destination resolves via
   cpg_utils.config.output_path(), so --access-level full lands files at
   cpg-common-main; --access-level test is a no-op smoke test of the wiring
   (source == destination).
4. Runs entirely in the analysis-runner driver container — no Batch job is
   submitted.
"""

import subprocess

from cloudpathlib import AnyPath
from cpg_utils.config import output_path


SRC_PREFIX = 'gs://cpg-common-test/references/illumina_microarray'
FILENAMES = (
    'GDA-8v1-0_D1_ClusterFile_info.bcf',
    'GDA-8v1-0_D1_ClusterFile_info.bcf.csi',
)


def main() -> None:
    for name in FILENAMES:
        dst = output_path(name)
        if AnyPath(dst).exists():
            print(f'{name} already at {dst}; skipping.')
            continue
        src = f'{SRC_PREFIX}/{name}'
        print(f'Copying {src} -> {dst}')
        subprocess.check_call(['gcloud', 'storage', 'cp', src, dst])


if __name__ == '__main__':
    main()
