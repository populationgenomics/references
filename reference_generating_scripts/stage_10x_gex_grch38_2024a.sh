#!/usr/bin/env bash
# Stage the FASTA, FASTA index and GTF of 10x Genomics' GRCh38-2024-A Cell
# Ranger reference in gs://cpg-common-test, from which the `tenx_gex_grch38_2024a`
# source in references.py is copied into cpg-common-main. The STAR index in the
# tarball is not staged: wf-single-cell reads only these three files.
#
# Requires: curl, tar, gcloud, ~16 GB of free disk under BATCH_TMPDIR (or TMPDIR),
# and write access to gs://cpg-common-test (test access on cpg-common).
#
# Usage: stage_10x_gex_grch38_2024a.sh [gs://destination/prefix]

set -euo pipefail

URL='https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCh38-2024-A.tar.gz'
# MD5 of the 11,481,868,421-byte tarball last modified 2024-03-08.
EXPECTED_MD5='a7b5b7ceefe10e435719edc1a8b8b2fa'
DEST="${1:-gs://cpg-common-test/references/10x/refdata-gex-GRCh38-2024-A}"
TOP='refdata-gex-GRCh38-2024-A'
FILES=(fasta/genome.fa fasta/genome.fa.fai genes/genes.gtf.gz)

# BATCH_TMPDIR is the storage-backed scratch on Hail Batch; elsewhere use TMPDIR.
WORK="$(mktemp -d "${BATCH_TMPDIR:-${TMPDIR:-/tmp}}/10x-gex.XXXXXX")"
trap 'rm -rf "${WORK}"' EXIT

curl --fail --location --silent --show-error --output "${WORK}/ref.tar.gz" "${URL}"

if command -v md5sum >/dev/null; then
    actual_md5="$(md5sum "${WORK}/ref.tar.gz" | cut -d' ' -f1)"
else
    actual_md5="$(md5 -q "${WORK}/ref.tar.gz")"
fi
if [[ "${actual_md5}" != "${EXPECTED_MD5}" ]]; then
    echo "MD5 mismatch: expected ${EXPECTED_MD5}, got ${actual_md5}" >&2
    exit 1
fi

tar -xzf "${WORK}/ref.tar.gz" -C "${WORK}" "${FILES[@]/#/${TOP}/}"

for f in "${FILES[@]}"; do
    gcloud storage cp "${WORK}/${TOP}/${f}" "${DEST%/}/${f}"
done
