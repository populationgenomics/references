#!/usr/bin/env bash

# first identify the location of the MANE summary based on the version
MANE_VERSION=${1:-"1.4"}
MANE_URL="https://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/release_${MANE_VERSION}/MANE.GRCh38.v${MANE_VERSION}.summary.txt.gz"
LOCAL_SUMMARY_NAME="mane_{MANE_VERSION}.summary.txt.gz"
# download that file to a slightly simplified name
wget "${MANE_URL}" -O "${LOCAL_SUMMARY_NAME}"

# generate the JSON summary file
python3 reference_generating_scripts/mane_summary_parser.py \
    --input "${LOCAL_SUMMARY_NAME}" \
    --output "mane_{MANE_VERSION}.json" \
    --format json

# destination for the output
MAIN_OR_TEST=${2:-"main"}

# copy all the content to GCP - keep the same file names, copy in bulk
DESTINATION="gs://cpg-common-${MAIN_OR_TEST}/references/mane_${ENSEMBL_VERSION}"
DESTINATION_JSON="${DESTINATION}/mane_{MANE_VERSION}.json"
DESTINATION_SUMMARY="${DESTINATION}/${LOCAL_SUMMARY_NAME}"

gcloud storage cp "mane_{MANE_VERSION}.json" "${DESTINATION_JSON}"
gcloud storage cp "${LOCAL_SUMMARY_NAME}" "${DESTINATION_SUMMARY}"
