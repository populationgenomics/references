#!/usr/bin/env bash

# This script is used to manage the generation of reference files for the alpha missense dataset.

# destination for the output
MAIN_OR_TEST=${1:-"main"}

WRITE_OUTPUT_TO="gs://cpg-common-${MAIN_OR_TEST}/references/alphamissense"

# wget the raw data
AM_ZENODO="https://zenodo.org/records/8208688/files/AlphaMissense_hg38.tsv.gz"

# Download the raw data
wget -q -O alphamissense_38.tsv.gz $AM_ZENODO

# copy it up to GCP
TSV_DESTINATION="${WRITE_OUTPUT_TO}/alphamissense_38.tsv.gz"
gcloud storage cp --do-not-decompress alphamissense_38.tsv.gz "$TSV_DESTINATION"

# convert the TSV to a HT
python3 reference_generating_scripts/alphamissense_formatting.py \
    --am_tsv alphamissense_38.tsv.gz \
    --ht_out alphamissense_38.ht

# copy that up to GCP
TABLE_DESTINATION="${WRITE_OUTPUT_TO}/alphamissense_38.ht"
gcloud --quiet storage cp -r alphamissense_38.ht "$TABLE_DESTINATION"

# compress the output, write as a single file (easier to localise)
tar -czf alphamissense_38.ht.tar.gz alphamissense_38.ht

TAR_DESTINATION="${WRITE_OUTPUT_TO}/alphamissense_38.ht.tar.gz"
gcloud --quiet storage cp alphamissense_38.ht.tar.gz "$TAR_DESTINATION"
