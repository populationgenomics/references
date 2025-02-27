#!/usr/bin/env bash

# This script is used to manage the generation of reference files for the alpha missense dataset.

WRITE_OUTPUT_TO="gs://cpg-common-test-tmp/references/alphamissense"

# wget the raw data
AM_ZENODO="https://zenodo.org/records/8208688/files/AlphaMissense_hg38.tsv.gz"

# Download the raw data
wget -O AlphaMissense_hg38.tsv.gz $AM_ZENODO

# copy it up to GCP
TSV_DESTINATION="${WRITE_OUTPUT_TO}/AlphaMissense_hg38.tsv.gz"
gcloud storage cp AlphaMissense_hg38.tsv.gz $TSV_DESTINATION

# convert the TSV to a HT
python3 reference_generating_scripts/alphamissense_formatting.py \
    --am_tsv AlphaMissense_hg38.tsv.gz \
    --ht_out output.ht

# copy that up to GCP
TABLE_DESTINATION="${WRITE_OUTPUT_TO}/AlphaMissense_hg38.ht"
gcloud storage cp -r output.ht $TABLE_DESTINATION
