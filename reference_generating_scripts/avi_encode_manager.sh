#!/usr/bin/env bash

set -e

# This script is used to manage the generation of reference files for the AlphaGenome AVI SNV dataset.

# destination for the output
MAIN_OR_TEST=${1:-"main"}

WRITE_OUTPUT_TO="gs://cpg-common-${MAIN_OR_TEST}/references/alphagenome/avi_snv"

# localise the data
gcloud storage cp --do-not-decompress gs://cpg-common-main/references/alphagenome/avi_scores_snvs_tabix.zip ${BATCH_TMPDIR}/avi_scores_snvs_tabix.zip

# create unzip location, unzip into it, delete the original archive
mkdir ${BATCH_TMPDIR}/avi_data
unzip -d ${BATCH_TMPDIR}/avi_data ${BATCH_TMPDIR}/avi_scores_snvs_tabix.zip && rm ${BATCH_TMPDIR}/avi_scores_snvs_tabix.zip

# convert the TSV to a HT
python3 reference_generating_scripts/avi_snv_formatting.py \
    --avi_tsv ${BATCH_TMPDIR}/avi_data/alphagenome_variant_impact_score_snvs.tsv.gz \
    --ht_out ${BATCH_TMPDIR}/avi_snv.ht

# copy that up to GCP
TABLE_DESTINATION="${WRITE_OUTPUT_TO}/avi_snv.ht"
gcloud --quiet storage cp -r ${BATCH_TMPDIR}/avi_snv.ht "$TABLE_DESTINATION"
