#!/usr/bin/env bash

# This script is used to manage the generation of a BED file from the ensembl GFF3

# destination for the output

# first identify the location of the ensembl GFF3 based on its version
ENSEMBL_VERSION=${1:-"113"}
GFF3_URL="https://ftp.ensembl.org/pub/release-${ENSEMBL_VERSION}/gff3/homo_sapiens/Homo_sapiens.GRCh38.${ENSEMBL_VERSION}.gff3.gz"
LOCAL_OUTPUT_GFF3="GRCh38.gff3.gz"
LOCAL_OUTPUT_BED="GRCh38.bed"
MERGED_OUTPUT_BED="merged_GRCh38.bed"

# get that
wget "${GFF3_URL}" -O "${LOCAL_OUTPUT_GFF3}"

# parse that gff3 into a BED file
python3 reference_generating_scripts/generate_bed_from_ensembl.py \
    --gff3 "${LOCAL_OUTPUT_GFF3}" \
    --unmerged_output unmerged.bed \
    --merged_output merged.bed

## sort and merge the result - local script in lieu of using BedTools (not installed here)
## -k1,1V: Version sort to get chr ordering correct
## -k2,2n: numerical sort on col 2
sort -k1,1V -k2,2n unmerged.bed > "${LOCAL_OUTPUT_BED}"
sort -k1,1V -k2,2n merged.bed > "${MERGED_OUTPUT_BED}"

# copy all the content to GCP - keep the same file names, copy in bulk
WRITE_OUTPUT_TO="gs://cpg-common-main/references/ensembl_${ENSEMBL_VERSION}"
#GFF3_DESTINATION="${WRITE_OUTPUT_TO}/${LOCAL_OUTPUT_GFF3}"
#UNMERGED_BED_DESTINATION="${WRITE_OUTPUT_TO}/${LOCAL_OUTPUT_BED}"
#MERGED_BED_DESTINATION="${WRITE_OUTPUT_TO}/${MERGED_OUTPUT_BED}"

gcloud storage cp --do-not-decompress \
    "${LOCAL_OUTPUT_GFF3}" "${LOCAL_OUTPUT_BED}" "${MERGED_OUTPUT_BED}" \
    "${WRITE_OUTPUT_TO}"
