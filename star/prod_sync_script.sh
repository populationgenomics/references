#!/usr/bin/env bash

# a script for syncing tested STAR reference to the production bucket
# takes a single argument for the STAR reference data to copy across
STAR_VERSION=$1

# check that the source location exists
echo "Checking gs://cpg-common-main-tmp/references/star/${STAR_VERSION}/hg38"
if gcloud storage ls "gs://cpg-common-main-tmp/references/star/${STAR_VERSION}/hg38" > /dev/null; then
    echo "Source location exists, continuing"
else
    echo "Source location vacant, exiting"
    exit 1
fi

# check that the destination directory doesn't already exist
echo "Checking gs://cpg-common-main/references/star/${STAR_VERSION}/hg38"
if gcloud storage ls "gs://cpg-common-main/references/star/${STAR_VERSION}/hg38" > /dev/null; then
    echo "Target location already exists, will not replace"
    exit 1
else
    echo "Target location vacant, copying"
    gsutil rsync -r -m "gs://cpg-common-main-tmp/references/star/${STAR_VERSION}/hg38" "gs://cpg-common-main/references/star/${STAR_VERSION}/hg38"
fi
