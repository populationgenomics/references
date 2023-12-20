#!/usr/bin/env bash

# check that the source location exists
echo "Checking gs://cpg-common-main-tmp/references/star/hg38"
if gcloud storage ls "gs://cpg-common-main-tmp/references/star/hg38" > /dev/null; then
    echo "Source location exists, continuing"
else
    echo "Source location vacant, exiting"
    exit 1
fi

# check that the destination directory doesn't already exist
echo "Checking gs://cpg-common-main/references/star/hg38"
if gcloud storage ls "gs://cpg-common-main/references/star/hg38" > /dev/null; then
    echo "Target location already exists, will not replace"
    exit 1
else
    echo "Target location vacant, copying"
    gsutil -m rsync -r "gs://cpg-common-main-tmp/references/star/hg38" "gs://cpg-common-main/references/star/hg38"
fi
