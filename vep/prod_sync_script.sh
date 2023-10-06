#!/usr/bin/env bash

# a script for syncing tested VEP data to the production bucket
# takes a single argument for the VEP data to copy across
VERSION=$1

# check that the source location exists
echo "Checking gs://cpg-common-test/references/vep/${VERSION}/mount"
if gcloud storage ls "gs://cpg-common-test/references/vep/${VERSION}/mount" > /dev/null; then
    echo "Source location exists, continuing"
else
    echo "Source location vacant, exiting"
    exit 1
fi

# check that the destination directory doesn't already exist
echo "Checking gs://cpg-common-main/references/vep/${VERSION}/mount"
if gcloud storage ls "gs://cpg-common-main/references/vep/${VERSION}/mount" > /dev/null; then
    echo "Target location already exists, will not replace"
    exit 1
else
    echo "Target location vacant, copying"
    gcloud storage rsync --recursive "gs://cpg-common-test/references/vep/${VERSION}/mount" "gs://cpg-common-main/references/vep/${VERSION}/mount"
fi
