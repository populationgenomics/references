#!/usr/bin/env bash

# Promote EGT-derived sample-less BCF + CSI from cpg-common-test to cpg-common-main.
# Run via analysis-runner at --access-level full --dataset common.

set -ex

gcloud storage cp \
    gs://cpg-common-test/references/illumina_microarray/GDA-8v1-0_D1_ClusterFile_info.bcf \
    gs://cpg-common-main/references/illumina_microarray/GDA-8v1-0_D1_ClusterFile_info.bcf

gcloud storage cp \
    gs://cpg-common-test/references/illumina_microarray/GDA-8v1-0_D1_ClusterFile_info.bcf.csi \
    gs://cpg-common-main/references/illumina_microarray/GDA-8v1-0_D1_ClusterFile_info.bcf.csi
