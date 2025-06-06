#!/usr/bin/env bash

#  copy filtered ancestry stage pruned_variants.ht frmo test to main
set -ex

gcloud storage cp -r \
gs://cpg-common-test/references/ancestry/gnomad_pruned_variants.ht \
gs://cpg-common-main/references/ancestry/gnomad_pruned_variants.ht

gcloud storage cp -r \
gs://cpg-common-test/references/ancestry/gnomad_pruned_variants.ht \
gs://cpg-common-main/references/ancestry/tenk10K_pruned_variants.ht