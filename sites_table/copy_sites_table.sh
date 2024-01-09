#!/usr/bin/env bash

set -ex

gsutil -m cp -r gs://cpg-hgdp-1kg-main/sites_table/v1-0/pruned_variants.ht gs://cpg-common-main/references/ancestry/
