#!/usr/bin/env bash

set -ex
curl -L -A "Mozilla/5.0" "https://s3.amazonaws.com/igv.org.genomes/hg38/annotations/cytoBandIdeo.txt.gz" | gcloud --billing-project cpg-common storage cp - "gs://cpg-common-main/references/igv_org_genomes/hg38/annotations/cytoBandIdeo.txt.gz"
curl -L -A "Mozilla/5.0" "https://s3.amazonaws.com/igv.org.genomes/hg38/hg38_alias.tab" | gcloud --billing-project cpg-common storage cp - "gs://cpg-common-main/references/igv_org_genomes/hg38/hg38_alias.tab"
curl -L -A "Mozilla/5.0" "https://s3.amazonaws.com/igv.org.genomes/hg38/refGene.sorted.txt.gz" | gcloud --billing-project cpg-common storage cp - "gs://cpg-common-main/references/igv_org_genomes/hg38/refGene.sorted.txt.gz"
