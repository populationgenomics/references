# Reference data

Populating common reference resources for bioinformatics analysis. Assuming the Google Cloud infrastructure, but the structure can be replicated on other cloud provider.

Assuming the destination prefix and the project name are set:

```sh
export PREFIX=gs://cpg-reference
export PROJECT=cpg-common
```

## VEP

While VEP for Hail Query [is in progress](https://github.com/hail-is/hail/pull/12428), we provide a temporary solution [described here](vep/README.md). Hopefully to be deprecated soon.

## Broad

The Broad's common hg38 resources.

```sh
gsutil -u $PROJECT -m rsync -d -r \
gs://gcp-public-data--broad-references/hg38/v0 \ 
$PREFIX/hg38/v0
```

### GATK-SV

```sh
gsutil -u $PROJECT -m rsync -d -r \
gs://gatk-sv-resources-public/hg38/v0/sv-resources \ 
$PREFIX/hg38/v0/sv-resources
```

### GnomAD QC

```sh
gsutil -u $PROJECT -m rsync -d -r \
gs://gnomad-public-requester-pays/resources/grch38
$PREFIX/gnomad/v0
```

### Seqr annotation

```sh
gsutil -u $PROJECT -m rsync -d -r \
gs://seqr-reference-data/GRCh38/all_reference_data/combined_reference_data_grch38.ht
$PREFIX/seqr/combined_reference_data_grch38.ht

gsutil -u $PROJECT -m rsync -d -r \
gs://seqr-reference-data/GRCh38/clinvar/clinvar.GRCh38.2022-09-17.ht
$PREFIX/seqr/clinvar.GRCh38.2022-09-17.ht
```

### Variant calling validation

```sh
gsutil -u $PROJECT -m rsync -d -r \
gs://gcp-public-data--gnomad/resources/grch38/syndip $PREFIX/validation/syndip

gsutil -u $PROJECT -m rsync -d -r \
gs://gcp-public-data--gnomad/resources/grch38/na12878 $PREFIX/validation/na12878
```
