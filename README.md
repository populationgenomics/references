# Reference data for bioinformatics workflows

* Set prefix

```sh
export PREFIX=gs://cpg-reference
```

* `vep`: [preparing VEP data for Hail Query](vep/README.md)

* `broad`: sync the Broad's hg38 resources:

```sh
gsutil -u cpg-common rsync -r \
gs://gcp-public-data--broad-references/hg38/v0 \ 
$PREFIX/hg38/v0
```

* 'broad.sv': sync the GATK-SV resources

```sh
gsutil -u cpg-common rsync -r \
gs://gatk-sv-resources-public/hg38/v0/sv-resources \ 
$PREFIX/hg38/v0/sv-resources
```

* `gnomad`: sync resources from the GnomAD QC pipeline

```sh
gsutil -u cpg-common rsync -r \
gs://gnomad-public-requester-pays/resources/grch38
$PREFIX/gnomad/v0
```

* `seqr`: sync Seqr annotation resources

```sh
gsutil -u cpg-common rsync -r \
gs://seqr-reference-data/GRCh38/all_reference_data/combined_reference_data_grch38.ht
$PREFIX/seqr/combined_reference_data_grch38.ht

gsutil -u cpg-common rsync -r \
gs://seqr-reference-data/GRCh38/clinvar/clinvar.GRCh38.2022-09-17.ht
$PREFIX/seqr/clinvar.GRCh38.2022-09-17.ht
```

* `validation`: sync validation resources

```sh
gsutil -u cpg-common rsync -r \
gs://gcp-public-data--gnomad/resources/grch38/syndip/syndip.b38_20180222.bed \
$PREFIX/validation/syndip/syndip.b38_20180222.bed

gsutil -u cpg-common rsync -r \
gs://gcp-public-data--gnomad/resources/grch38/syndip/full.38.20180222.vcf.gz \
$PREFIX/validation/syndip/full.38.20180222.vcf.gz

gsutil -u cpg-common rsync -r \
gs://gcp-public-data--gnomad/resources/grch38/na12878/HG001_GRCh38_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_nosomaticdel_noCENorHET7.bed \
$PREFIX/validation/na12878/HG001_GRCh38_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_nosomaticdel_noCENorHET7.bed

gsutil -u cpg-common rsync -r \
gs://gcp-public-data--gnomad/resources/grch38/na12878/HG001_GRCh38_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_PGandRTGphasetransfer.vcf.gz \
$PREFIX/validation/na12878/HG001_GRCh38_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_PGandRTGphasetransfer.vcf.gz
```
