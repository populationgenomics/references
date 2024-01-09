# HGDP + 1KG reference sites table

This sites table is a representative subset of variants, generated from the publicly-available HGDP + 1KG dataset. Site selection is based off of the following criteria: sites are biallelic, autosomal, have an allele frequency above 1%, have a call rate above 99%, and have an inbreeding coefficient > -0.25 (this selection criteria is based predominantly off of the [gnomAD selection criteria](https://macarthurlab.org/2018/10/17/gnomad-v2-1/), with the exception that variants can be present in exomes or genomes).

This table can be used to generate a PCA for other highly divergent datasets, as tests conducted at the CPG showed little difference between sites selected from HGDP + 1KG compared to sites selected using HGDP + 1KG + other diverse cohorts.

The code used to generate this table can be found [here](https://github.com/populationgenomics/cpg-methods/tree/main/scripts/hgdp_1kg_marker_selection; see Hail batch [423327](https://batch.hail.populationgenomics.org.au/batches/423327/jobs/2) for log info).

The sites table can be found in `gs://cpg-common-main/references/ancestry/` and was copied from the `hgdp-1kg-main` bucket using the following command:

```bash
analysis-runner --dataset hgdp-1kg --access-level full --output-dir "references/ancestry" --description "copy sites table" sites_table/copy_sites_table.sh 
```

