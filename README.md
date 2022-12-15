# Reference data

Populating common reference resources for bioinformatics analysis. Assuming the Google Cloud infrastructure, but the structure can be replicated for other cloud providers.

Script `references.py` describes sources of reference resources. Each object of class `Source` specifies one source location to pull data from (either a GCS bucket or an HTTP URL), along with an optional map of keys/files to expand this source as a section in the finalised config. E.g.

```py
Source(
    'liftover_38_to_37',
    src='gs://hail-common/references/grch38_to_grch37.over.chain.gz',
    dst='liftover/grch38_to_grch37.over.chain.gz',
)
```

Is expanded into flat

```toml
liftover_38_to_37 = "gs://cpg-common-main/references/liftover/grch38_to_grch37.over.chain.gz"
```

Whereas

```py
Source(
    'gatk_sv',
    src='gs://gatk-sv-resources-public/hg38/v0/sv-resources',
    dst='hg38/v0/sv-resources',
    files=dict(
        wham_include_list_bed_file='resources/v1/wham_whitelist.bed',
        primary_contigs_list='resources/v1/primary_contigs.list',
    )
)
```

Is expanded into a section

```toml
[gatk_sv]
wham_include_list_bed_file = "gs://cpg-common-main/references/hg38/v0/sv-resources/resources/v1/wham_whitelist.bed"
primary_contigs_list = "gs://cpg-common-main/references/hg38/v0/sv-resources/resources/v1/primary_contigs.list"
```
