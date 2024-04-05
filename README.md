# Reference data

The aim of this repository is to track common reference resources for bioinformatics pipelines. The `references.py` script describes sources of reference resources, and the [CI workflow](https://github.com/populationgenomics/references/blob/main/.github/workflows/changed.yaml) uses it to populate the CPG references bucket and write a TOML template with fully qualified paths to resources for the analysis runner to use.

## Usage in analysis scripts

In analysis scripts that are run in the [analysis runner](https://github.com/populationgenomics/analysis-runner) environment, paths can be retrieved using the [reference_path](https://github.com/populationgenomics/cpg-utils/blob/main/cpg_utils/hail_batch.py#L252) helper function. For example, you can run the following code to retrieve the path to the GRCh38 reference fasta file:

```python
from cpg_utils.hail_batch import reference_path
path = reference_path('broad/ref_fasta')
```

In this case, `path` would resolve into `CloudPath('gs://cpg-common-main/references/hg38/v0/dragen_reference/Homo_sapiens_assembly38_masked.fasta')`

## Adding new sources

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

The script assumes the Google Cloud infrastructure, but the structure can be replicated for other cloud providers.

## Transfer Type

When adding a new source `transfer_cmd` can be specified to indicate the type of transfer that should be used to bring the resource(s) into our reference bucket. Without a specified `transfer_cmd` the config entries will still be populated, but no new transfer will be actioned. This can still be useful if the resource is already in the reference bucket, but the config entry is missing.

The transfer commands are actioned in CI using appropriate credentials, and the following commands are supported:

* `gcs_rsync`: Uses a recursive, non-destructive, `gcloud storage rsync` to copy the source to the destination. i.e. it doesn't delete files in the destination that are not in the source.
* `gcs_cp_r`: Uses a recursive `gcloud storage cp` to copy the source to the destination
* `curl`: uses a `curl` & `gcloud storage cp` to pull resources from an HTTP URL
