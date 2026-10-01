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
* `gcs_cp_single`: Copies one object with `gcloud storage cp`
* `curl`: uses a `curl` & `gcloud storage cp` to pull one object from an HTTP URL. A failed download removes the partial object so the next run retries it.
* `curl_with_user_agent`: `curl` with a browser user agent, for hosts that refuse curl's default

A source with `files` transfers exactly those entries and nothing else under `src`: each is copied on its own to `dst/<suffix>`, and any missing one is enough to schedule the source again. For a `gs://` source a directory-like entry (`.ht`, `.mt`, `.vds`) is rsynced and anything else is copied as a single object, whatever `transfer_cmd` says; for an HTTP source `src` must be the directory URL (ending in `/`) and each entry is curled. Only listed files reach the config, so only listed files are copied and paid for; the rest of an upstream prefix is never pulled. A source without `files` copies `src` whole with its `transfer_cmd`. A single-object `src` ending in `/` without `files` is rejected on import.

A file name holding a `{placeholder}` (such as `gatk_sv`'s `shard-{shard}.tar.gz`) stands for a family of files the consumer expands with `.format()`. Its parent folder is rsynced and checked for instead, once however many templates share it. The placeholder may only appear in the last path component.

A source scheduled because an entry is missing copies only the missing entries, so a listed entry that has since disappeared upstream does not fail the run while our copy exists. A source whose `src` or `dst` changed is copied in full.
