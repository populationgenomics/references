# GWAS summary statistics

Public GWAS summary statistics for clinical biomarkers, converted to one GRCh38
layout so any project can use them (allele-frequency lookups against 1000 Genomes or
OurDNA, colocalisation with TenK10K, ...).

| File | What it does |
| --- | --- |
| `files.csv` | The 336 source files: URL, published MD5 (or where to find it), genome build and the evidence for it, study, biomarker, ancestry, N |
| `download.py` | Downloads each file unchanged into `gs://cpg-common-main-tmp/gwas_sumstats/original/` and checks its MD5 |
| `format.py` | Converts each download to GRCh38 and the standard columns, writing to `gs://cpg-common-main/references/gwas_sumstats/v1/`; then writes `manifest.csv` there and `manifest.html` to `gs://cpg-common-main-web/gwas_sumstats/` |

Every source is the authors' own file, never the GWAS Catalog's harmonised copy, so
every change to the data is made, and counted, by `format.py`. Older versions of the
harmoniser lifted GRCh37 positions with an off-by-one error for some variants
([sumstats-harmoniser#52](https://github.com/gwas-catalog/sumstats-harmoniser/pull/52)).
GRCh37 and NCBI36 sources are lifted by `format.py`; GRCh38 sources are used as
published. `build_evidence` in `files.csv` says how each build is known and, for
GRCh38 sources, whether the data were made natively on GRCh38 (sequencing or a
GRCh38 imputation panel) where the GWAS Catalog metadata says so. Three studies
(UKB-WGS 2025, Timsina 2026, Wei 2024) publish no rsID column, so their `rsid` is NA.

Objects in `-main-tmp` are deleted after 8 days, so run `format.py` within a week of
`download.py`. Each format job keeps its stats (row counts, drop reasons, and the
download's MD5 and time) in `v1/stats/`, so `--manifest` can be rerun at any time. A
rerun of either script skips files already written.

## Running

Both scripts need `--access-level full` (writes to main), so they run from `main`
after merge.

Prerequisites: the Google Cloud CLI (`gcloud`, https://cloud.google.com/sdk/docs/install)
and full access to the `common` dataset in analysis-runner. Who holds that is not
listed in `common`'s `members.yaml`; ask Software Platforms if analysis-runner
refuses the job.

1. Merge the PR, then wait for the **Deploy resources and config** action on `main` to
   go green (GitHub, Actions tab). It copies the NCBI36 chain file and adds
   `liftover_36_to_38` to the references config; `format.py` fails until it has run.

2. From a clone of this repo:

   ```bash
   gcloud auth login
   python3 -m venv .venv
   source .venv/bin/activate
   pip install analysis-runner
   ```

3. Download:

   ```bash
   analysis-runner --dataset common --access-level full --output-dir gwas_sumstats --description "Download biomarker GWAS summary statistics" python3 gwas_sumstats/download.py
   ```

   The analysis-runner job only submits a second batch, one job per file, and
   finishes within minutes. Its log (open the link analysis-runner prints) ends with
   the link to that second batch,
   `https://batch.hail.populationgenomics.org.au/batches/<id>`: that page shows when
   every download is done. If some failed, run the same command again; files already
   downloaded are skipped.

4. Format, within 8 days of the download, and check the batch page the same way.
   If the 8 days have passed, the originals are gone from `-main-tmp`: run step 3
   again, then this step, which skips files already formatted and lists any it
   could not submit. Step 3 downloads every original no longer in `-main-tmp`, so
   add `--only <file_id> ...` to it to fetch just the files still to format.

   ```bash
   analysis-runner --dataset common --access-level full --output-dir gwas_sumstats --description "Format biomarker GWAS summary statistics" python3 gwas_sumstats/format.py
   ```

5. Manifest, once every format job has finished:

   ```bash
   analysis-runner --dataset common --access-level full --output-dir gwas_sumstats --description "Biomarker GWAS summary statistics manifest" python3 gwas_sumstats/format.py --manifest
   ```

`--output-dir` is required by analysis-runner but unused: outputs go to the buckets
listed above.

Add `--only <file_id> ...` to the download or format command to run a subset. The
manifest command takes no `--only`: it always indexes every formatted file, so rerun
it as is after a partial download or format. It fails, listing them, if any
`files.csv` row is not formatted; `--allow-missing` writes it anyway and lists the
missing files in `manifest.html`. `--local` runs on your
machine instead of Batch (see Running locally, below). After correcting a row in
`files.csv`, rerun `format.py --only <file_id> --force` to redo that file: without
`--force`, files already formatted are skipped.

## Running locally

To download and format a few files on your machine, without the cloud. Use a folder
outside the repo (here `~/gwas_sumstats_local`) so nothing gets committed by
accident. From the repo root:

```bash
mkdir -p ~/gwas_sumstats_local/ref
python3 -m venv ~/gwas_sumstats_local/venv
source ~/gwas_sumstats_local/venv/bin/activate
pip install polars==1.34.0 pysam==0.23.3 pyliftover==0.4.1 numpy
curl -fL -o ~/gwas_sumstats_local/ref/GRCh38.fasta https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta
curl -fL -o ~/gwas_sumstats_local/ref/GRCh38.fasta.fai https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta.fai
curl -fL -o ~/gwas_sumstats_local/ref/hg19ToHg38.over.chain.gz https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz
curl -fL -o ~/gwas_sumstats_local/ref/hg18ToHg38.over.chain.gz https://hgdownload.soe.ucsc.edu/goldenPath/hg18/liftOver/hg18ToHg38.over.chain.gz
```

Then, for one file (Wheeler 2017 is small, 18 MB, and NCBI36, so it exercises the
liftover):

```bash
python3 gwas_sumstats/download.py --local --out ~/gwas_sumstats_local/original --only 2017_Wheeler_PLoSMed_HbA1c_SAS_GCST007951
python3 gwas_sumstats/format.py --local --originals ~/gwas_sumstats_local/original --out ~/gwas_sumstats_local/v1 --stats ~/gwas_sumstats_local/stats --web ~/gwas_sumstats_local/web --fasta ~/gwas_sumstats_local/ref/GRCh38.fasta --chain-grch37 ~/gwas_sumstats_local/ref/hg19ToHg38.over.chain.gz --chain-ncbi36 ~/gwas_sumstats_local/ref/hg18ToHg38.over.chain.gz --only 2017_Wheeler_PLoSMed_HbA1c_SAS_GCST007951
```

- The fasta is Google's public copy of the Broad GRCh38 reference (3.2 GB, about
  10 minutes; UCSC contig names `chr1`). The cloud jobs use the Ensembl 113 primary
  assembly from the references config (`ensembl_113/unmasked_reference`, contigs
  `1`); `format.py` accepts either, and the primary chromosomes are the same
  sequence. Ensembl's own FTP is much slower, and copying the references bucket's
  copy to a laptop is billed as egress.
- The chain files are the same UCSC files the references config points at
  (`liftover_37_to_38`, `liftover_36_to_38`). Each is needed only if a selected
  file is on that build; `format.py --local` says which flag is missing.
- `python3 gwas_sumstats/format.py --check-columns` needs none of these: it reads
  the start of each source straight from its URL.

## Tests

`test_format.py` checks the transforms (OR to beta, P values, strand flips, liftover)
on small made-up inputs. From the repo root:

```bash
uv run --no-project --with pytest --with polars==1.34.0 --with pysam==0.23.3 --with pyliftover==0.4.1 --with numpy pytest gwas_sumstats
```

## Output

One `<file_id>_GRCh38_formatted.tsv.gz` (bgzipped, tabix-indexed) per source file,
with these columns:

`chromosome  base_pair_location  effect_allele  other_allele  beta  standard_error
effect_allele_frequency  p_value  neg_log_10_p_value  rsid  n  z`

- `chromosome` uses GWAS-SSF codes: 1-22, X=23, Y=24, MT=25.
- Positions are GRCh38, 1-based. GRCh37 and NCBI36 sources are lifted with the UCSC
  chain files (`liftover_37_to_38`, `liftover_36_to_38`).
- Both alleles are checked against GRCh38; SNVs reported on the other strand are
  complemented. A file whose A/T and C/G SNPs mostly fail this check is on the wrong
  build, and formatting stops rather than writing it.
- `beta` is per effect allele; odds ratios become ln(OR), and a source with only z
  and SE gets beta = z x SE. A source with no effect size keeps `beta` as NA (none
  in the current list).
- Indels coded without sequence (`D`/`I`, `Y`/`Z`) cannot be placed on GRCh38 and
  are dropped (`dropped_non_acgt_allele` in the manifest).
- `effect_allele_frequency` is the study's own frequency. The Wheeler 2017 HbA1c
  files give only HapMap reference panel frequencies, so theirs is NA.
- `p_value` keeps the source text, so values below the float range survive. Only
  numeric text in (0, 1] is kept. A P of exactly 0 (an underflowed top hit) becomes
  NA, counted as `p_value_zero`, with beta and SE kept. Anything else (`<1e-300`,
  `-0.01`) becomes NA as `p_value_unparseable`, and a file with more than 1% of
  those fails.
- Space-separated files keep empty fields in place (a missing rsID is two spaces),
  so their rows are not shifted.
- `n` is per variant where the source gives it, otherwise the study total.
- Nothing is filtered on frequency or P value. `manifest.csv` counts every row dropped
  and why.

`file_id` is `<year>_<first author>_<journal>_<biomarker>_<ancestry>_<accession>`,
where the accession is the GWAS Catalog study ID or the source's dataset ID.

## Not included

Sources that need a login or an access application are not in `files.csv`: UKB-PPP
proteins (Synapse), ThyroidOmics (click-through terms), tau PET (ADNI), and the
dbGaP-controlled and on-request studies.

The Chen 2020 trans-ethnic blood-count files come from the authors' site
(https://www.mhi-humangenetics.org/en/resources), not the GWAS Catalog: the Catalog
holds only their MR-MEGA results (P value, no effect size), while the authors also
publish a fixed-effect GWAMA meta-analysis of the same data with beta and SE.

The Chen 2021 trans-ancestry glycaemic files (GCST90002229, GCST90002235,
GCST90002241) are left out: they hold MR-MEGA Bayes factors, with no association P
value or effect size. The ancestry-specific files from the same paper are included.

The Wheeler 2017 trans-ethnic HbA1c file (GCST004903) is left out for the same reason:
it holds MANTRA Bayes factors. The ancestry-specific files from the same paper are
included.
