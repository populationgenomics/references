# GWAS summary statistics

Public GWAS summary statistics for clinical biomarkers, converted to one GRCh38
layout so any project can use them (allele-frequency lookups against 1000 Genomes or
OurDNA, colocalisation with TenK10K, ...).

| File | What it does |
| --- | --- |
| `files.csv` | The 336 source files: URL, published MD5 (or where to find it), genome build and the evidence for it, study, biomarker, ancestry, N |
| `download.py` | Downloads each file unchanged into `gs://cpg-common-main-tmp/gwas_sumstats/original/` and checks its MD5 |
| `format.py` | Converts each download to GRCh38 and the standard columns, writing to `gs://cpg-common-main/references/gwas_sumstats/v1/`; then writes `manifest.csv` there and `manifest.html` to `gs://cpg-common-main-web/gwas_sumstats/` |

Objects in `-main-tmp` are deleted after 8 days, so run `format.py` within a week of
`download.py`. Each format job keeps its stats (row counts, drop reasons, and the
download's MD5 and time) in `v1/stats/`, so `--manifest` can be rerun at any time. A
rerun of either script skips files already written.

## Running

Both scripts need `--access-level full` (writes to main), so they run from `main`
after merge.

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

4. Format, within 8 days of the download, and check the batch page the same way:

   ```bash
   analysis-runner --dataset common --access-level full --output-dir gwas_sumstats --description "Format biomarker GWAS summary statistics" python3 gwas_sumstats/format.py
   ```

5. Manifest, once every format job has finished:

   ```bash
   analysis-runner --dataset common --access-level full --output-dir gwas_sumstats --description "Biomarker GWAS summary statistics manifest" python3 gwas_sumstats/format.py --manifest
   ```

`--output-dir` is required by analysis-runner but unused: outputs go to the buckets
listed above.

Add `--only <file_id> ...` to any of these to run a subset. `--local` runs on your
machine instead of Batch (see each script's docstring). After correcting a row in
`files.csv`, rerun `format.py --only <file_id> --force` to redo that file: without
`--force`, files already formatted are skipped.

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
- `beta` is per effect allele; odds ratios become ln(OR). Sources with no effect
  size (the Chen 2020 trans-ethnic blood-count files) keep `beta` as NA.
- Indels coded without sequence (`D`/`I`, `Y`/`Z`) cannot be placed on GRCh38 and
  are dropped (`dropped_non_acgt_allele` in the manifest).
- `effect_allele_frequency` is the study's own frequency. The Wheeler 2017 HbA1c
  files give only HapMap reference panel frequencies, so theirs is NA.
- `p_value` keeps the source text, so values below the float range survive.
- `n` is per variant where the source gives it, otherwise the study total.
- Nothing is filtered on frequency or P value. `manifest.csv` counts every row dropped
  and why.

`file_id` is `<year>_<first author>_<journal>_<biomarker>_<ancestry>_<accession>`,
where the accession is the GWAS Catalog study ID or the source's dataset ID.

## Not included

Sources that need a login or an access application are not in `files.csv`: UKB-PPP
proteins (Synapse), ThyroidOmics (click-through terms), tau PET (ADNI), and the
dbGaP-controlled and on-request studies.

The Chen 2021 trans-ancestry glycaemic files (GCST90002229, GCST90002235,
GCST90002241) are left out: they hold MR-MEGA Bayes factors, with no association P
value or effect size. The ancestry-specific files from the same paper are included.

The Wheeler 2017 trans-ethnic HbA1c file (GCST004903) is left out for the same reason:
it holds MANTRA Bayes factors. The ancestry-specific files from the same paper are
included.
