# Exome capture probesets/ design files

Sequencing libraries for exome datasets are constructed with capture probesets.
Exome capture probsets (or designs) can markedly differ in both the oligo probes
and the targeted regions of the genome, so we need to have access to exome designs suitable
for each exome data cohort (and sub-cohort cf. Mackenzie's mission).

[UCSC Genome browser](https://genome.ucsc.edu/index.html) has 'tracks' for many of the most common
exome probesets.

## Get exome probesets

use `download_ucsc_exomes.py` to retrieve the available exomes, and convert them from BigBed to Bed format.

Currently, the genome version `hg38`, and the api and download urls are hard coded constants.

** NB:**  
Some sets have two files 'target' and 'probes'. We should take Targets for use in workflows.
I have not added any padding to the bed files.


This script also produces `csc_exome_probsets_source_entry.txt`. This text file contains a `Source`
description to cut and paste into ../references.py, and lists a `bed` and `intervals_list` file for each exome probeset.


Note that the bed files need to bed converted with `picard` to 'interval_list' files. This will be done on  `gs://cpg-common-test`
using `../reference_generating_scripts/convert_bed_to_interval_list.py` or some modification of that script.

## Transfer exome probesets to `gs://cpg-common-test`

hg38 BEDs (default UCSC output, plus the Twist VCGS custom BED below):

gcloud storage cp *_hg38.bed gs://cpg-common-test/references/exome-probesets/hg38/

hg19 BEDs (intermediate inputs to the liftover job — currently only the Agilent
CRE v1 files below):

gcloud storage cp *_hg19.bed gs://cpg-common-test/references/exome-probesets/hg19/

## Run interval_conversion on main

```
analysis-runner \
    --dataset common \
    --access-level full \
    --description  'ref beds to interval_lists' \
    --output-dir cpg-common-main/references \
    convert_bedfiles_to_interval_lists.py

```

## Twist VCGS custom

Fetches Twist Comprehensive Exome + VCGS custom content (hg38) — the design used
for Mackenzie's Twist cohort. Output goes alongside the script; upload to
`cpg-common-test` afterwards.

```
bash get_twist_vcgs_custom.sh
gcloud storage cp Twist_VCGS_Exome_Covered_Targets_hg38.bed \
    gs://cpg-common-test/references/exome-probesets/hg38/
```

The matching hg38 `interval_list` is produced by the standard
`convert_bedfiles_to_interval_lists.py` AR invocation above.

## Agilent CRE v1 (hg19 → hg38 liftover)

Agilent SureSelect Clinical Research Exome v1 (design id `S06588914`) is only
distributed by UCSC at hg19, so a manual liftover step is required before the
hg38 BEDs land in references.

### 1. Download hg19 BEDs

```
python download_hg19_agilent_cre_v1.py
gcloud storage cp Agilent_ClinicalResearchExome_v1_*_hg19.bed \
    gs://cpg-common-test/references/exome-probesets/hg19/
```

### 2. Bootstrap hg19 sequence dictionary (one-time)

Once `hg19.fa.gz` has been transferred by CI (via the `hg19_fasta` Source in
`references.py`), generate the matching picard `.dict` alongside it. The job is
idempotent: it no-ops if `hg19.dict` is already present.

```
analysis-runner \
    --dataset common \
    --access-level full \
    --description  'Bootstrap hg19 picard sequence dictionary' \
    --output-dir cpg-common-main/references \
    generate_hg19_dict.py
```

### 3. Liftover + convert

Chains `picard BedToIntervalList → LiftOverIntervalList → IntervalListToBed` per
`*_hg19.bed` in `cpg-common-test`, writes the hg38 BED + interval_list to the
paths declared in `references.py` under `exome_probesets`, and logs the count of
rejected intervals so the operator can sanity-check liftover loss before
promoting the BEDs.

```
analysis-runner \
    --dataset common \
    --access-level full \
    --description  'liftover hg19 exome beds to hg38' \
    --output-dir cpg-common-main/references \
    liftover_and_convert_hg19_bedfiles.py
```

### Sanity-check post-liftover

Agilent's portal ships hg38 directly. Compare interval counts and total bp of
the lifted `_hg38.bed` outputs against the portal baseline at
`/Users/jossch/Downloads/S06588914/S06588914_{Regions,Covered}.bed` before
promoting the BEDs into ICA. Those files are not the source of truth — only a
reference for diffing.

## One-shot maintenance

### Rename liftover_37_to_38 chain typo

PR #96 landed the `liftover_37_to_38` chain at a typo'd path
(`grch37_to_grch38over.chain.gz` — missing a `.` before `over`). A follow-up
references PR corrects the Source `dst` to `grch37_to_grch38.over.chain.gz`,
and CI auto-curls the 228 KB chain to the corrected path on merge — leaving
the typo'd blob as an orphan in `cpg-common-main`. Run this script once
post-merge to delete it. The script is idempotent: it renames if only the
typo'd path exists, deletes the typo'd path if both exist, and no-ops if only
the corrected path is present.

```
analysis-runner \
    --dataset common \
    --access-level full \
    --description 'Rename liftover_37_to_38 chain typo' \
    --output-dir references/liftover-rename \
    rename_liftover_chain_typo.sh
```