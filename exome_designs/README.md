# Exome capture probesets/ design files

Sequencing libraries for exome datasets are constructed with capture probesets.
Exome capture probsets (or designs) can markedly differ in both the oligo probes
and the targetd regions of the genome, so we need to have access to exome designs suitable
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
gcloud storage cp *.bed gs://cpg-common-test/references/exome-probesets/hg38/ 

## TODO

the next stages will include:
convert to `intervals_list`
copy folder on `test` to `main`