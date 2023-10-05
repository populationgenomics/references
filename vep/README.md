# Preparing a new VEP release

## Context

The following steps describe how to prepare reference data for obtaining and new version of VEP, 
and preparing the corresponding cache data. This is based on our desire to side-step the ability
of Hail to run VEP directly on MatrixTable data, which is currently limited to VEP version 95.

This approach is less streamlined than the process in [previous_process](previous_process/README.md), 
but it's cheaper, more flexible, and can be used to prepare VEP versions not available on Bioconda.

At time of writing (5/10/2023) the Ensembl cache is avaialble from a FTP server which only delivers 
~200kb/s, which necessitates using non-preeemptible VMs to avoid job termination over the >24hr 
runtime. This expensive process can be side-stepped using an asynchronus download facilitated by 
Globus, a secondary provider for accessing the same data.

This process describes a cheap, asynchronous download of the cache data, unwrapping locally, and 
syncing the contents to a GCP bucket where testing can take place prior to transfer of the relevant 
cache content into main/production.

## VEP Image

Create a new image in the populationgenomics `Images` repository, either using a new name 
(`ensembl-vep-VERSION`) or replacing the current build instructions (`vep`/`ensembl-vep`).

## Obtaining Cache files

1. Choose the VEP version you want to use. Make sure it's available on both the ensembl-vep
   DockerHub, and the Ensembl FTP server

 - DockerHub: e.g. https://hub.docker.com/r/ensemblorg/ensembl-vep/tags
 - FTP: e.g. https://ftp.ensembl.org/pub/release-110/variation/indexed_vep_cache/

2. Create/Log in to an account on [Globus](https://www.globus.org/) 
3. Download Globus Connect Personal and install it on your computer
4. Navigate to the Collection `Shared EMBL-EBI public endpoint`
5. Navigate to the Path `/gridftp/ensemblorg/pub/release-110/variation/indexed_vep_cache/` 
   (substituting for the version required)
6. Scroll to the Homo Sapiens cache file (~20GB in size. Choose the regular version, rather
   than the RefSeq or Merged versions)
7. Start transferring this file to a local directory. Once initiated this will continue
   in the background while your computer is on and the Globus app is running. You will
   receive an email once the transfer is completed.

## Processing Local Cache

1. Open the cache tarball 

```commandline
tar xzf homo_sapiens_vep_110_GRCh38.tar.gz
```

2. The resulting directory should be called `homo_sapiens`, containing a subdirectory called
   `VERSION_GRCh38`. Within that subdirectory there will be a directory for each contig, and
   a number of LRGs. Each of these folders should contain a range of region-specific cache files,
   as well as a single `all_vars.gz` and corresponding index. The `all_vars.gz` files are a result
   of downloading the indexed version of the vep cache, and result in much faster VEP runtimes.

3. Sync the data to a test location in GCP (after first checking that no data exists in the target
   location). This requires the `gcloud` library to be installed locally (see team-docs). Follow
   this structure to prevent any code changes when VEP cache is accessed by pipeline stages.
   The `cpg-common-test` bucket is typically writeable by staff.

```commandline
gcloud storage rsync --recursive homo_sapiens gs://cpg-common-test/references/vep/VERSION/mount/vep/homo_sapiens
```

4. Copy in additional files relating to the LOFTEE plugin and general running of VEP. The VEP
   stage uses cloudfuse to mount this whole path as a readable directory, so the relative location of
   the files is important to the successful running of vep

```bash
gcloud storage cp "gs://cpg-common-test/references/vep/110/mount/AlphaMissense_hg38.tsv.gz*" \\
   gs://cpg-common-test/references/vep/110/mount/gerp_conservation_scores.homo_sapiens.GRCh38.bw \\
   "gs://cpg-common-test/references/vep/110/mount/human_ancestor.fa.gz*" \\
   gs://cpg-common-test/references/vep/110/mount/loftee.sql \\
   gs://cpg-common-test/references/vep/VERSION/mount/
```

## Testing

This data should then be ready to test. The new VEP cache location can be chosen instead of the 
current standard VEP mount, the new VEP image can be used instead of the previous standard one,
and any other changes required can be actioned. Examples of other actions:

- Modifying the command sent to VEP (e.g. adding in new Plugins & corresponding data, or modifying output options) [link](https://github.com/populationgenomics/production-pipelines/blob/main/cpg_workflows/jobs/vep.py#L239)
- Updating the schema required when parsing the JSON output of VEP into a Hail Table [link](https://github.com/populationgenomics/production-pipelines/blob/main/cpg_workflows/query_modules/vep.py). 
  This can be trial and error, this schema definition will fail when an input is not specified, 
  but will be silent if a specified field is not provided.
- Updating AnnotateCohort to pull annotations from slightly different locations, or obtaining 
  values from a list instead of a single value (related to schema changes, [link](https://github.com/populationgenomics/production-pipelines/blob/main/cpg_workflows/query_modules/seqr_loader.py#L76))
