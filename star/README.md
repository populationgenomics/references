# Building a new STAR reference

## Context

The following steps describe how to prepare a reference package for the [STAR aligner](https://github.com/alexdobin/STAR). This is required for running the rare disease RNA sequencing analysis pipeline.

## STAR Image

A Dockerfile for generating a STAR image exists in the populationgenomics `Images` repository (`images/star/Dockerfile`). At the time of writing (15-12-2023), the current version of STAR being containerised is 2.7.10b. To update the version of STAR within the CPG artifact registry, a pull request will need to be made in the `Images` repository to change the defaut value of the `VERSION` variable in the Dockerfile, as well as updating the file `images.toml` in the root directory of the `Images` repository to reflect the new version of STAR. Once the changes have been merged into main, the new STAR image will be automatically built and pushed to artifact registry.

## Building the STAR reference

1. Update the file `star.toml` by supplying the new version of STAR to the `version` parameter.
2. (Optional) If a newer Gencode version is desired (currently version 44 at the time of writing), also update the `gencode_version` varialbe.
2. From the root directory of the `references` repository, run `analysis-runner` to start the script on Hail Batch:

```bash
STAR_VERSION="2.7.10b"  # Update to your current version number
analysis-runner \
    --dataset common \
    --description "Generate STAR reference" \
    --output-dir "references/star/${STAR_VERSION}/hg38 \
    --access-level standard \
    --config star/star.toml \
    --image australia-southeast1-docker.pkg.dev/cpg-common/images/cpg_workflows:latest \
    python3 star/build_star_reference.py
```

This will generate the reference package at `gs://cpg-common-main-test/references/star/<STAR_VERSION>/hg38`. This can be used for testing prior to pushing to production. Due to it's location within the `cpg-common-main-test` bucket, it will be deleted after a week.

## Production

Once testing is complete, the new STAR reference can be moved into production. The companion script
[prod_sync_script.sh](prod_sync_script.sh) can be used to move the data from the temporary test location into
the main bucket. The script takes one argument - the STAR version to push. It
does a minimal check that the source location `gs://cpg-common-main-test/references/star/${STAR_VERSION}/hg38`
exists, and the target location `gs://cpg-common-main/references/star/${STAR_VERSION}/hg38` does not. It 
then initiates the rsync between the two locations. This should be run using analysis-runner, at the
`full` access level:

```bash
STAR_VERSION="2.7.10b"  # Update to your current version number
analysis-runner \
    --dataset common \
    --access-level full \
    --description "STAR: Sync ${STAR_VERSION} to production" \
    --output-dir "references/star/${STAR_VERSION}/hg38" \
    bash star/prod_sync_script.sh ${STAR_VERSION}
```
