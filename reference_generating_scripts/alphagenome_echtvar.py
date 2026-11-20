#!/usr/bin/env python3

"""
Build echtvar archives of AlphaGenome's precomputed SNV scores.

DeepMind publishes three tables of scores for every single-base change of GRCh38, and
download_alphagenome.py copied its three zips to the references bucket. echtvar
annotates a VCF with such values from an archive, and this script builds one archive
for each table given in --archives:

    score      alphagenome_avi            avi_scores_snvs_tabix.zip
    breakdown  alphagenome_avi_breakdown  avi_feature_importances_snvs_tabix.zip
    splicing   alphagenome_splicing       combined_splicing_snvs_tabix.zip

The conversion of a table's rows, the encode configuration and the check of a finished
archive are commands of the case annotation pipeline's tools image
(https://github.com/populationgenomics/case-annotation-pipeline, docs/annotation.md,
"Building an AlphaGenome archive"). This script runs them on the whole tables:

1. In the driver: read the zip's directory and the table's tabix index, which give
   each chromosome's byte range in the zip and its number of rows. The zip's
   generation is recorded and every later read asks for that generation.
2. One job per chromosome: read the chromosome's range, convert its rows with
   `case-annotation alphagenome-vcf` to a BCF, and encode the BCF with `echtvar
   encode` and the configuration of `case-annotation alphagenome-config`. The job
   fails unless it read, and echtvar encoded, as many rows as the index counts.
3. One merge job: copy the entries of the chromosomes' archives into one archive
   (alphagenome_echtvar_steps.py says why that is the archive of one encode), and
   compare every entry with its checksum.
4. One check job: read rows from several places of every chromosome, annotate their
   alleles with the merged archive and compare with `case-annotation
   alphagenome-check`. Only this job publishes.

One encode over all chromosomes works on one core, and for the largest table it was
estimated at 35 hours. The chromosomes' jobs run side by side instead.

Each archive is published into a folder of its name under the output directory, with
its SHA-256, the configuration that echtvar embedded in it, the report of the check,
a receipt (the source object with its generation and size, the rows of each
chromosome, the tools image and versions, the commit of this repository) and, for
the two tables that DeepMind offers for non-commercial use only, TERMS.txt. A run is
refused if a folder that it would publish into holds a file. Read DeepMind's terms
before building the breakdown or the splicing archive:
https://deepmind.google.com/science/alphagenome/terms

Run the score archive first and look at what it published, then the other two:

    analysis-runner --dataset common --access-level test \
        --output-dir references/alphagenome/echtvar \
        --description "AlphaGenome AVI score echtvar archive" \
        python3 reference_generating_scripts/alphagenome_echtvar.py \
        --archives score

    analysis-runner --dataset common --access-level test \
        --output-dir references/alphagenome/echtvar \
        --description "AlphaGenome breakdown and splicing echtvar archives" \
        python3 reference_generating_scripts/alphagenome_echtvar.py \
        --archives breakdown splicing

With --access-level test the archives go to gs://cpg-common-test/<output-dir>/, and
with full access to gs://cpg-common-main/<output-dir>/. The zips are read from
cpg-common-main at every level. To print the jobs without submitting them, add a
--config file to the analysis-runner call that sets dry_run = true under [hail],
which cpg_utils passes on to the batch.

The jobs of the chromosomes run on spot VMs: Hail Batch starts a reclaimed job again,
and the job then loses one chromosome's work. The merge and the check run on standard
VMs, because they come last and hold the largest disks. --no-spot runs every job on a
standard VM.
"""

import dataclasses
import json
import logging
import math
import re
import subprocess
from argparse import ArgumentParser
from pathlib import Path

from cpg_utils import to_path
from cpg_utils.config import config_retrieve, dataset_path, output_path, reference_path
from cpg_utils.hail_batch import get_batch

import alphagenome_echtvar_steps as steps

# Where download_alphagenome.py put DeepMind's zips, under DeepMind's names.
SOURCE_FOLDER = 'gs://cpg-common-main/references/alphagenome'

# The case annotation pipeline's tools image, by digest, and the commit of the
# pipeline that it was built from, which is the image's tag. Change both together.
# The image holds case-annotation, bcftools 1.21 and echtvar 0.2.4.
TOOLS_IMAGE = (
    'australia-southeast1-docker.pkg.dev/cpg-common/images-dev/case-annotation-tools'
    '@sha256:fff349e6d2cc21833ec541eca4703103c7fd9c7f925188473de66845a862225b'
)
TOOLS_COMMIT = '46cd586250925d71fe525580d32c0d01796279be'

# The notice that DeepMind's terms give as an example for output under its Output
# Terms of Use, word for word.
NOTICE = (
    'By using this information, you agree to AlphaGenome Output Terms of Use found at '
    'https://deepmind.google.com/science/alphagenome/output-terms'
)

GIB = 2**30
# The two estimates behind the disks, each above what sample rows of the three tables
# came to. A chromosome's BCF, in bytes for each byte of the chromosome's range in the
# zip: 1.2 to 1.5. An archive, likewise: 0.5 to 0.7 for rows of a few regions, and up
# to 1.0 for rows from many separate places, which echtvar stores less densely.
BCF_PER_SOURCE_BYTE = 1.6
ARCHIVE_PER_SOURCE_BYTE = 1.2
# Beside a job's files, in GiB.
SPARE_STORAGE = 10
# Seconds after which Hail Batch stops a job. The largest chromosome of the slowest
# table is estimated at 5 hours, from rates measured on a laptop under emulation.
ENCODE_TIMEOUT = 24 * 3600
MERGE_TIMEOUT = 12 * 3600
CHECK_TIMEOUT = 6 * 3600


@dataclasses.dataclass(frozen=True)
class Archive:
    """
    One of DeepMind's tables and the archive that is built from it.

    `name` is the archive's name, which is the table's name in the pipeline's
    commands. `title` is the name of the download on DeepMind's page. `notice` is
    True for the two tables that the page lists under "Downloadable artifacts for
    non-commercial use only".
    """

    name: str
    zip: str
    title: str
    notice: bool


ARCHIVES = {
    'score': Archive(
        'alphagenome_avi', 'avi_scores_snvs_tabix.zip', 'AVI SNV scores', False
    ),
    'breakdown': Archive(
        'alphagenome_avi_breakdown',
        'avi_feature_importances_snvs_tabix.zip',
        'AVI SNV feature importance scores',
        True,
    ),
    'splicing': Archive(
        'alphagenome_splicing',
        'combined_splicing_snvs_tabix.zip',
        'AlphaGenome 1 SNV merged splicing scores',
        True,
    ),
}


def storage_gib(needed_bytes: float) -> int:
    """The storage to request for a job whose files take `needed_bytes`, in GiB."""
    return SPARE_STORAGE + math.ceil(needed_bytes / GIB)


def encode_storage(span: dict) -> int:
    """
    GiB for the job of one chromosome. The range is read as a stream, so the disk
    holds the BCF and the chromosome's archive.
    """
    return storage_gib((BCF_PER_SOURCE_BYTE + ARCHIVE_PER_SOURCE_BYTE) * span['bytes'])


def merge_storage(plan: dict) -> int:
    """GiB for the merge job: the chromosomes' archives and the merged one."""
    return storage_gib(2 * ARCHIVE_PER_SOURCE_BYTE * plan['member_bytes'])


def check_storage(plan: dict) -> int:
    """GiB for the check job: the merged archive."""
    return storage_gib(ARCHIVE_PER_SOURCE_BYTE * plan['member_bytes'])


def repository_commit() -> str:
    """The commit of this repository that the analysis-runner checked out."""
    here = Path(__file__).resolve().parent
    return subprocess.check_output(
        ['git', 'rev-parse', 'HEAD'], cwd=here, text=True
    ).strip()


def refuse_existing(folder: str) -> None:
    """
    An archive's folder is written by one run. A second run would publish among the
    first one's files, and a file that it does not write again would stay beside
    its archive.

    Raises:
        ValueError: if an object lies under `folder`
    """
    if to_path(folder).exists():
        held = sorted(str(path) for path in to_path(folder).iterdir())
        raise ValueError(
            f'{folder} already holds {", ".join(held[:5])}. A run publishes into an '
            f'empty folder only: give another --output-dir to the analysis-runner, '
            f'or leave this archive out of --archives'
        )


def make_plan(
    archive: Archive,
    source_folder: str,
    reference_index: str,
    check_places: int,
    tools_image: str,
    tools_commit: str,
    commit: str,
) -> dict:
    """
    The plan of one archive: where the rows of each chromosome lie in DeepMind's zip
    (alphagenome_echtvar_steps.locate), and what the receipt says of the run.

    Args:
        archive: the archive to build
        source_folder: the folder of DeepMind's zips
        reference_index: the .fai of the reference
        check_places: how many places of a chromosome the check reads rows from
        tools_image: the tools image, by digest
        tools_commit: the commit of the pipeline that the image was built from
        commit: the commit of this repository
    """
    source = steps.open_source(f'{source_folder.rstrip("/")}/{archive.zip}')
    plan = steps.locate(source, check_places)
    for span in plan['chromosomes']:
        # The jobs' commands and file names hold the name as the index gives it.
        if not re.fullmatch(r'[0-9A-Za-z_.]+', span['name']):
            raise ValueError(
                f'{source.described.url}: the index has a sequence '
                f'{span["name"]!r}, which is not the name of a chromosome'
            )
    plan['table'] = {
        'name': archive.name,
        'title': archive.title,
        'notice': NOTICE if archive.notice else None,
    }
    plan['reference_index'] = dataclasses.asdict(
        steps.open_source(reference_index).described
    )
    plan['tools'] = {'image': tools_image, 'pipeline_commit': tools_commit}
    plan['script'] = {
        'repository': 'populationgenomics/references',
        'script': f'reference_generating_scripts/{Path(__file__).name}',
        'commit': commit,
        'analysis_runner': config_retrieve(['workflow', 'ar-guid']),
        'access_level': config_retrieve(['workflow', 'access_level']),
    }
    return plan


def encode_job(batch, plan: dict, span: dict, inputs: dict, image: str, spot: bool):
    """
    The job that encodes one chromosome's archive.

    Args:
        batch: the batch
        plan: the archive's plan
        span: the chromosome's entry in the plan
        inputs: the plan, this script's steps and the reference index, as inputs
        image: the tools image
        spot: False to run the job on a standard VM
    """
    table, chromosome, rows = plan['table']['name'], span['name'], span['rows']
    job = batch.new_bash_job(
        f'Encode {table} {chromosome}', attributes={'tool': 'echtvar'}
    )
    job.image(image)
    job.cpu(4)
    job.memory('standard')
    job.storage(f'{encode_storage(span)}Gi')
    job.timeout(ENCODE_TIMEOUT)
    if not spot:
        job.spot(False)
    job.command(
        f"""
set -euo pipefail
mkdir -p "$BATCH_TMPDIR/work"
cd "$BATCH_TMPDIR/work"
case-annotation alphagenome-config --table {table} > config.json
python3 {inputs['steps']} rows --plan {inputs['plan']} --chromosome {chromosome} \\
    | case-annotation alphagenome-vcf --table {table} \\
        --reference-index {inputs['reference_index']} \\
    | bcftools view --no-version --threads 2 -Ob -o {chromosome}.bcf
if ! echtvar encode {job.archive} config.json {chromosome}.bcf 2> encode.log; then
    cat encode.log >&2
    exit 1
fi
cat encode.log >&2
# echtvar ends with the number of records that it encoded.
if ! grep -q -F 'wrote {rows} total variants' encode.log; then
    echo 'error: echtvar did not encode the {rows} rows of {chromosome}' >&2
    exit 1
fi
"""
    )
    return job


def merge_job(batch, plan: dict, encoded: list, inputs: dict, image: str):
    """
    The job that merges the chromosomes' archives, in the table's order, and writes
    the SHA-256 of the merged archive under the name that the archive is published
    with.
    """
    table = plan['table']['name']
    job = batch.new_bash_job(f'Merge {table}', attributes={'tool': 'echtvar'})
    job.image(image)
    job.cpu(2)
    job.memory('standard')
    job.storage(f'{merge_storage(plan)}Gi')
    job.timeout(MERGE_TIMEOUT)
    job.spot(False)
    archives = ' '.join(str(one.archive) for one in encoded)
    job.command(
        f"""
set -euo pipefail
python3 {inputs['steps']} merge --out {job.archive} {archives}
digest=$(sha256sum < {job.archive} | cut -d' ' -f1)
printf '%s  %s\\n' "$digest" {table}.zip > {job.sha256}
"""
    )
    return job


def check_job(batch, plan: dict, merged, inputs: dict, image: str):
    """
    The job that checks the merged archive against rows of the table and, when the
    check has passed, hands on what is published.
    """
    table = plan['table']['name']
    job = batch.new_bash_job(f'Check {table}', attributes={'tool': 'echtvar'})
    job.image(image)
    job.cpu(2)
    job.memory('standard')
    job.storage(f'{check_storage(plan)}Gi')
    job.timeout(CHECK_TIMEOUT)
    job.spot(False)
    terms = f'--terms {job.terms}' if plan['table']['notice'] else ''
    job.command(
        f"""
set -euo pipefail
mkdir -p "$BATCH_TMPDIR/work"
cd "$BATCH_TMPDIR/work"
digest=$(sha256sum < {merged.archive} | cut -d' ' -f1)
if [ "$digest" != "$(cut -d' ' -f1 {merged.sha256})" ]; then
    echo 'error: the archive is not the one that the merge job wrote' >&2
    exit 1
fi
python3 {inputs['steps']} sample --plan {inputs['plan']} > sample.tsv
case-annotation alphagenome-vcf --table {table} \\
        --reference-index {inputs['reference_index']} < sample.tsv \\
    | bcftools annotate --no-version -x INFO -Ob -o sample.bcf
echtvar anno -e {merged.archive} sample.bcf sample.annotated.bcf
bcftools index sample.annotated.bcf
case-annotation alphagenome-check --table {table} \\
    --annotated sample.annotated.bcf < sample.tsv > check.txt
cat check.txt
python3 {inputs['steps']} records --plan {inputs['plan']} \\
    --archive {merged.archive} --sha256 "$digest" --check check.txt \\
    --config {job.config} --receipt {job.receipt} {terms}
cp check.txt {job.check}
cp {merged.sha256} {job.sha256}
mv {merged.archive} {job.archive}
"""
    )
    return job


def add_archive(
    batch, plan: dict, plan_path: str, steps_path: str, image: str, spot: bool
) -> None:
    """
    Adds the jobs of one archive to the batch, and what the last of them publishes.

    Args:
        batch: the batch
        plan: the archive's plan
        plan_path: where the plan lies, for the jobs to read
        steps_path: where alphagenome_echtvar_steps.py lies, for the jobs to run
        image: the tools image
        spot: False to run the chromosomes' jobs on standard VMs
    """
    table = plan['table']['name']
    inputs = {
        'plan': batch.read_input(plan_path),
        'steps': batch.read_input(steps_path),
        'reference_index': batch.read_input(plan['reference_index']['url']),
    }
    encoded = [
        encode_job(batch, plan, span, inputs, image, spot)
        for span in plan['chromosomes']
    ]
    merged = merge_job(batch, plan, encoded, inputs, image)
    checked = check_job(batch, plan, merged, inputs, image)
    published = {
        f'{table}.zip': checked.archive,
        f'{table}.zip.sha256': checked.sha256,
        f'{table}.config.json': checked.config,
        f'{table}.check.txt': checked.check,
        'receipt.json': checked.receipt,
    }
    if plan['table']['notice']:
        published['TERMS.txt'] = checked.terms
    for name, resource in published.items():
        batch.write_output(resource, output_path(f'{table}/{name}'))


def main(
    archives: list[str],
    tools_image: str,
    tools_commit: str,
    source_folder: str,
    reference_index: str,
    check_places: int,
    spot: bool,
) -> None:
    """
    Plans the archives given and submits their jobs.

    Args:
        archives: which archives to build, by the keys of ARCHIVES
        tools_image: the case annotation pipeline's tools image, by digest
        tools_commit: the commit of the pipeline that the image was built from
        source_folder: the folder of DeepMind's zips
        reference_index: the .fai of the reference that VCFs are called against
        check_places: how many places of a chromosome the check reads rows from
        spot: False to run every job on a standard VM
    """
    if len(set(archives)) != len(archives):
        raise ValueError(f'--archives has an archive twice: {" ".join(archives)}')
    if not re.fullmatch(r'[^\s@]+@sha256:[0-9a-f]{64}', tools_image):
        raise ValueError(
            f'--tools-image {tools_image} is not an image by its digest, '
            f'<image>@sha256:<digest>'
        )
    if not re.fullmatch(r'[0-9a-f]{40}', tools_commit):
        raise ValueError(f'--tools-commit {tools_commit} is not a full commit id')
    chosen = [ARCHIVES[key] for key in archives]
    for archive in chosen:
        refuse_existing(output_path(archive.name))

    commit = repository_commit()
    run_tmp = dataset_path(
        f'alphagenome_echtvar/{config_retrieve(["workflow", "ar-guid"])}',
        category='tmp',
    )
    steps_path = f'{run_tmp}/{Path(steps.__file__).name}'
    with to_path(steps_path).open('w') as handle:
        handle.write(Path(steps.__file__).read_text())

    batch = get_batch(name=f'AlphaGenome echtvar archives: {" ".join(archives)}')
    for archive in chosen:
        plan = make_plan(
            archive,
            source_folder,
            reference_index,
            check_places,
            tools_image,
            tools_commit,
            commit,
        )
        plan_path = f'{run_tmp}/{archive.name}.plan.json'
        with to_path(plan_path).open('w') as handle:
            json.dump(plan, handle, indent=2)
        disks = [encode_storage(span) for span in plan['chromosomes']]
        logging.info(
            f'{archive.name}: {plan["rows"]:,} rows on {len(plan["chromosomes"])} '
            f'chromosomes, {plan["member_bytes"]:,} bytes in generation '
            f'{plan["source"]["generation"]} of {plan["source"]["url"]}; plan at '
            f'{plan_path}; disks of {min(disks)} to {max(disks)} GiB for the '
            f'chromosomes, {merge_storage(plan)} GiB for the merge and '
            f'{check_storage(plan)} GiB for the check; publishing to '
            f'{output_path(archive.name)}'
        )
        add_archive(batch, plan, plan_path, steps_path, tools_image, spot)
    batch.run(wait=False)


if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO)
    parser = ArgumentParser(description=__doc__.split('\n\n')[0])
    parser.add_argument(
        '--archives',
        nargs='+',
        required=True,
        choices=list(ARCHIVES),
        help='Which archives to build',
    )
    parser.add_argument(
        '--tools-image',
        default=TOOLS_IMAGE,
        help='The case annotation pipeline tools image, as <image>@sha256:<digest>',
    )
    parser.add_argument(
        '--tools-commit',
        default=TOOLS_COMMIT,
        help='The commit of the pipeline that the tools image was built from',
    )
    parser.add_argument(
        '--source-folder',
        default=SOURCE_FOLDER,
        help="The folder that holds DeepMind's three zips under DeepMind's names",
    )
    parser.add_argument(
        '--reference-index',
        help="The .fai of the reference. Default: reference_path('broad/ref_fasta') "
        "with '.fai' appended",
    )
    parser.add_argument(
        '--check-places',
        type=int,
        default=8,
        help="How many of a chromosome's 16 kb windows the check reads rows from, "
        "beside the chromosome's first rows",
    )
    parser.add_argument(
        '--no-spot',
        action='store_true',
        help='Run every job on a standard VM. By default the jobs of the '
        'chromosomes run on spot VMs',
    )
    args = parser.parse_args()
    main(
        archives=args.archives,
        tools_image=args.tools_image,
        tools_commit=args.tools_commit,
        source_folder=args.source_folder,
        reference_index=args.reference_index
        or f'{reference_path("broad/ref_fasta")}.fai',
        check_places=args.check_places,
        spot=not args.no_spot,
    )
