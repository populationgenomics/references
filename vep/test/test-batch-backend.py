#!/usr/bin/env python3

import click
import hail as hl
from cpg_utils import to_path
from cpg_workflows.batch import get_batch
from cpg_workflows.jobs.vep import add_vep_jobs


@click.command()
@click.argument('vep_version')
def main(vep_version: str):
    """
    Run VEP in parallel using Picard tools intervals as partitions.
    """
    vcf_path = f'gs://cpg-common-main/references-test/vep/test/sample.vcf.gz'
    out_ht_path = f'gs://cpg-common-main/references-test/vep/test/batch/sample-vep.ht'

    b = get_batch(f'Test VEP with Batch Backend, VEP v{vep_version}')
    add_vep_jobs(
        b=b,
        vcf_path=to_path(vcf_path),
        tmp_prefix=to_path('gs://cpg-common-main/references-test-tmp/vep/test'),
        out_path=to_path(out_ht_path),
        scatter_count=2,
    )
    b.run(wait=True)
    ht = hl.read_table(out_ht_path)
    ht.describe()
    ht.show()


if __name__ == '__main__':
    main()
