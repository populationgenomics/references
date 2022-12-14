#!/usr/bin/env python3

import click
import hail as hl
from cpg_workflows.batch import get_batch


@click.command()
@click.argument('vep_version')
def main(vep_version: str):
    """
    Submit `test-dataproc-script.py` to cluster.
    """
    from analysis_runner import dataproc

    mt_path = f'gs://cpg-common-main/references-test/vep/test/sample.vcf.mt'
    out_mt_path = (
        f'gs://cpg-common-main/references-test/vep/test/dataproc/sample-vep.vcf.mt'
    )

    b = get_batch(f'Test VEP with Spark backend, VEP v{vep_version}')
    dataproc.hail_dataproc_job(
        b,
        f'test-dataproc-script.py {mt_path} {out_mt_path}',
        max_age='4h',
        job_name=f'Test VEP {vep_version}',
        init=[f'gs://cpg-common-main/references/vep/{vep_version}/dataproc/init.sh'],
        worker_machine_type='n1-highmem-8',
        worker_boot_disk_size=200,
        secondary_worker_boot_disk_size=200,
        num_secondary_workers=20,
        num_workers=2,
    )
    b.run(wait=True)
    mt = hl.read_matrix_table(out_mt_path)
    mt.describe()
    mt.show()


if __name__ == '__main__':
    main()
