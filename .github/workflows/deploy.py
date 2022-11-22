"""
* Read sources described in `references.py`
* Transfer files to the local bucket
* Prepare a TOML config with finalised paths
"""

import sys
import toml
import os
import click
from cpg_utils import to_path

from references import SOURCES, GENOME_BUILD

PROJECT = 'cpg-common'
PREFIX = 'gs://cpg-reference'


@click.command()
@click.option('--dry-run/--no-dry-run', is_flag=True, default=True)
@click.option('--project', default=PROJECT)
@click.option('--prefix', default=PREFIX)
def main(dry_run: bool, project, prefix):
    def _cmd(cmd):
        print(cmd, file=sys.stderr)
        if not dry_run:
            os.system(cmd)
    
    d = {'genome_build': GENOME_BUILD}
    for source in SOURCES:
        dst_path = to_path(prefix) / source.dst
        if source.src:
            _cmd(source.transfer_cmd(project, dst_path))
        if not source.files:
            d[source.name] = str(dst_path)
        else:
            d[source.name] = {
                k: str(dst_path / suffix)
                for k, suffix in source.files.items()
            }
    
    print(toml.dumps(d))


if __name__ == '__main__':
    main()
