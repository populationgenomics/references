"""
Copy references in TOML and prepare a ready TOML
"""

import sys
import toml
import os
import click
from cpg_utils import to_path


PROJECT = 'cpg-common'
PREFIX = 'gs://cpg-reference'
TOML_PATH = '../../references.toml'


@click.command()
@click.option('--dry-run/--no-dry-run', is_flag=True, default=True)
@click.option('--project', default=PROJECT)
@click.option('--prefix', default=PREFIX)
@click.option('--toml-path', default=TOML_PATH)
def main(dry_run: bool, project, prefix, toml_path):
    def _cmd(cmd):
        print(cmd, file=sys.stderr)
        if not dry_run:
            os.system(cmd)
    
    d = toml.load(toml_path)
    ready_d = {}
    for section_k, section in d.items():
        if not isinstance(section, dict):
            ready_d[section_k] = section
            continue
        dst_path = to_path(prefix) / section.pop('destination')
        if source := section.pop('source', None):
            if source.startswith('gs://'):
                _cmd(
                    f'gsutil -u {project} -m rsync -d -r '
                    f'{source} {dst_path}'
                )
            if source.startswith('https://'):
                _cmd(
                    f'curl {source} -o tmp && '
                    f'gsutil -u {project} cp tmp {dst_path}'
                )
        if not section:
            ready_d[section_k] = str(dst_path)
        else:
            ready_d[section_k] = {
                k: str(dst_path / suffix)
                for k, suffix in section.items()
            }
    
    print(toml.dumps(ready_d))


if __name__ == '__main__':
    main()
