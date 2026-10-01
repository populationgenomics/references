"""
Transfer one reference source
"""

import argparse
import subprocess
import sys

from prep_matrix import gcs_file_exists
from references import SOURCES


def parser(args: list[str]) -> argparse.Namespace:
    """Parse arguments for the script

    Args:
        args (list[str]): _description_

    Returns:
        argparse.Namespace: _description_
    """
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        '--references-prefix',
        required=True,
        help='Prefix for the references path',
    )
    # --gcp-project
    parser.add_argument(
        '--gcp-project',
        required=True,
        help='GCP project',
    )

    parser.add_argument(
        'name',
        help='Name of the reference source to transfer',
    )
    return parser.parse_args(args)


def main(name: str, references_prefix: str, gcp_project: str) -> None:
    """Copy the entries of one source that are missing from the bucket.

    Entries already present are left alone, so a listed entry that has since
    vanished upstream is harmless while our copy exists, and reference data is
    never overwritten in place.

    Args:
        name (str): Name of the reference source to transfer
        references_prefix (str): Prefix for the references path
        gcp_project (str): GCP project
    """
    source = {s.name: s for s in SOURCES}[name]
    if source.transfer_cmd and source.src:
        for transfer_cmd, src, dst in source.transfers(references_prefix):
            if gcs_file_exists(dst):
                print(f'{dst} exists, skipping')
                continue
            cmd = transfer_cmd(src=src, dst=dst, project=gcp_project)
            print(cmd)
            # bash for `set -o pipefail` in the curl commands; check so a failed copy
            # fails the job instead of deploying a config that points at it.
            subprocess.run(cmd, shell=True, check=True, executable='/bin/bash')


if __name__ == '__main__':
    args = parser(sys.argv[1:])
    main(
        name=args.name,
        references_prefix=args.references_prefix,
        gcp_project=args.gcp_project,
    )
