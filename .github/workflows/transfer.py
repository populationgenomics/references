"""
Transfer one reference source
"""

import argparse
import os
import subprocess
import sys

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
    """Main function for the script

    Args:
        name (str): Name of the reference source to transfer
        references_prefix (str): Prefix for the references path
        gcp_project (str): GCP project
    """
    source = {s.name: s for s in SOURCES}[name]
    if source.transfer_cmd and source.src:
        cmd = source.transfer_cmd(
            src=source.src,
            dst=os.path.join(references_prefix, source.dst),
            project=gcp_project,
        )
        print(cmd)
        subprocess.run(cmd, shell=True)


if __name__ == '__main__':
    args = parser(sys.argv[1:])
    main(
        name=args.name,
        references_prefix=args.references_prefix,
        gcp_project=args.gcp_project,
    )
