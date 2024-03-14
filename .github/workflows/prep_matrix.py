"""
Prepare test matrix (to transfer references in parallel)
"""

import argparse
import sys
from os.path import join

from cloudpathlib import AnyPath

from references import SOURCES as NEW_SOURCES

try:
    # copied into place by the github action
    from references_before import SOURCES as OLD_SOURCES
except ImportError:
    OLD_SOURCES = []


def generate_matrix(references_prefix: str) -> dict:
    """Generate matrix for transferring references in parallel

    Args:
        references_prefix (str): References prefix

    Returns:
        dict: {"include": [<list of transfers>]}
    """
    transfers = {}
    for source in NEW_SOURCES:
        old_sources_d = {s.name: s for s in OLD_SOURCES}
        dst_path = join(references_prefix, source.dst)
        is_changed = (
            source.name not in old_sources_d
            or source.src != old_sources_d[source.name].src
            or source.dst != old_sources_d[source.name].dst
            or not AnyPath(dst_path).exists()
        )
        if is_changed and source.src and source.transfer_cmd:
            print(f'{source.name} has changed, will transfer', file=sys.stderr)
            transfers[source.name] = {'src': source.src, 'dst': dst_path}
        else:
            print(
                f'{source.name} has not changed since previous revision',
                file=sys.stderr,
            )

    if not transfers:
        return {}
    return {
        'include': [
            {
                'name': name,
                'src': data['src'],
                'dst': data['dst'],
            }
            for name, data in transfers.items()
        ]
    }


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--references-prefix', help='Prefix for references')
    return parser.parse_args()


def print_matrix(matrix: dict):
    print(str(matrix).replace(' ', ''), end='', file=sys.stderr)
    print(str(matrix).replace(' ', ''), end='')


if __name__ == '__main__':
    args = parse_args()
    matrix = generate_matrix(references_prefix=args.references_prefix)
    print_matrix(matrix)
