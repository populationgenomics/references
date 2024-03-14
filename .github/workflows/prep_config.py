"""
Prepare ready contig TOML
"""

import os
import toml
import argparse

from references import SOURCES, GENOME_BUILD

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--references-prefix', type=str)
    return parser.parse_args()

def main(references_prefix: str) -> dict
    d: dict[str, str | dict[str, str]] = {'genome_build': GENOME_BUILD}

    for source in SOURCES:
        dst_path = os.path.join(references_prefix, source.dst)
        if not source.files:
            d[source.name] = str(dst_path)
        else:
            d[source.name] = {
                k: os.path.join(dst_path, suffix) for k, suffix in source.files.items()
            }

    return {'references': d}

if __name__ == '__main__':
    args = parse_args()
    print(toml.dumps(main(references_prefix=args.references_prefix)))
