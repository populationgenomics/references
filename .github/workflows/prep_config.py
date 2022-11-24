"""
Prepare ready contig TOML
"""

import os
import toml
from references import SOURCES, GENOME_BUILD

REFERENCES_PREFIX = os.environ['REFERENCES_PREFIX']

d = {'genome_build': GENOME_BUILD}

for source in SOURCES:
    dst_path = os.path.join(REFERENCES_PREFIX, source.dst)
    if not source.files:
        d[source.name] = str(dst_path)
    else:
        d[source.name] = {
            k: os.path.join(dst_path, suffix) for k, suffix in source.files.items()
        }

print(toml.dumps({'references': d}))
