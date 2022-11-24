import os
import sys
import subprocess
from references import SOURCES

REFERENCES_PREFIX = os.environ['REFERENCES_PREFIX']
NAME = sys.argv[1]

source = {s.name: s for s in SOURCES}[NAME]
if source.transfer_cmd:
    cmd = source.transfer_cmd(
        src=source.src,
        dst=os.path.join(REFERENCES_PREFIX, source.dst),
    )
    print(cmd)
    subprocess.run(cmd, shell=True)
