import os
import sys
import subprocess
from references import SOURCES

PROJECT = os.environ['PROJECT']
REFERENCES_PREFIX = os.environ['REFERENCES_PREFIX']
NAME = sys.argv[1]

source = {s.name: s for s in SOURCES}[NAME]
cmd = source.transfer_cmd(PROJECT, REFERENCES_PREFIX)
if cmd:
    print(cmd)
    subprocess.run(cmd, shell=True)
