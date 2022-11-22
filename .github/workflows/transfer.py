import os
import sys
from references import SOURCES

PROJECT = os.environ['PROJECT']
NAME = sys.argv[1]
SRC = sys.argv[2]
DST = sys.argv[3]

source = {s.name: s for s in SOURCES}[NAME]
cmd = None
if source.src.startswith('gs://'):
    if source.files:
        cmd = f'gsutil -u {PROJECT} -m rsync -d -r {SRC} {DST}'
    else:
        cmd = f'gsutil -u {PROJECT} cp {SRC} {DST}'
if source.src.startswith('https://'):
    cmd = (
      f'curl {SRC} -o tmp; gsutil -u {PROJECT} cp tmp {DST}'
    )
if cmd:
    print(cmd)
    os.system(cmd)
