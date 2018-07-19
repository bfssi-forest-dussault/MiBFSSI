import os
from pathlib import Path

DEPENDENCIES = [
    'ncbi-genome-download',
    'bbmap.sh',
    'bbmerge.sh',
    'sendsketch.sh',
    'fuse.sh',
    'makeblastdb',
    'blastn',
    'qualimap',
    'snippy',
    'nullarbor.pl',
    'nice',
    'spades.py',
    'quast.py'
]

BBDUK_ADAPTERS = Path(os.path.join(os.path.dirname(__file__), 'resources/adapters.fa'))
