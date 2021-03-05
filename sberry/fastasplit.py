#!/usr/bin/env python3

import sys,re
from Bio import SeqIO

def insert_newlines(inString, every=80):
    return inString if every <= 0 else '\n'.join(inString[i:i+every] for i in range(0, len(inString), every))


#usage statement
if len(sys.argv) != 3:
    print("Usage: python3 fastasplit.py <input.fasta> <output_path>")
    sys.exit()


for sequence in SeqIO.parse(sys.argv[1], "fasta"):
    outDir = sys.argv[2].rstrip("/")
    fileName = re.sub(r'[\\/:"*?<>|]+', "", str(sequence.id))
    pathName = f'{outDir}/{fileName}.fa'
    with open(pathName, 'w') as myFile:
        myFile.write(f'>{sequence.id}\n')
        myFile.write(f'{insert_newlines(str(sequence.seq))}\n')

