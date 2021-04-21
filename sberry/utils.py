import sys
from datetime import datetime


rctable = str.maketrans("ACTGNactgn", "TGACNtgacn")
def reverse_complement(seq):
    return seq.translate(rctable)[::-1]


# input is the path to a fasta index 
# (i.e., .fai file generated with samtools faidx)
# output is a tuple -> assembly_size, sequence_number, assembly_N50
def faidx_stats(fname):
    seqlen_list = []
    with open(fname,'r') as fai:
        for line in fai:
            _,seqlen,_=line.split('\t',2)
            seqlen_list.append(int(seqlen))
    seqlen_list.sort(reverse=True)
    
    assembly_size = sum(seqlen_list)
    assembly_nseq = len(seqlen_list)
    assembly_largest = seqlen_list[0] if assembly_nseq>0 else 0
    assembly_n50 = 0
    partial_sum = 0
    for seqlen in seqlen_list:
        partial_sum += seqlen
        if 2*partial_sum >= assembly_size:
            assembly_n50 = seqlen
            break
    return assembly_size,assembly_nseq,assembly_n50
    

def insert_newlines(inString, every=80):
    return inString if every <= 0 else '\n'.join(inString[i:i+every] for i in range(0, len(inString), every))


def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)

def datetime_now():
    return datetime.now().strftime("%Y-%m-%d %H:%M:%S")

def print_error(msg):
    eprint(f'[{datetime_now()}] error: {msg}')

def print_warning(msg):
    eprint(f'[{datetime_now()}] warning: {msg}')

def print_status(msg):
    eprint(f'[{datetime_now()}] {msg}')

