import sys
from datetime import datetime


rctable = str.maketrans("ACTGNactgn", "TGACNtgacn")
def reverse_complement(seq):
    return seq.translate(rctable)[::-1]


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

