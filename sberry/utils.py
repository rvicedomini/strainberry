import sys
from datetime import datetime


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

