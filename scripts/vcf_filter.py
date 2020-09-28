#!/usr/bin/env python3

import sys,argparse

def main( argv = None ):
    parser=argparse.ArgumentParser()
    parser.add_argument('-q','--qual', dest='minQual', type=float, default=50, help='keep only variants above this QUAL threshold')
    opt=parser.parse_args()
    for line in sys.stdin:
        if line.startswith('#'):
            sys.stdout.write(line)
            continue
        rec=line.rstrip().split('\t') # CHROM  POS  ID  REF  ALT  QUAL  FILTER  INFO  FORMAT  SAMPLE
        qual=float(rec[5])
        fltr=rec[6]
        if qual < opt.minQual or fltr == 'dn':
            continue
        frmt=rec[8].split(':')
        sample=rec[9].split(':')
        for i,val in enumerate(frmt):
            # round GQ/UQ value so that whatshap does not complain
            if val == 'GQ' or val == 'UQ':
                sample[i]=f'{int(float(sample[i]))}'
        sys.stdout.write('\t'.join(rec[0:9])+'\t'+':'.join(sample)+'\n')
    return 0

# Check if the program is not being imported
if __name__ == "__main__":
    sys.exit(main())
