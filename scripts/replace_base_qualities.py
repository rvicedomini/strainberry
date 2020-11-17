#!/usr/bin/env python3

import sys,os,argparse
import pysam

def main(argv=None):

    parser = argparse.ArgumentParser()
    parser.add_argument('-b','--bam', dest='bamFile', required=True, help='input BAM file (use - for reading from stdin)')
    parser.add_argument('-o','--output', dest='outFile', required=True, help='output BAM files')
    parser.add_argument('-t','--threads', dest='threads', type=int, default=int(4), help='number of threads for compressing/decompressing BAM files')
    opt = parser.parse_args()

    validParameters=True
    if opt.bamFile != "-" and not os.path.isfile(opt.bamFile):
        validParameters=False
        print(f'-b|--bam file "{opt.bamFile}" does not exist.',file=sys.stderr)
    if not validParameters:
        return 1

    with pysam.AlignmentFile(opt.bamFile,'rb',threads=opt.threads) as inBam:
        with pysam.AlignmentFile(opt.outFile,'wb',template=inBam,threads=opt.threads) as outBam:
            for read in inBam:
                if read.query_sequence is not None:
                    read.query_qualities=[20]*read.query_length
                outBam.write(read)

    return 0

# Check if the program is not being imported
if __name__ == "__main__":
    sys.exit(main())

