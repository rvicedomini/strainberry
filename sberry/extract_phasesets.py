#!/usr/bin/env python3

import sys, errno, os, argparse, gzip, re

from collections import defaultdict
from Bio import SeqIO


def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)

def insert_newlines(inString, every=80):
    return inString if every <= 0 else '\n'.join(inString[i:i+every] for i in range(0, len(inString), every))

def mask_sequence(mseq,beg,end):
    for i in range(beg,end):
        mseq[i]='N'
    return mseq

def overlap_length(a,b):
    return min(a[1],b[1])-max(a[0],b[0])


def main( argv = None ):

    # GET PARAMETERS
    parser = argparse.ArgumentParser()
    parser.add_argument('-i','--vcf', dest='vcfFile', required=True, help='phased VCF file produced by whatshap polyphase')
    parser.add_argument('-f','--fasta', dest='fastaFile', required=True, help='FASTA file of the reference sequences')
    parser.add_argument('-l','--list', dest='listFile', help='list of phasesets to take into account')
    parser.add_argument('-m','--min-length', dest='minLength', type=int, default=0, help='consider only sequences longer than this threshold')
    parser.add_argument('-p','--prefix', dest='outPrefix', default='out', help='prefix of output files')
    parser.add_argument('--complement', dest='doComplement', action='store_true', help='extract non-phased regions')
    options = parser.parse_args()

    # VALIDATE PARAMETERS
    validParameters = True
    if not os.path.isfile(options.vcfFile):
        validParameters = False
        eprint( "VCF file \"{}\" does not exist.".format(options.vcfFile))
    if not os.path.isfile(options.fastaFile):
        validParameters = False
        eprint("FASTA file \"{}\" does not exist.".format(options.fastaFile))
    if (options.listFile is not None) and (not os.path.isfile(options.listFile)):
        validParameters = False
        eprint(f'List file "{options.listFile}" does not exist.')
    if not validParameters:
        return 1

    # Read gzipped VCF file
    varDict = defaultdict(dict)
    varDens = defaultdict(dict)
    with gzip.open(options.vcfFile,'rt') as vcfFile:
        for line in vcfFile:
            if line.startswith('#'):
                continue
            #CHROM  POS  ID  REF  ALT  QUAL  FILTER  INFO  FORMAT  unknown
            rec = line.rstrip().split('\t')
            chrom = rec[0]
            pos   = int(rec[1])-1
            vinfo  = { k:v for k,v in zip(rec[8].split(':'),rec[9].split(':')) }
            if ('PS' in vinfo) and (vinfo['PS']!='.'):
                pid = vinfo['PS']
                psDict = varDict[chrom]
                psDens = varDens[chrom]
                if pid in psDict:
                    psBeg = min(pos,psDict[pid][0])
                    psEnd = max(pos+len(rec[3]),psDict[pid][1])
                    psDict[pid] = (psBeg,psEnd)
                    psDens[pid] += 1
                else:
                    psDict[pid] = (pos,pos+len(rec[3]))
                    psDens[pid] = 1

    psSubDict = defaultdict(set)
    if options.listFile is not None:
        with open(options.listFile,'r') as listFile:
            for line in listFile:
                psid = line.rstrip().split('\t')[0]
                chrom,pid = psid.rsplit('_',1)
                psSubDict[chrom].add(pid)
    else:
        for chrom in varDict:
            for pid in varDict[chrom]:
                psSubDict[chrom].add(pid)

    with open(options.fastaFile,'r') as fastaFile:
        with open(f'{options.outPrefix}.fa','w') as outFasta, open(f'{options.outPrefix}.len.tsv','w') as outLenTsv, open(f'{options.outPrefix}.cov.tsv','w') as covFile, open(f'{options.outPrefix}.dens.tsv','w') as densFile:
            for rec in SeqIO.parse(fastaFile,'fasta'):
                psDict = varDict[rec.id]
                psDens = varDens[rec.id]
                if options.doComplement:
                    masked = rec.seq.tomutable()
                    for ps,reg in psDict.items():
                        if ps in psSubDict[rec.id]:
                            masked = mask_sequence(masked,reg[0],reg[1])
                    maskedseq = str(masked)
                    p = re.compile(r'[^nN]+')
                    sumlen=0
                    for match in p.finditer(maskedseq):
                        seq_beg=match.start()
                        seq_end=match.end()
                        seq_len=seq_end-seq_beg
                        if seq_len >= options.minLength:
                            seq_id=f'{rec.id}_{seq_beg+1}'
                            outFasta.write(f'>{seq_id}  id={rec.id}  region={seq_beg+1}-{seq_end}  len={seq_len}\n')
                            outFasta.write(f'{ insert_newlines(maskedseq[match.start():match.end()]) }\n')
                            outLenTsv.write(f'{seq_id}\t{seq_len}\n')
                            sumlen+=seq_len
                    covFile.write(f'{rec.id}\t{sumlen}\t{len(rec.seq)}\t{float(sumlen)/len(rec.seq):.2f}\n')
                else:
                    sumlen=0
                    for ps,reg in psDict.items():
                        if ps in psSubDict[rec.id]:
                            seq_beg=reg[0]
                            seq_end=reg[1]
                            seq_len=seq_end-seq_beg
                            snp_den=psDens[ps]
                            if seq_len >= options.minLength:
                                seq_id=f'{rec.id}_{ps}'
                                outFasta.write(f'>{seq_id}  id={rec.id}  region={seq_beg+1}-{seq_end}  len={seq_len}\n')
                                outFasta.write(f'{ insert_newlines(str(rec.seq[reg[0]:reg[1]])) }\n')
                                outLenTsv.write(f'{seq_id}\t{seq_len}\n')
                                densFile.write(f'{seq_id}\t{100.0*snp_den/seq_len:.4f}\n')
                                sumlen+=seq_len
                    covFile.write(f'{rec.id}\t{sumlen}\t{len(rec.seq)}\t{float(sumlen)/len(rec.seq):.2f}\n')


# Check if the program is not being imported
if __name__ == "__main__":
    sys.exit(main())
