import sys,os,multiprocessing
from collections import defaultdict

import pysam, vcf
from sberry.utils import *

#import matplotlib.pyplot as plt
#import seaborn as sns
#import pandas as pd

class CondRead:
    def __init__(self):
        self.plist = []
        self.seq = []

class PhaseSet:
    def __init__(self):
        self.id = None
        self.plist = []
        self.refseq = []
        self.haplo = defaultdict(list)
        self.creads = defaultdict(CondRead)

def hamming_ratio(str_a,str_b):
    return sum(a!=b for a,b in zip(str_a,str_b))/len(str_a)

def average_contig_hratio(bamfile,ctg,ps_list):
    with pysam.AlignmentFile(bamfile,'rb') as bam:
        for phaseset in ps_list:
            ps_start=phaseset.plist[0]-1 # pysam pileup require 0-based indices
            ps_end=phaseset.plist[-1]
            ps_positions=set(phaseset.plist)
            for pileupcol in bam.pileup(reference=ctg, start=ps_start, stop=ps_end, min_base_quality=20, ignore_overlaps=False, stepper='all'):
                pos=pileupcol.reference_pos
                if pos not in ps_positions: # FIXME: this could be done more efficiently
                    continue
                for pileupread in pileupcol.pileups:
                    alignment=pileupread.alignment
                    nuc = alignment.query_sequence[pileupread.query_position] if not pileupread.is_del and not pileupread.is_refskip else '-'
                    read_name=alignment.query_name
                    read_overlap=min(ps_end,alignment.reference_end)-max(ps_start,alignment.reference_start)
                    assert read_overlap > 0
                    if read_overlap >= min(3000,ps_end-ps_start):
                        if len(phaseset.creads[read_name].plist) == 0 or phaseset.creads[read_name].plist[-1] != pos:
                            phaseset.creads[read_name].plist.append(pos)
                            phaseset.creads[read_name].seq.append(nuc)
    ctg_hratio=defaultdict(list)
    for phaseset in ps_list:
        creads=phaseset.creads
        for read_name,cread in creads.items():
            cread_positions = set(cread.plist)
            assert len(cread_positions) == len(cread.plist)
            ps_indices = [ i for i,pos in enumerate(phaseset.plist) if pos in cread_positions ]
            ps_haplotypes = { h:''.join(htype[i] for i in ps_indices) for h,htype in phaseset.haplo.items() }
            cread_seq = ''.join(cread.seq)
            hratio = min(hamming_ratio(cread_seq,hap_seq) for hap_seq in ps_haplotypes.values())
            ctg_hratio[phaseset.id].append(hratio)
    return (ctg,ctg_hratio)


def average_hratio(bamfile,vcffile,nproc=1):

    # load SNV dictionary of called positions
    print_status(f'Loading SNVs')
    psDict=defaultdict(lambda:defaultdict(PhaseSet))
    with open(vcffile,'rb') as fp:
        vcf_reader=vcf.Reader(fp)
        for rec in vcf_reader:
            call=rec.samples[0]
            if call.phased:
                phaseset = psDict[rec.CHROM][call.data.PS]
                phaseset.id = f'{call.data.PS}'
                phaseset.plist.append(rec.POS-1)
                phaseset.refseq.append(rec.REF)
                gtypes = call.gt_bases.split('|')
                for i,gt in enumerate(gtypes):
                    phaseset.haplo[i].append(gt)
    contigs = list(psDict.keys())

    results=[]
    def get_result(res):
        ctg,ctg_hratio=res
        for ps_id,ps_list in ctg_hratio.items():
            nreads = len(ps_list)
            avg_hratio = sum(ps_list)/nreads if nreads > 0 else 0.0
            results.append( (ctg,ps_id,avg_hratio,nreads) )

    print_status(f'Processing {sum(len(ps) for ps in psDict.values())} phased regions')
    pool = multiprocessing.Pool(processes=nproc,maxtasksperchild=1)
    for ctg in contigs:
        phasesets = list(psDict[ctg].values())
        pool.apply_async(average_contig_hratio, args=(bamfile,ctg,phasesets,), callback=get_result)
    pool.close()
    pool.join()

    tot_nreads=0
    weighted_sum=0
    for contig,psid,hratio,nreads in results:
        weighted_sum+=hratio*nreads
        tot_nreads+=nreads

    return (weighted_sum/nreads,nreads)
    
#    print(f'Writing results to {opt.prefix}.hratio.tsv',file=sys.stderr)
#    with open(f'{opt.prefix}.hratio.tsv','w') as out:
#        out.write('#contig\tPS\tavg-hratio\tnreads\n')
#        for rec in results:
#            out.write('\t'.join(rec))
#            out.write('\n')

#    df = pd.DataFrame(data={'type':'best_haplo','hratio':best_hratios})
#    #df = df.append(pd.DataFrame(data={'type':'reference','hratio':ref_hratios}))
#    
#    sns.set_theme('paper')
#    sns.displot(df, x='hratio', hue='type', fill=True, legend=False, kind='kde', cut=0)
#    plt.axvline(x=statistics.mean(best_hratios),color=sns.color_palette()[0])
#    plt.axvline(x=statistics.median(best_hratios),color=sns.color_palette()[1])
#    plt.title(f'Kernel density plot of best Hamming ratios ({opt.prefix})')
#    plt.xlabel('Hamming ratio')
#    plt.ylabel('Density')
#    plt.tight_layout()
#    plt.savefig(f'{opt.prefix}.hratios.pdf')

