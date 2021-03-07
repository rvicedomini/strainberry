import sys,os,multiprocessing
from collections import defaultdict

import pysam
from sberry.utils import *
from sberry.phaseset import PhasesetCollection

#import matplotlib.pyplot as plt
#import seaborn as sns
#import pandas as pd


class CondensedRead:
    def __init__(self):
        self.positions = []
        self.seq = []


def hamming_ratio(str_a,str_b):
    assert len(str_a)==len(str_b) and len(str_a)>0
    return sum(a!=b for a,b in zip(str_a,str_b))/len(str_a)


def average_contig_hratio( bamfile, ctg, ps_list, min_read_overlap=3000):

    with pysam.AlignmentFile(bamfile,'rb') as bam:
        for phaseset in ps_list:
            ps_start=phaseset.start()
            ps_end=phaseset.end()
            ps_positions=set(phaseset.positions)
            phaseset.creads=defaultdict(CondensedRead)
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
                    if read_overlap >= min(min_read_overlap,ps_end-ps_start):
                        if len(phaseset.creads[read_name].positions) == 0 or phaseset.creads[read_name].positions[-1] != pos:
                            phaseset.creads[read_name].positions.append(pos)
                            phaseset.creads[read_name].seq.append(nuc)
    
    ctg_hratio=defaultdict(list)
    for phaseset in ps_list:
        creads=phaseset.creads
        for read_name,cread in creads.items():
            cread_positions = set(cread.positions)
            assert len(cread_positions) == len(cread.positions)
            ps_indices = [ i for i,pos in enumerate(phaseset.positions) if pos in cread_positions ]
            ps_haplotypes = { h:''.join(htype[i] for i in ps_indices) for h,htype in phaseset.haplo.items() }
            cread_seq = ''.join(cread.seq)
            hratio = min(hamming_ratio(cread_seq,hap_seq) for hap_seq in ps_haplotypes.values())
            ctg_hratio[phaseset.psid].append(hratio)
    
    return (ctg,ctg_hratio)


def average_hratio( bamfile:str, psc:PhasesetCollection, 
        snv_dens:float = None, nproc: int = 1 ):

    references = psc.references()

    results=[]
    def get_result(res):
        ctg,ctg_hratio=res
        for ps_id,hr_list in ctg_hratio.items():
            nreads = len(hr_list)
            avg_hratio = sum(hr_list)/nreads if nreads > 0 else 0.0
            results.append( (ctg,ps_id,avg_hratio,nreads) )
    
    pool = multiprocessing.Pool(processes=nproc,maxtasksperchild=1)
    for ctg in references:
        ctg_phasesets = psc.phasesets(ctg)
        pool.apply_async(average_contig_hratio, args=(bamfile,ctg,ctg_phasesets,), callback=get_result)
    pool.close()
    pool.join()

    tot_nreads=0
    weighted_sum=0
    for contig,psid,hratio,nreads in results:
        weighted_sum+=hratio*nreads
        tot_nreads+=nreads
    
    num_phasesets = psc.size()
    avg_hratio = weighted_sum/tot_nreads if tot_nreads > 0 else 1.0
    return (avg_hratio,num_phasesets,tot_nreads)
    
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

