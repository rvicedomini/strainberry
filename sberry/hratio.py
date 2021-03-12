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

class AverageHammingRatio:
    
    def __init__(self):
        self.ahrdict = {}

    def _hamming_ratio(self,str_a,str_b):
        assert len(str_a)==len(str_b) and len(str_a)>0
        return sum(a!=b for a,b in zip(str_a,str_b))/len(str_a)

    def _average_contig_hratio(self, bamfile, ctg, ps_list, min_read_overlap=3000):
        with pysam.AlignmentFile(bamfile,'rb') as bam:
            for phaseset in ps_list:
                ps_start=phaseset.start()
                ps_end=phaseset.end()
                ps_positions=set(phaseset.positions)
                phaseset.creads=defaultdict(CondensedRead)
                for pileupcol in bam.pileup(reference=ctg, start=ps_start, stop=ps_end, min_base_quality=0, ignore_overlaps=False, stepper='all'):
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
        ps_hratios=[]
        for phaseset in ps_list:
            ps_hratio_sum=0
            ps_hratio_count=0
            creads=phaseset.creads
            for read_name,cread in creads.items():
                cread_positions = set(cread.positions)
                assert len(cread_positions) == len(cread.positions)
                ps_indices = [ i for i,pos in enumerate(phaseset.positions) if pos in cread_positions ]
                ps_haplotypes = { h:''.join(htype[i] for i in ps_indices) for h,htype in phaseset.haplo.items() }
                cread_seq = ''.join(cread.seq)
                ps_hratio_sum += min(self._hamming_ratio(cread_seq,hap_seq) for hap_seq in ps_haplotypes.values())
                ps_hratio_count += 1
            ps_hratio_avg = ps_hratio_sum/ps_hratio_count if ps_hratio_count > 0 else 1.0
            ps_hratios.append((phaseset.psid,ps_hratio_avg,ps_hratio_count))
        return (ctg,ps_hratios)

    def _aggregate_contig_hratios(self,result):
        ctg,ctg_ps_hratios=result
        self.ahrdict[ctg]={}
        for psid,ps_hratio_avg,ps_hratio_count in ctg_ps_hratios:
            self.ahrdict[ctg][psid]=(ps_hratio_avg,ps_hratio_count)

    def reference_average_hratio(self,reference):
        ref_hratio_avg=1.0
        ref_hratio_count=0
        if reference in self.ahrdict:
            for ps_hratio,ps_count in self.ahrdict[reference].values():
                if ps_count == 0:
                    continue
                ref_hratio_avg = (ref_hratio_avg*ref_hratio_count + ps_hratio*ps_count)/(ref_hratio_count+ps_count)
                ref_hratio_count = ref_hratio_count + ps_count
        return (ref_hratio_avg,ref_hratio_count)

    def phaseset_average_hratio(self,reference,psid):
        if reference in self.ahrdict and psid in self.ahrdict[reference]:
            return self.ahrdict[reference][psid]
        return (1.0,0)

    def compute( self, bamfile:str, psc:PhasesetCollection, nproc:int = 1 ):

        self.ahrdict = {}
        references = psc.references()
        
        pool = multiprocessing.Pool(processes=nproc,maxtasksperchild=1)
        for ctg in references:
            ctg_phasesets = psc.phasesets(ctg)
            pool.apply_async(self._average_contig_hratio, args=(bamfile,ctg,ctg_phasesets,), callback=self._aggregate_contig_hratios)
        pool.close()
        pool.join()

        hratio_avg=1.0
        hratio_size=0
        for reference in self.ahrdict:
            ref_hratio,ref_size=self.reference_average_hratio(reference)
            if ref_size == 0:
                continue
            hratio_avg = (hratio_avg*hratio_size + ref_hratio*ref_size)/(hratio_size+ref_size)
            hratio_size = hratio_size + ref_size
        
        num_phasesets = psc.size()
        return (hratio_avg,num_phasesets,hratio_size)
    
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

