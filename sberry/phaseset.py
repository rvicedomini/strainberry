from typing import List
from collections import defaultdict

import vcf


class Phaseset:

    def __init__(self):
        self.reference = None
        self.psid = None
        self.positions = []
        self.haplo = defaultdict(list)
        self.nreads_tagged = 0 # number of tagged primary alignments covering the phaseset
        self.nreads_mapped = 0 # number of primary alignments covering the phaseset

    def id(self):
        return self.psid
    
    def start(self):
        return self.positions[0]

    def end(self):
        return self.positions[-1]+1

    def density(self):
        return len(self.positions)/(self.end()-self.start())

    def __len__(self):
        return self.end()-self.start() if len(self.positions)>0 else 0

    def __repr__(self):
        return f'Phaseset( ref={self.reference}, id={self.psid}, |positions|={len(self.positions)}, density={self.density():.4f} tagged={self.nreads_tagged} mapped={self.nreads_mapped} )'


class PhasedPosition:

    def __init__(self,record):
        call = record.samples[0]
        self.reference = record.CHROM
        self.psid = f'{call.data.PS}'
        self.position = record.POS-1
        self.gtypes = call.gt_bases.split('|')


def _PhasesetDict():
    return defaultdict(Phaseset)


class PhasesetCollection:

    def __init__(self):
        self.psdict=defaultdict(_PhasesetDict)

    def __contains__(self, contig):
        return (contig in self.psdict)

    def load_from_vcf(self,vcffile,min_density=0.0,min_length=0):
        # load phased positions from the input VCF file
        phased_positions = []
        with open(vcffile,'r') as fp:
            vcf_reader=vcf.Reader(fp)
            phased_positions = sorted( (PhasedPosition(rec) for rec in vcf_reader if rec.samples[0].phased),
                    key=lambda pp:(pp.reference,pp.psid,pp.position) )
        # load internal dictionary
        for pp in phased_positions:
            phaseset=self.psdict[pp.reference][pp.psid]
            phaseset.reference=pp.reference
            phaseset.psid=pp.psid
            phaseset.positions.append(pp.position)
            for i,gt in enumerate(pp.gtypes):
                phaseset.haplo[i].append(gt)
        # remove phasets with too few SNVs
        for ref,psd in self.psdict.items():
            for ps_id in [ps_id for ps_id,phaseset in psd.items() if len(phaseset) < min_length or 100.0 * phaseset.density() < min_density]:
                del psd[ps_id]
        # remove references with no remaining phasesets
        for ref in [ref for ref,psd in self.psdict.items() if len(psd)==0]:
            del self.psdict[ref]

    def update(self,other):
        return self.psdict.update(other.psdict)

    def remove_reference(self,reference):
        del self.psdict[reference]

    def size(self):
        return sum(len(ps) for ps in self.psdict.values())

    def references(self):
        return list(self.psdict.keys())

    def phasesets(self,reference:str) -> List[Phaseset]:
        return list(self.psdict[reference].values())

    def phaseset(self,reference:str,psid:str) -> Phaseset:
        return self.psdict[reference].get(psid) if reference in self.psdict else None


