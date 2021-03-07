from typing import List
from collections import defaultdict,namedtuple

import vcf


class Phaseset:

    def __init__(self):
        self.reference = None
        self.psid = None
        self.positions = []
        self.haplo = defaultdict(list)
    
    def start(self):
        return self.positions[0]

    def end(self):
        return self.positions[-1]+1

    def density(self):
        return len(self.positions)/(self.end()-self.start())

    def __repr__(self):
        return f'Phaseset( ref={self.reference}, id={self.psid}, |positions|={len(self.positions)}, density={self.density():.4f} )'


class PhasedPosition:

    def __init__(self,record):
        call = record.samples[0]
        self.reference = record.CHROM
        self.psid = f'{call.data.PS}'
        self.position = record.POS-1
        self.gtypes = call.gt_bases.split('|')


class PhasesetCollection:

    def __init__(self):
        self.psdict=defaultdict(lambda:defaultdict(Phaseset))

    def load_from_vcf(self,vcffile,min_density=0.0):
        # load phased positions from the input VCF file
        phased_positions = []
        with open(vcffile,'rb') as fp:
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
        # remove phaset with too few SNVs
        for ref,psd in self.psdict.items():
            for ps_id in [ ps_id for ps_id,phaseset in psd.items() if 100.0 * phaseset.density() < min_density ]:
                del psd[ps_id]

    def size(self):
        return sum(len(ps) for ps in self.psdict.values())

    def references(self):
        return list(self.psdict.keys())

    def phasesets(self,reference:str) -> List[Phaseset]:
        return list(self.psdict.get(reference).values())

    def phaseset(self,reference:str,psid:str) -> Phaseset:
        return self.psdict[reference].get(psid) if reference in self.psdict else None

