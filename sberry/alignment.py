from dataclasses import dataclass,field
from enum import IntEnum


# DOVETAIL_SUFFIX < REF_CONTAINED < DOVETAIL_PREFIX is important for adding edges
class MappingType(IntEnum):
    DOVETAIL_SUFFIX=0
    QUERY_PREFIX=1
    QUERY_SUFFIX=2
    REF_PREFIX=3
    REF_SUFFIX=4
    INTERNAL=5
    QUERY_CONTAINED=6
    REF_CONTAINED=7
    DOVETAIL_PREFIX=8


def classify_mapping(seq1,seq2,o=50,r=0.05):
    b1,e1,l1 = seq1
    b2,e2,l2 = seq2
    left_overhang = min(b1,b2)
    right_overhang = min(l1-e1,l2-e2)
    maplen = max(e1-b1,e2-b2)
    oh_threshold = min(o,maplen*r)
    if b1 <= oh_threshold and right_overhang > oh_threshold: 
        return MappingType.QUERY_PREFIX
    elif b2 <= oh_threshold and right_overhang > oh_threshold:
        return MappingType.REF_PREFIX
    elif left_overhang > oh_threshold and l1-e1 <= oh_threshold:
        return MappingType.QUERY_SUFFIX
    elif left_overhang > oh_threshold and l2-e2 <= oh_threshold:
        return MappingType.REF_SUFFIX
    elif left_overhang > oh_threshold or right_overhang > oh_threshold:
        return MappingType.INTERNAL
    elif b1<=b2 and l1-e1<=l2-e2:
        return MappingType.QUERY_CONTAINED
    elif b1>=b2 and l1-e1>=l2-e2:
        return MappingType.REF_CONTAINED
    elif b1<=b2:
        return MappingType.DOVETAIL_PREFIX
    return MappingType.DOVETAIL_SUFFIX


def overlap_length(arange,brange):
    abeg,aend,alen=arange
    bbeg,bend,blen=brange
    if abeg > bbeg:
        return overlap_length(brange,arange)
    obeg = max(abeg,bbeg)
    oend = min(aend,bend)
    return oend-obeg


@dataclass
class ReadAlignment:
    query: str
    query_length: int
    query_start: int
    query_end: int
    strand: str
    reference: str
    reference_length: int
    reference_start: int
    reference_end: int
    matches: int
    mapping_length: int
    mapq: int
    tags: dict = field(repr=False)

    def identity(self):
        return (100.0*self.matches)/self.mapping_length if self.mapping_length!=0 else 0.0

    def has_tag(self,tag):
        return tag in self.tags

    def get_tag(self,tag):
        return self.tags[tag] if self.has_tag(tag) else None


def parse_paf(infile):
    for line in infile:
        line=line.strip()
        if not line or line.startswith('#'):
            continue
        col=line.split('\t')
        qid=col[0]
        qlen,qbeg,qend=[int(x) for x in col[1:4]]
        strand=col[4]
        tid=col[5] # contig_26_710192_h1-1:11849-328803
        tlen,tbeg,tend=[int(x) for x in col[6:9]]
        nmatches,maplen,mapq=[int(x) for x in col[9:12]]
        tags = { k:v for k,v in map(lambda s:s.split(':',2)[::2],col[12:]) }
        yield ReadAlignment(qid,qlen,qbeg,qend,strand,tid,tlen,tbeg,tend,nmatches,maplen,mapq,tags)


