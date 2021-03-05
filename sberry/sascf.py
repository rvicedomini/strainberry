#!/usr/bin/env python3
#TODO put code related to paf parsing in a different module

import sys,os,argparse,gzip
import networkx as nx
from operator import itemgetter
from collections import defaultdict
from dataclasses import dataclass,field
from statistics import median
from enum import Enum
from Bio import SeqIO


def eprint(*args, **kwargs):
   print(*args, file=sys.stderr, **kwargs)

# TODO: possibly provide a better implementation
rctable = str.maketrans("ACTGNactgn", "TGACNtgacn")
def reverse_complement(seq):
    return seq.translate(rctable)[::-1]

def insert_newlines(inString, every=80):
    return inString if every <= 0 else '\n'.join(inString[i:i+every] for i in range(0, len(inString), every))


# DOVETAIL_SUFFIX < REF_CONTAINED < DOVETAIL_PREFIX is important for adding edges
class MappingType(Enum):
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


def fix_neighbors(asmGraph):
    for a in asmGraph.nodes:
        toremove = set( b for b in asmGraph.nodes[a]['inedges'] if not asmGraph.has_edge(a,b) )
        asmGraph.nodes[a]['inedges'] -= toremove
        toremove = set( b for b in asmGraph.nodes[a]['outedges'] if not asmGraph.has_edge(a,b) )
        asmGraph.nodes[a]['outedges'] -= toremove
    return asmGraph

def remove_simple_transitive_edges(asmGraph,asmContigs,npContigs):
    # TODO: the following is a naive implementation, check for unwanted behaviors...
    # TODO: maybe I can avoid code duplication, just looking at etype edge property
    # remove transitive edges a->c if there exists a->b->c
    removeEdgeSet=set()
    for a,nbrs in asmGraph.adj.items():
        for b in nbrs:
            ab_atype=asmGraph[a][b]['etype'][a]
            if a!=b and asmGraph[a][b]['etype'][b] == MappingType.DOVETAIL_PREFIX:
                for c in asmGraph.nodes[b]['outedges']:
                    bc_ctype=asmGraph[b][c]['etype'][c]
                    if b!=c and asmGraph.has_edge(a,c) and asmGraph[a][c]['etype'][c] == bc_ctype and asmGraph[a][c]['etype'][a] == ab_atype and len(asmGraph.nodes[b]['outedges'])==1 and len(asmGraph.nodes[b]['inedges'])==1: # and c in asmContigs:
                        removeEdgeSet.add((a,c))
                        asmGraph[a][b]['nreads']+=asmGraph[a][c]['nreads']
                        asmGraph[b][c]['nreads']+=asmGraph[a][c]['nreads']
                        #asmGraph[a][c]['ctglink'].add(b)
                        #asmGraph[a][c]['nreads']+=asmGraph[a][b]['nreads']+asmGraph[b][c]['nreads']
                        #asmGraph[a][c]['label']=f'{asmGraph[a][c]["nreads"]}/{len(asmGraph[a][c]["ctglink"])}'
                        #asmGraph.nodes[a]['length']+=asmGraph.nodes[b]['length']
            elif a!=b and asmGraph[a][b]['etype'][b] == MappingType.DOVETAIL_SUFFIX:
                for c in asmGraph.nodes[b]['inedges']:
                    bc_ctype=asmGraph[b][c]['etype'][c]
                    if b!=c and asmGraph.has_edge(a,c) and asmGraph[a][c]['etype'][c] == bc_ctype and asmGraph[a][c]['etype'][a] == ab_atype and len(asmGraph.nodes[b]['outedges'])==1 and len(asmGraph.nodes[b]['inedges'])==1: # and c in asmContigs:
                        removeEdgeSet.add((a,c))
                        asmGraph[a][b]['nreads']+=asmGraph[a][c]['nreads']
                        asmGraph[b][c]['nreads']+=asmGraph[a][c]['nreads']
                        #asmGraph[a][c]['ctglink'].add(b)
                        #asmGraph[a][c]['nreads']+=asmGraph[a][b]['nreads']+asmGraph[b][c]['nreads']
                        #asmGraph[a][c]['label']=f'{asmGraph[a][c]["nreads"]}/{len(asmGraph[a][c]["ctglink"])}'
                        #asmGraph.nodes[a]['length']+=asmGraph.nodes[b]['length']
    asmGraph.remove_edges_from(list(removeEdgeSet))
    return fix_neighbors(asmGraph)

def remove_weak_edges(asmGraph,minReads=10,minReadFrac=0.9):
    removeEdgeSet=set()
    for a in asmGraph.nodes:
        nreads=sum(asmGraph[a][b]['nreads'] for b in asmGraph.nodes[a]['inedges'] if asmGraph.has_edge(a,b))
        for b in asmGraph.nodes[a]['inedges']:
            asmGraph[a][b]['rfrac']=asmGraph[a][b]['nreads']/float(nreads)
            asmGraph[a][b]['label']=f'{asmGraph[a][b]["rfrac"]:.2f}/{asmGraph[a][b]["nreads"]}/{len(asmGraph[a][b]["ctglink"])}/{int(median(asmGraph[a][b]["overlaps"]))}'
            if asmGraph[a][b]['rfrac'] < minReadFrac or nreads < minReads:
                removeEdgeSet.add((a,b))
        nreads=sum(asmGraph[a][b]['nreads'] for b in asmGraph.nodes[a]['outedges'] if asmGraph.has_edge(a,b))
        for b in asmGraph.nodes[a]['outedges']:
            asmGraph[a][b]['rfrac']=asmGraph[a][b]['nreads']/float(nreads)
            asmGraph[a][b]['label']=f'{asmGraph[a][b]["rfrac"]:.2f}/{asmGraph[a][b]["nreads"]}/{len(asmGraph[a][b]["ctglink"])}/{int(median(asmGraph[a][b]["overlaps"]))}'
            if asmGraph[a][b]['rfrac'] < minReadFrac or nreads < minReads:
                removeEdgeSet.add((a,b))
    asmGraph.remove_edges_from(list(removeEdgeSet))
    return fix_neighbors(asmGraph)
    
def build_scaffolds(asmGraph,ctgSeq):
    visited=set()
    out_scaffolds=[]
    for a in asmGraph.nodes:
        if a not in visited:
            visited.add(a)
            aseq=ctgSeq[a]
            outedges=asmGraph.nodes[a]['outedges']
            if len(outedges)==1 and list(outedges)[0] not in visited:
                for b in asmGraph.nodes[a]['outedges']:
                    edata=asmGraph[a][b]
                    atype=edata['etype'][a]
                    btype=edata['etype'][b]
                    bseq=build_scaffold_from(b,btype,asmGraph,ctgSeq,visited)
                    ab_gap=int(median(edata["overlaps"]))
                    ab_gap=-ab_gap if ab_gap<0 else 0
                    aseq+=(ab_gap*'N') + (bseq if atype!=btype else reverse_complement(bseq))
            inedges=asmGraph.nodes[a]['inedges']
            if len(inedges)==1 and list(inedges)[0] not in visited:
                for b in asmGraph.nodes[a]['inedges']:
                    edata=asmGraph[a][b]
                    atype=edata['etype'][a]
                    btype=edata['etype'][b]
                    bseq=build_scaffold_from(b,btype,asmGraph,ctgSeq,visited)
                    ab_gap=int(median(edata["overlaps"]))
                    ab_gap=-ab_gap if ab_gap<0 else 0
                    aseq=(bseq if atype!=btype else reverse_complement(bseq)) + (ab_gap*'N') + aseq
            out_scaffolds.append(aseq)
    return out_scaffolds

# TODO: AVOID RECURSION => MAKE IT ITERATIVE; also there are maybe too many reverse complement operations
def build_scaffold_from(a,atype,asmGraph,ctgSeq,visited):
    # assuming this function is called on an unvisited node
    visited.add(a)
    aseq=ctgSeq[a]
    if atype == MappingType.DOVETAIL_PREFIX:
        outedges=asmGraph.nodes[a]['outedges']
        if len(outedges) == 1 and list(outedges)[0] not in visited:
            for b in asmGraph.nodes[a]['outedges']:
                edata=asmGraph[a][b]
                atype=edata['etype'][a]
                btype=edata['etype'][b]
                bseq=build_scaffold_from(b,btype,asmGraph,ctgSeq,visited)
                ab_gap=int(median(edata["overlaps"]))
                ab_gap=-ab_gap if ab_gap<0 else 0
                aseq+=(ab_gap*'N') + (bseq if atype!=btype else reverse_complement(bseq))
    elif atype == MappingType.DOVETAIL_SUFFIX:
        inedges=asmGraph.nodes[a]['inedges']
        if len(inedges)==1 and list(inedges)[0] not in visited:
            for b in asmGraph.nodes[a]['inedges']:
                edata=asmGraph[a][b]
                atype=edata['etype'][a]
                btype=edata['etype'][b]
                bseq=build_scaffold_from(b,btype,asmGraph,ctgSeq,visited)
                ab_gap=int(median(edata["overlaps"]))
                ab_gap=-ab_gap if ab_gap<0 else 0
                aseq=(bseq if atype!=btype else reverse_complement(bseq)) + (ab_gap*'N') + aseq
    return aseq


def estimate_n50(asmGraph):
    components = nx.connected_components(asmGraph)
    pathlengths = []
    totsum=0
    for comp in components:
        pl=sum( asmGraph.nodes[n]['length'] for n in comp)
        totsum += pl
        pathlengths.append(pl)
    pathlengths.sort(reverse=True)
    partsum=0
    for pl in pathlengths:
        partsum+=pl
        if partsum >= 0.5*totsum:
            return pl
    return pl


def getRefPos(ctgname):
    if not ctgname.startswith("sberry|"):
        refid,pos=ctgname.rsplit('_',1)
        return (refid,int(pos))
    refid,pos,_,_=ctgname.split('|',1)[1].rsplit('_',3)
    return (refid,int(pos))


def main( argv = None ):

    # GET PARAMETERS
    parser = argparse.ArgumentParser()
    parser.add_argument('--paf', dest='pafFile', required=True, help='gzipped paf file')
    parser.add_argument('--reference', dest='refFile', required=True, help='reference of the mapping')
    parser.add_argument('--hap2ref', dest='hap2ref', help='contig to strain mapping')
    parser.add_argument('--prefix', dest='outPrefix', default='out', help='prefix of output files')
    parser.add_argument('--min-mapq', dest='minMapq', type=int, default=int(40), help='minimum MAPQ value to consider a read alignment during scaffolding')
    parser.add_argument('--min-reads', dest='minReads', type=int, default=int(10), help='minimum number of reads matching the prefix/suffix of a contig')
    parser.add_argument('--min-read-frac', dest='minReadFrac', type=float, default=float(0.75), help='minimum read fraction to define good edges')
    opt = parser.parse_args()

    # VALIDATE PARAMETERS
    validParameters = True
    if not os.path.isfile(opt.pafFile):
        validParameters = False
        eprint(f'--paf file "{opt.pafFile}" does not exist.')
    if not os.path.isfile(opt.refFile):
        validParameters = False
        eprint(f'--reference file "{opt.refFile}" does not exist.')
    if not validParameters:
        return 1

    # retrieve identifiers and put not-phased/assembled contigs in two different sets
    ctgSeq={ ctg.id:str(ctg.seq) for ctg in SeqIO.parse(opt.refFile,'fasta') }
    ctgDict={ ctg_id:len(ctg_seq) for ctg_id,ctg_seq in ctgSeq.items() }
    npContigs={ ctg:ctglen for ctg,ctglen in ctgDict.items() if not ctg.startswith('sberry|') }
    asmContigs={ ctg:ctglen for ctg,ctglen in ctgDict.items() if ctg.startswith('sberry|') }

    psDict=defaultdict(list)
    psAdjSet=set()
    for ctg in asmContigs.keys():
        ps=getRefPos(ctg)
        refid,pos=ps
        psDict[refid].append(ps)
    for refid in psDict:
        psDict[refid].sort()
        psList=psDict[refid]
        for i,ps_i in enumerate(psList):
            ref_i,pos_i=ps_i
            if i>0:
                ps_j=psList[i-1]
                ref_j,pos_j=ps_j
                if ref_i==ref_j:
                    psAdjSet.add( (ps_i,ps_j) )
                    psAdjSet.add( (ps_j,ps_i) )

    fragList=sorted(map(getRefPos,ctgDict.keys()))
    fragAdjSet=set()
    refEnds=defaultdict(set)
    for i,frag in enumerate(fragList):
        ref_i=frag[0]
        if i>0:
            ref_j=fragList[i-1][0]
            if ref_i!=ref_j:
                refEnds[ref_i].add(frag)
                refEnds[ref_j].add(fragList[i-1])
            else:
                fragAdjSet.add( (frag,fragList[i-1]) )
                fragAdjSet.add( (fragList[i-1],frag) )
        if i==len(fragList)-1:
            refEnds[ref_i].add(frag)

    #for ref in refEnds:
    #    print(f'{ref} -> {refEnds[ref]}')

    # create graph and contig->node_id mapping
    asmGraph=nx.Graph(style="filled")
    asmGraph.add_nodes_from(list(npContigs.keys()),ctgtype='notphased')
    asmGraph.add_nodes_from(list(asmContigs.keys()),ctgtype='assembled')
    nx.set_node_attributes(asmGraph,npContigs,name='length')
    nx.set_node_attributes(asmGraph,asmContigs,name='length')
    for n in asmGraph.nodes:
        asmGraph.nodes[n]['outedges']=set()
        asmGraph.nodes[n]['inedges']=set()

    print(f'connected-components: {len(asmGraph.nodes)}')
    print(f'N50-before: {estimate_n50(asmGraph)}')

    # give colors to contigs according to their strain label
    # TODO: right now it only works with <= 12 different strain labels
    if opt.hap2ref is not None:
        ctg2strain={}
        straincol={}
        nstrains=1
        with open(opt.hap2ref,'r') as f:
            for line in f:
                cols=line.rstrip('\n').split('\t')
                n,strain=cols[0:2]
                ctg2strain[n]=strain
                if strain not in straincol:
                    straincol[strain]=nstrains
                    nstrains+=1
            for n in asmGraph.nodes:
                strain = ctg2strain[n] if n in ctg2strain else 'none'
                asmGraph.nodes[n]['style']='filled'
                asmGraph.nodes[n]['fillcolor']=f'/set312/{straincol[strain]}'
    else:
        for n in asmGraph.nodes:
            if n in asmContigs:
                asmGraph.nodes[n]['style']='filled'
                asmGraph.nodes[n]['fillcolor']=f'/set312/1'

    # parse alignments to find dovetail mappings
    read2contigs=defaultdict(list)
    with gzip.open(opt.pafFile,'rt') as pafFile:
        for aln in parse_paf(pafFile):
            if aln.get_tag('tp') == 'P' and aln.mapq >= opt.minMapq:
                qry = (aln.query_start,aln.query_end,aln.query_length) if aln.strand == '+' else (aln.query_length-aln.query_end,aln.query_length-aln.query_start,aln.query_length)
                ref = (aln.reference_start,aln.reference_end,aln.reference_length)
                maptype = classify_mapping(qry,ref)
                if maptype == MappingType.DOVETAIL_PREFIX:
                    read2contigs[aln.query].append((aln.reference,MappingType.DOVETAIL_SUFFIX,aln))
                if maptype == MappingType.DOVETAIL_SUFFIX:
                    read2contigs[aln.query].append((aln.reference,MappingType.DOVETAIL_PREFIX,aln))
                # TODO: handle case where a read have alignments different to DOVETAIL
                #if maptype == MappingType.SECOND_CONTAINED: # read spanning a whole contig

    # add edges between contigs linked by dovetail read mappings
    for read in read2contigs:
        if len(read2contigs[read])==2: # read linking exactly two contigs with dovetail alignment
            fst,snd=read2contigs[read]
            a,atype,aaln=fst
            b,btype,baln=snd
            if a==b and atype==btype:
                continue
            # avoid linking certain type of contigs
            afrag=getRefPos(a)
            aref=afrag[0]
            bfrag=getRefPos(b)
            bref=bfrag[0]
            is_reflink=(afrag,bfrag) in psAdjSet|fragAdjSet
            is_singlink=(a in asmContigs and len(refEnds[bref])==1) or (b in asmContigs and len(refEnds[aref])==1)
            is_endlink=(afrag in refEnds[aref] and bfrag in refEnds[bref])
            if is_reflink or is_singlink or is_endlink:
                if not asmGraph.has_edge(a,b):
                    asmGraph.add_edge(a,b,nreads=0,etype={a:atype,b:btype},label="0",overlaps=[],color="black") # TODO: assuming no self-loops or loops of size 2, fix that
                    asmGraph[a][b]['ctglink']=set()
                    if a in asmContigs and b in asmContigs and (getRefPos(a),getRefPos(b)) in psAdjSet:
                        asmGraph[a][b]['color']="red"
                    if atype == MappingType.DOVETAIL_SUFFIX:
                        asmGraph.nodes[a]['outedges'].add(b)
                    elif atype == MappingType.DOVETAIL_PREFIX:
                        asmGraph.nodes[a]['inedges'].add(b)
                    if btype == MappingType.DOVETAIL_SUFFIX:
                        asmGraph.nodes[b]['outedges'].add(a)
                    elif btype == MappingType.DOVETAIL_PREFIX:
                        asmGraph.nodes[b]['inedges'].add(a)
                arange=(aaln.query_start,aaln.query_end,aaln.query_length)
                brange=(baln.query_start,baln.query_end,baln.query_length)
                asmGraph[a][b]['overlaps'].append(overlap_length(arange,brange))
                asmGraph[a][b]['nreads']+=1
                asmGraph[a][b]['label']=f'{asmGraph[a][b]["nreads"]}'

    nx.nx_agraph.to_agraph(asmGraph).write(f'{opt.outPrefix}.raw.dot')
    
    asmGraph=remove_simple_transitive_edges(asmGraph,asmContigs,npContigs)
    nx.nx_agraph.to_agraph(asmGraph).write(f'{opt.outPrefix}.simplified.dot')
    
    asmGraph=remove_weak_edges(asmGraph,opt.minReads,opt.minReadFrac)
    nx.nx_agraph.to_agraph(asmGraph).write(f'{opt.outPrefix}.final.dot')
    
    eprint(f'connected-components: {len(list(nx.connected_components(asmGraph)))}')
    eprint(f'N50-after: {estimate_n50(asmGraph)}')

    scaffolds=build_scaffolds(asmGraph,ctgSeq)
    with open(f'{opt.outPrefix}.fa','w') as of:
        for i,scaff in enumerate(scaffolds,1):
            of.write(f'>scaffold_{i}\n')
            of.write(f'{insert_newlines(scaff)}\n')

    return 0


# Check if the program is not being imported
if __name__ == "__main__":
    sys.exit(main())

