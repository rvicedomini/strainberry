#!/usr/bin/env python3

import sys,os,gzip,subprocess
import networkx as nx

from collections import defaultdict
from statistics import median
from enum import Enum
from Bio import SeqIO

from sberry.utils import *
from sberry.alignment import *


def fix_neighbors(asmGraph):
    for a in asmGraph.nodes:
        toremove = set( b for b in asmGraph.nodes[a]['inedges'] if not asmGraph.has_edge(a,b) )
        asmGraph.nodes[a]['inedges'] -= toremove
        toremove = set( b for b in asmGraph.nodes[a]['outedges'] if not asmGraph.has_edge(a,b) )
        asmGraph.nodes[a]['outedges'] -= toremove
    return asmGraph

def remove_simple_transitive_edges(asmGraph):
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
                    if b!=c and asmGraph.has_edge(a,c) and asmGraph[a][c]['etype'][c] == bc_ctype and asmGraph[a][c]['etype'][a] == ab_atype and len(asmGraph.nodes[b]['outedges'])==1 and len(asmGraph.nodes[b]['inedges'])==1:
                        removeEdgeSet.add((a,c))
                        asmGraph[a][b]['nreads']+=asmGraph[a][c]['nreads']
                        asmGraph[b][c]['nreads']+=asmGraph[a][c]['nreads']
            elif a!=b and asmGraph[a][b]['etype'][b] == MappingType.DOVETAIL_SUFFIX:
                for c in asmGraph.nodes[b]['inedges']:
                    bc_ctype=asmGraph[b][c]['etype'][c]
                    if b!=c and asmGraph.has_edge(a,c) and asmGraph[a][c]['etype'][c] == bc_ctype and asmGraph[a][c]['etype'][a] == ab_atype and len(asmGraph.nodes[b]['outedges'])==1 and len(asmGraph.nodes[b]['inedges'])==1:
                        removeEdgeSet.add((a,c))
                        asmGraph[a][b]['nreads']+=asmGraph[a][c]['nreads']
                        asmGraph[b][c]['nreads']+=asmGraph[a][c]['nreads']
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
    out_desc=[]
    for a in asmGraph.nodes:
        if a not in visited:
            visited.add(a)
            aseq=ctgSeq[a]
            adesc=[(a,0)]
            outedges=asmGraph.nodes[a]['outedges']
            if len(outedges)==1 and list(outedges)[0] not in visited:
                for b in asmGraph.nodes[a]['outedges']:
                    edata=asmGraph[a][b]
                    atype=edata['etype'][a]
                    btype=edata['etype'][b]
                    bseq,bdesc=build_scaffold_from(b,btype,asmGraph,ctgSeq,visited)
                    ab_gap=int(median(edata["overlaps"]))
                    ab_gap=-ab_gap if ab_gap<0 else 0
                    aseq+=(ab_gap*'N') + (bseq if atype!=btype else reverse_complement(bseq))
                    adesc.extend( bdesc if atype!=btype else ((b,1-rc) for b,rc in reversed(bdesc)) )
            inedges=asmGraph.nodes[a]['inedges']
            if len(inedges)==1 and list(inedges)[0] not in visited:
                for b in asmGraph.nodes[a]['inedges']:
                    edata=asmGraph[a][b]
                    atype=edata['etype'][a]
                    btype=edata['etype'][b]
                    bseq,bdesc=build_scaffold_from(b,btype,asmGraph,ctgSeq,visited)
                    ab_gap=int(median(edata["overlaps"]))
                    ab_gap=-ab_gap if ab_gap<0 else 0
                    aseq=(bseq if atype!=btype else reverse_complement(bseq)) + (ab_gap*'N') + aseq
                    adesc=(bdesc if atype!=btype else [(b,1-rc) for b,rc in reversed(bdesc)]) + adesc
            out_scaffolds.append(aseq)
            out_desc.append(adesc)
    return out_scaffolds,out_desc

# TODO: AVOID RECURSION => MAKE IT ITERATIVE; also there are maybe too many reverse complement operations
def build_scaffold_from(a,atype,asmGraph,ctgSeq,visited):
    # assuming this function is called on an unvisited node
    visited.add(a)
    aseq=ctgSeq[a]
    adesc=[(a,0)]
    if atype == MappingType.DOVETAIL_PREFIX:
        outedges=asmGraph.nodes[a]['outedges']
        if len(outedges) == 1 and list(outedges)[0] not in visited:
            for b in asmGraph.nodes[a]['outedges']:
                edata=asmGraph[a][b]
                atype=edata['etype'][a]
                btype=edata['etype'][b]
                bseq,bdesc=build_scaffold_from(b,btype,asmGraph,ctgSeq,visited)
                ab_gap=int(median(edata["overlaps"]))
                ab_gap=-ab_gap if ab_gap<0 else 0
                aseq+=(ab_gap*'N') + (bseq if atype!=btype else reverse_complement(bseq))
                adesc.extend( bdesc if atype!=btype else ((b,1-rc) for b,rc in reversed(bdesc)) )
    elif atype == MappingType.DOVETAIL_SUFFIX:
        inedges=asmGraph.nodes[a]['inedges']
        if len(inedges)==1 and list(inedges)[0] not in visited:
            for b in asmGraph.nodes[a]['inedges']:
                edata=asmGraph[a][b]
                atype=edata['etype'][a]
                btype=edata['etype'][b]
                bseq,bdesc=build_scaffold_from(b,btype,asmGraph,ctgSeq,visited)
                ab_gap=int(median(edata["overlaps"]))
                ab_gap=-ab_gap if ab_gap<0 else 0
                aseq=(bseq if atype!=btype else reverse_complement(bseq)) + (ab_gap*'N') + aseq
                bdesc=(bdesc if atype!=btype else [(b,1-rc) for b,rc in reversed(bdesc)]) + adesc
    return aseq,adesc


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



class Scaffolder(object):

    def __init__(self,fastafile,bamfile,outdir,is_ont=False,min_mapq=40,min_nreads=10,min_read_frac=0.75,nproc=1):
        self.SAMTOOLS_BIN='samtools'
        self.MINIMAP2_BIN='minimap2'
        # input parameters
        self.fasta=fastafile
        self.bam=bamfile
        self.is_ont=is_ont
        self.min_mapq=min_mapq # minimum MAPQ value to consider a read alignment during scaffolding
        self.min_nreads=min_nreads # minimum number of reads matching the prefix/suffix of a contig
        self.min_read_frac=min_read_frac # minimum read fraction to define good edges
        self.nproc=nproc
        # output paths
        self.scaffold_file=os.path.join(outdir,'assembly.scaffolds.fa')
        self.scaffold_info_file=os.path.join(outdir,'assembly.scaffolds.info.tsv')
        self.scaffold_dir=os.path.join(outdir,'30-scaffolds')
        self.pafgz=os.path.join(self.scaffold_dir,'alignment.paf.gz')
        self.prefix=os.path.join(self.scaffold_dir,'scaffolds')
        # create output directories
        os.makedirs(self.scaffold_dir,exist_ok=True)

    
    # Map input reads to input assembly
    def _align_reads(self):
        #print_status('mapping reads to strain-separated contigs')
        with open(self.pafgz,'wb') as pafgz:
            samtools_fastq = subprocess.Popen([self.SAMTOOLS_BIN,'fastq',self.bam], stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
            minimap2_cmd = [ self.MINIMAP2_BIN, '-cx', 'map-ont' if self.is_ont else 'map-pb', '-t', f'{self.nproc}', self.fasta, '-' ]
            minimap2 = subprocess.Popen(minimap2_cmd, stdin=samtools_fastq.stdout, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
            gzpipe = subprocess.Popen(['gzip','-c'], stdin=minimap2.stdout, stdout=pafgz, stderr=subprocess.DEVNULL)
            gzpipe.communicate()
            if gzpipe.returncode != 0:
                print_error(f'read mapping to strain-separated contigs failed!')
                return False
        return True
    
    # Parse paf alignment
    def _parse_paf(self):
        # parse alignments to find dovetail mappings
        read2contigs=defaultdict(list)
        with gzip.open(self.pafgz,'rt') as pafFile:
            for aln in parse_paf(pafFile):
                if aln.get_tag('tp') == 'P' and aln.mapq >= self.min_mapq:
                    qry = (aln.query_start,aln.query_end,aln.query_length) if aln.strand == '+' else (aln.query_length-aln.query_end,aln.query_length-aln.query_start,aln.query_length)
                    ref = (aln.reference_start,aln.reference_end,aln.reference_length)
                    maptype = classify_mapping(qry,ref)
                    if maptype == MappingType.DOVETAIL_PREFIX:
                        read2contigs[aln.query].append((aln.reference,MappingType.DOVETAIL_SUFFIX,aln))
                    if maptype == MappingType.DOVETAIL_SUFFIX:
                        read2contigs[aln.query].append((aln.reference,MappingType.DOVETAIL_PREFIX,aln))
                    # TODO: handle case where a read have alignments different to DOVETAIL
        return read2contigs

    def _getRefPos(self,ctg_id):
        ctg_info = self.ctgInfo[ctg_id]
        if ctg_info['phased'] == 'true':
            return (ctg_info['reference'],int(ctg_info['ps']))
        return (ctg_info['reference'],int(ctg_info['start']))

    # Run the scaffolding procedure
    def run(self):

        # map reads to input (contig) assembly
        if not self._align_reads():
            return False

        # retrieve identifiers and put not-phased/assembled contigs in two different sets
        self.ctgInfo={}
        self.ctgSeq={}
        self.npContigs={}  #TODO: remove this and use only self.ctgInfo
        self.asmContigs={} #TODO: remove this and use only self.ctgInfo
        for rec in SeqIO.parse(self.fasta,'fasta'):
            self.ctgInfo[rec.id] = { key:val for key,val in (x.split('=') for x in rec.description.split() if '=' in x) }
            self.ctgSeq[rec.id] = str(rec.seq)
            if self.ctgInfo[rec.id]['phased'] == 'true':
                self.asmContigs[rec.id] = len(rec.seq)
            else:
                self.npContigs[rec.id] = len(rec.seq)

        #ctgSeq={ ctg.id:str(ctg.seq) for ctg in SeqIO.parse(self.fasta,'fasta') }
        #ctgDict={ ctg_id:len(ctg_seq) for ctg_id,ctg_seq in ctgSeq.items() }
        #npContigs={ ctg:ctglen for ctg,ctglen in ctgDict.items() if not ctg.startswith('sberry|') }
        #asmContigs={ ctg:ctglen for ctg,ctglen in ctgDict.items() if ctg.startswith('sberry|') }

        psDict=defaultdict(list)
        psAdjSet=set()
        for ctg in self.asmContigs.keys():
            ps=self._getRefPos(ctg)
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

        fragList=sorted(map(self._getRefPos,self.ctgInfo.keys()))
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

        # create graph and contig->node_id mapping
        asmGraph=nx.Graph(style="filled")
        asmGraph.add_nodes_from(list(self.npContigs.keys()),ctgtype='notphased')
        asmGraph.add_nodes_from(list(self.asmContigs.keys()),ctgtype='assembled')
        nx.set_node_attributes(asmGraph,self.npContigs,name='length')
        nx.set_node_attributes(asmGraph,self.asmContigs,name='length')
        for n in asmGraph.nodes:
            asmGraph.nodes[n]['outedges']=set()
            asmGraph.nodes[n]['inedges']=set()

        #print_status(f'connected-components: {len(asmGraph.nodes)}')
        #print_status(f'N50-before: {estimate_n50(asmGraph)}')

        for n in asmGraph.nodes:
            if n in self.asmContigs:
                asmGraph.nodes[n]['style']='filled'
                asmGraph.nodes[n]['fillcolor']=f'/set312/1'
        
        read2contigs = self._parse_paf()

        # add edges between contigs linked by dovetail read mappings
        for read in read2contigs:
            if len(read2contigs[read])==2: # read linking exactly two contigs with dovetail alignment
                fst,snd=read2contigs[read]
                a,atype,aaln=fst
                b,btype,baln=snd
                if a==b and atype==btype:
                    continue
                # avoid linking certain type of contigs
                afrag=self._getRefPos(a)
                aref=afrag[0]
                bfrag=self._getRefPos(b)
                bref=bfrag[0]
                is_reflink=(afrag,bfrag) in psAdjSet|fragAdjSet
                is_singlink=(a in self.asmContigs and len(refEnds[bref])==1) or (b in self.asmContigs and len(refEnds[aref])==1)
                is_endlink=(afrag in refEnds[aref] and bfrag in refEnds[bref])
                if is_reflink or is_singlink or is_endlink:
                    if not asmGraph.has_edge(a,b):
                        asmGraph.add_edge(a,b,nreads=0,etype={a:atype,b:btype},label="0",overlaps=[],color="black") # FIXME: assuming no self-loops or loops of size 2
                        asmGraph[a][b]['ctglink']=set()
                        if a in self.asmContigs and b in self.asmContigs and (self._getRefPos(a),self._getRefPos(b)) in psAdjSet:
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

        nx.nx_agraph.to_agraph(asmGraph).write(f'{self.prefix}.raw.dot')
        
        asmGraph=remove_simple_transitive_edges(asmGraph)
        nx.nx_agraph.to_agraph(asmGraph).write(f'{self.prefix}.simplified.dot')
        
        asmGraph=remove_weak_edges(asmGraph,self.min_nreads,self.min_read_frac)
        nx.nx_agraph.to_agraph(asmGraph).write(f'{self.prefix}.final.dot')
        
        #print_status(f'connected-components: {len(list(nx.connected_components(asmGraph)))}')
        #print_status(f'N50-after: {estimate_n50(asmGraph)}')

        scaffolds,scf_info=build_scaffolds(asmGraph,self.ctgSeq)
        with open(self.scaffold_file,'w') as scf_fh, open(self.scaffold_info_file,'w') as info_fh:
            for i,scf in enumerate(scaffolds):
                scf_id = f'scaffold_{i+1}'
                scf_fh.write(f'>{scf_id}\n')
                scf_fh.write(f'{insert_newlines(scf)}\n')
                for ctg_id,rev in scf_info[i]:
                    record = [ scf_id, ctg_id, '+' if rev==0 else '-' ]
                    record.extend(f'{key}={val}' for key,val in self.ctgInfo[ctg_id].items())
                    info_fh.write('\t'.join(record) + '\n')

        return True

