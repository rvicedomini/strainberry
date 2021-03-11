import sys,os,shutil,gzip
import subprocess,multiprocessing
from collections import defaultdict

import pysam

from sberry.utils import *
from sberry.phaseset import Phaseset,PhasesetCollection


# TODO:
# output logging information from longshot processes
# check if files exist

class LongshotReadSeparator(object):

    def __init__(self,fastafile,bamfile,outdir,rmtemp=False,snv_dens=0.1,ps_length=3000,min_nreads=10,qual=50,nproc=1):
        self.LONGSHOT_BIN='longshot'
        self.SAMTOOLS_BIN='samtools'
        self.TABIX_BIN='tabix'
        # output results
        self.contigs={}
        self.psc=PhasesetCollection()
        self.is_phased=False
        self.hap_set=set()
        # input parameters
        self.fasta=fastafile
        self.bam=bamfile
        self.rmtemp=rmtemp
        self.snv_dens=snv_dens
        self.ps_length=ps_length
        self.min_nreads=min_nreads
        self.qual=qual
        self.nproc=nproc
        # output paths
        self.separation_dir=os.path.join(outdir,'10-separation')
        self.logs_dir=os.path.join(self.separation_dir,'logs')
        self.vcfs_dir=os.path.join(self.separation_dir,'vcfs')
        self.tagged_dir=os.path.join(self.separation_dir,'tagged')
        self.reads_dir=os.path.join(self.separation_dir,'fasta')
        # create output directories
        os.makedirs(self.logs_dir,exist_ok=True)
        os.makedirs(self.vcfs_dir,exist_ok=True)
        os.makedirs(self.tagged_dir,exist_ok=True)
        os.makedirs(self.reads_dir,exist_ok=True)

    def _run_longshot(self,contig):
        ctg_vcffile = os.path.join(self.vcfs_dir,f'{contig}.vcf')
        ctg_bamfile = os.path.join(self.tagged_dir,f'{contig}.bam')
        ctg_logfile = os.path.join(self.logs_dir,f'{contig}.log')
        # run longshot
        with open(ctg_logfile,'w') as log:
            longshot_cmd = [ self.LONGSHOT_BIN, '--region', contig, '--bam', self.bam, '--ref', self.fasta, '--out', ctg_vcffile, '--out_bam', ctg_bamfile, '--force_overwrite' ]
            longshot = subprocess.run(longshot_cmd, stdout=log, stderr=log)
            # index output bam file
            if longshot.returncode == 0 and os.path.isfile(ctg_bamfile):
                if subprocess.run(['samtools','index',ctg_bamfile]).returncode != 0:
                    print_error(f'cannot index bam file {ctg_bamfile}')
            return longshot.returncode

    def _update_extracted_reads(self,hap_nreads):
        self.hap_set.update( hap_id for hap_id,nreads in hap_nreads.items() if nreads >= self.min_nreads )

    def _extract_tagged_reads(self,ctg_bamfile,phaseset):
        hap_files={}
        hap_nreads=defaultdict(int)
        with pysam.AlignmentFile(ctg_bamfile,'rb') as bam:
            for read in bam.fetch( phaseset.reference, start=phaseset.start(), stop=phaseset.end() ):
                # TODO: is the following really necessary?
                if read.is_unmapped or read.is_secondary or read.is_supplementary:
                    continue
                if read.has_tag('PS') and f'{read.get_tag("PS")}' == phaseset.psid:
                    hp_tag=read.get_tag("HP")
                    hap_id=(phaseset.reference,phaseset.psid,f'{hp_tag}')
                    if not hap_id in hap_files:
                        hap_filename=f'{phaseset.reference}_{phaseset.psid}_h{hp_tag}.fa.gz'
                        hap_files[hap_id]=gzip.open(os.path.join(self.reads_dir,hap_filename),'wt')
                    hap_files[hap_id].write(f'>{read.query_name}\n{read.query_sequence}\n')
                    hap_nreads[hap_id]+=1
        for fh in hap_files.values():
            fh.close()
        return hap_nreads

    def phase_and_tag(self):
        # load contig IDs and lengths
        with open(f'{self.fasta}.fai','r') as faidx:
            for line in faidx:
                ctg_id,ctg_len,_=line.split('\t',2)
                self.contigs[ctg_id]=int(ctg_len)
        # run longshot in parallel for each contig
        with multiprocessing.Pool(processes=self.nproc,maxtasksperchild=1) as pool:
            results = [ pool.apply_async(self._run_longshot, args=(ctg,)) for ctg in self.contigs ]
            pool.close()
            if any( r.get() != 0 for r in results ):
                print_error(f'error running {self.LONGSHOT_BIN}')
                return False
        self.is_phased=True
        # load phasesets from phased VCFs
        for ctg in self.contigs:
            ctg_vcffile = os.path.join(self.vcfs_dir,f'{ctg}.vcf')
            ctg_psc = PhasesetCollection()
            ctg_psc.load_from_vcf(ctg_vcffile, min_density=self.snv_dens, min_length=self.ps_length)
            self.psc.update(ctg_psc)
        return True

    def separate_reads(self):
        if not self.is_phased:
            print_error(f'separate_reads called before phase_and_tag method')
            return False
        # extract tagged reads in parallel for each phaseset
        with multiprocessing.Pool(processes=self.nproc,maxtasksperchild=1) as pool:
            results = []
            for ctg_id in self.contigs:
                ctg_bamfile = os.path.join(self.tagged_dir,f'{ctg_id}.bam')
                for phaseset in self.psc.phasesets(ctg_id):
                    results.append( pool.apply_async(self._extract_tagged_reads, args=(ctg_bamfile,phaseset,), callback=self._update_extracted_reads, error_callback=eprint) )
            pool.close()
            pool.join()
        # possibly remove phased VCFs and tagged BAMs after read separation
        if self.rmtemp:
            print_status(f'temporary phased-vcf and tagged-bam directories will be deleted')
            shutil.rmtree(self.vcfs_dir,ignore_errors=True)
            shutil.rmtree(self.tagged_dir,ignore_errors=True)
        return True

