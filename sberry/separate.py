import sys,os,shutil,pathlib
import subprocess,multiprocessing

from Bio import SeqIO,bgzf

from sberry.utils import *


# TODO:
# output logging information from longshot processes
# check if files exist

class LongshotReadSeparator(object):

    def __init__(self,fastafile,bamfile,outdir,qual=50,nproc=1):
        self.LONGSHOT_BIN='longshot'
        self.TABIX_BIN='tabix'
        # input parameters
        self.fasta=fastafile
        self.bam=bamfile
        self.qual=qual
        self.nproc=nproc
        # output paths
        self.separation_dir=os.path.join(outdir,'10-separation')
        self.vcfs_dir=os.path.join(self.separation_dir,'vcfs')
        self.tagged_dir=os.path.join(self.separation_dir,'tagged')
        self.phased_vcf=os.path.join(self.separation_dir,'variants.phased.vcf.gz')
        # create output directories
        os.makedirs(self.vcfs_dir,exist_ok=True)
        os.makedirs(self.tagged_dir,exist_ok=True)

    def _run_longshot(self,contig):
        longshot_cmd = [ self.LONGSHOT_BIN, '--region', contig,
                '--bam', self.bam, '--ref', self.fasta,
                '--out', os.path.join(self.vcfs_dir,f'{contig}.vcf'),
                '--out_bam', os.path.join(self.tagged_dir,f'{contig}.bam'),
                '--force_overwrite' ]
        longshot = subprocess.run(longshot_cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        return longshot.returncode

    def run(self):
        # run longshot in parallel for each contig
        contigs = [ ctg.id for ctg in SeqIO.parse(self.fasta,'fasta') ]
        pool = multiprocessing.Pool(processes=self.nproc,maxtasksperchild=1)
        results = [ pool.apply_async(self._run_longshot, args=(ctg,)) for ctg in contigs ]
        pool.close()
        if any( r.get() != 0 for r in results ):
            print_error(f'error running {self.LONGSHOT_BIN}')
            return False
        # merge all vcf files in a single one
        first_header=True
        with bgzf.BgzfWriter(self.phased_vcf, "wb") as out:
            for path in pathlib.Path(self.vcfs_dir).glob('*.vcf'):
                with open(path,'r') as vcf:
                    for line in vcf:
                        if line.startswith('#') and not first_header:
                            continue
                        out.write(line)
                first_header=False
        # index phased bgzipped vcf file
        tabix = subprocess.run([self.TABIX_BIN,'-fp','vcf',self.phased_vcf])
        if tabix.returncode != 0:
            print_error(f'error indexing {self.phased_vcf}')
            return False

        return True

