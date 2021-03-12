import sys,os,gzip
import multiprocessing,subprocess
from collections import defaultdict

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

from sberry.utils import *

class WtdbgAssembler(object):
    
    def __init__(self,reffile,psc,hap_set,reads_dir,work_dir,is_ont=False,rmtemp=False,nproc=1):
        self.mindepth=5 # min edge depth for wtdbg2
        self.min_unphased_len=500 # minimum length of unphased regions
        # input parameters
        self.reffile=reffile
        self.psc=psc
        self.hap_set=hap_set
        self.reads_dir=reads_dir
        self.work_dir=work_dir
        self.is_ont=is_ont
        self.rmtemp=rmtemp
        self.nproc=nproc
        # output
        self.asm_list=[]
        self.asm_regions=defaultdict(list)
        self.wtdbg2_dir=os.path.join(self.work_dir,'wtdbg2')
        os.makedirs(self.wtdbg2_dir,exist_ok=True)
        self.phaseset_dir=os.path.join(self.work_dir,'phasesets')
        os.makedirs(self.phaseset_dir,exist_ok=True)
        self.trimmed_dir=os.path.join(self.work_dir,'trimmed')
        os.makedirs(self.trimmed_dir,exist_ok=True)
        self.unphased_file=os.path.join(self.work_dir,'unphased.contigs.fa')
        self.contig_file=os.path.join(self.work_dir,'assembly.contigs.fa')
        #self.contig_desc_file=os.path.join(self.work_dir,'assembly.contigs.desc.tsv')

    def _is_gz_empty(self,filename):
        try:
            with gzip.open(filename,'rb') as f:
                return len(f.read(1)) == 0
        except:
            #print_warning(f'{filename} is a bad gzip file')
            return True

    def _run_wtdbg2(self,hapid,fagz_file):
        contig,psid,hpid = hapid
        asm_id = f'{contig}_{psid}_h{hpid}'
        hapasm_dir = os.path.join(self.wtdbg2_dir,asm_id)
        os.makedirs(hapasm_dir,exist_ok=True)
        out_prefix = os.path.join(hapasm_dir,asm_id)
        with open(f'{out_prefix}.log','w') as logfile:
            # run wtdbg2
            gzpipe = subprocess.Popen(['gzip','-dc',fagz_file], stdout=subprocess.PIPE)
            wtdbg2 = subprocess.Popen(['wtdbg2', '-x','preset2', '-t','2', '-e',str(self.mindepth), '-l','1000', '-L','3000', '-S','1', '-R', '-o',out_prefix, '-i','-'],
                    stdin=gzpipe.stdout, stdout=logfile, stderr=subprocess.STDOUT)
            wtdbg2.communicate()
            if wtdbg2.returncode != 0 or not os.path.isfile(f'{out_prefix}.ctg.lay.gz') or self._is_gz_empty(f'{out_prefix}.ctg.lay.gz'):
                return (False,None)
            # run wtpoa-cns
            wtpoa_cns = subprocess.run(['wtpoa-cns', '-t','2', '-i',f'{out_prefix}.ctg.lay.gz', '-fo',f'{out_prefix}.ctg.fa'], 
                    stdout=logfile, stderr=subprocess.STDOUT)
            if wtpoa_cns.returncode != 0 or not os.path.isfile(f'{out_prefix}.ctg.fa') or os.stat(f'{out_prefix}.ctg.fa').st_size == 0:
                return (False,None)
            # run polishing
            mm2_preset = 'map-ont' if self.is_ont else 'map-pb'
            minimap2 = subprocess.Popen(['minimap2', '-t2', '-ax',mm2_preset, '-r2k', f'{out_prefix}.ctg.fa', fagz_file], stdout=subprocess.PIPE, stderr=logfile )
            samtools_sort = subprocess.Popen(['samtools','sort','-@2'], stdin=minimap2.stdout, stdout=subprocess.PIPE, stderr=logfile )
            samtools_view = subprocess.Popen(['samtools','view','-F0x900'], stdin=samtools_sort.stdout, stdout=subprocess.PIPE, stderr=logfile )
            wtpoa_cns = subprocess.Popen(['wtpoa-cns', '-t','2', '-d',f'{out_prefix}.ctg.fa', '-i','-', '-fo',f'{out_prefix}.cns.fa'], 
                    stdin=samtools_view.stdout, stdout=logfile, stderr=subprocess.STDOUT)
            wtpoa_cns.communicate()
            if wtpoa_cns.returncode != 0 or not os.path.isfile(f'{out_prefix}.cns.fa') or os.stat(f'{out_prefix}.cns.fa').st_size == 0:
                return (False,None)
        # TODO: get rid of the following, I think it might not be really necessary...
        # rename of contigs
        with open( os.path.join(self.wtdbg2_dir,f'{asm_id}.fa'), 'w' ) as ofh:
            for i,ctg in enumerate(SeqIO.parse(f'{out_prefix}.cns.fa','fasta'),start=1):
                ctg_rec = SeqRecord(ctg.seq, id=f'contig_{i}', description=f'phased=true reference={contig} ps={psid} hp={hpid}')
                SeqIO.write(ctg_rec,ofh,'fasta')
        return (True,hapid)

    def _update_assembled_haps(self,res):
        is_valid,hapid=res
        if is_valid:
            self.asm_list.append(hapid)

    def assemble_reads(self):
        with multiprocessing.Pool(processes=(1+self.nproc//2),maxtasksperchild=1) as pool:
            for hapid in self.hap_set:
                contig,psid,hpid=hapid
                hap_fagz = os.path.join(self.reads_dir,f'{contig}_{psid}_h{hpid}.fa.gz')
                pool.apply_async(self._run_wtdbg2, args=(hapid,hap_fagz,), callback=self._update_assembled_haps, error_callback=eprint)
            pool.close()
            pool.join()

    def _write_phaseset(self,phaseset,farecord):
        ps_file=os.path.join(self.phaseset_dir,f'{phaseset.reference}_{phaseset.psid}.fa')
        with open(ps_file,'w') as fh:
            ps_rec = SeqRecord(farecord.seq[phaseset.start():phaseset.end()], id=f'{phaseset.reference}_{phaseset.psid}')
            SeqIO.write(ps_rec,fh,'fasta')

    def split_asm_regions(self):
        fadict=SeqIO.to_dict(SeqIO.parse(self.reffile,'fasta'))
        asm_phasesets = set( (contig,psid) for contig,psid,_ in self.asm_list )
        # collect assemble phasesets and write them in output for trimming
        for ctg,psid in asm_phasesets:
            phaseset=self.psc.phaseset(ctg,psid)
            self.asm_regions[ctg].append(phaseset)
            self._write_phaseset(phaseset,fadict[ctg])
        # write unphased regions
        with open( self.unphased_file,'w') as fh:
            for ctg in fadict:
                ps_list=sorted(self.asm_regions[ctg],key=lambda ps:ps.start())
                i=0
                start=0
                while i < len(ps_list):
                    phaseset=ps_list[i]
                    if phaseset.start()-start >= self.min_unphased_len:
                        end=phaseset.start()
                        out_rec = SeqRecord(fadict[ctg].seq[start:end], id=f'{ctg}_{start+1}', description=f'phased=false reference={ctg} start={start+1} end={end}')
                        SeqIO.write(out_rec,fh,'fasta')
                    start=phaseset.end()
                    i+=1
                    while i < len(ps_list) and ps_list[i].start() <= start:
                        start=max(start,ps_list[i].end())
                        i+=1
                ctg_len=len(fadict[ctg])
                if ctg_len-start >= self.min_unphased_len:
                    out_rec = SeqRecord(fadict[ctg].seq[start:ctg_len], id=f'{ctg}_{start+1}', description=f'phased=false reference={ctg} start={start+1} end={ctg_len}')
                    SeqIO.write(out_rec,fh,'fasta')
        return True

    # TODO: in parallel
    def trim_assemblies(self):
        for contig,psid,hpid in self.asm_list:
            asm_id = f'{contig}_{psid}_h{hpid}'
            asm_file = os.path.join(self.wtdbg2_dir,f'{asm_id}.fa')
            ps_file = os.path.join(self.phaseset_dir,f'{contig}_{psid}.fa')
            nucmer_prefix = os.path.join(self.trimmed_dir,f'{contig}_{psid}_h{hpid}')
            with open(f'{nucmer_prefix}.1delta', 'w') as delta, open(f'{nucmer_prefix}.1coords', 'w') as coords, open(f'{nucmer_prefix}.err', 'w') as errlog:
                nucmer = subprocess.run(['nucmer', '-p',nucmer_prefix, ps_file, asm_file], stderr=errlog)
                delta_filter = subprocess.run(['delta-filter', '-1', f'{nucmer_prefix}.delta'], stdout=delta, stderr=errlog)
                show_coords = subprocess.run(['show-coords', '-rlcHT', f'{nucmer_prefix}.1delta'], stdout=coords, stderr=errlog)
                if nucmer.returncode != 0 or delta_filter.returncode != 0 or show_coords.returncode != 0:
                    print_error(f'failed trimming assembly {asm_file} with phaseset file {ps_file}')
                    print_error(f'see error log at: {nucmer_prefix}.err')
                    return False
            coord_dict={}
            with open(f'{nucmer_prefix}.1coords','r') as coords:
                for line in coords:
                    aln = line.rstrip().split('\t')
                    asm_ctg,asm_start,asm_end = aln[12],int(aln[2]),int(aln[3])
                    if asm_start > asm_end:
                        asm_start,asm_end = asm_end,asm_start
                    if not asm_ctg in coord_dict:
                        coord_dict[asm_ctg]=(asm_start,asm_end)
                    else:
                        ctg_start,ctg_end=coord_dict[asm_ctg]
                        coord_dict[asm_ctg]=(min(asm_start,ctg_start),max(asm_end,ctg_end))
            asm_dict=SeqIO.to_dict(SeqIO.parse(asm_file,'fasta'))
            trimmed_file=os.path.join(self.trimmed_dir,f'{asm_id}.fa')
            with open(trimmed_file,'w') as fh:
                for asm_ctg,asm_region in coord_dict.items():
                    asm_start,asm_end=asm_region
                    asm_rec = asm_dict[asm_ctg]
                    out_rec = SeqRecord(asm_rec.seq[asm_start-1:asm_end], id=asm_ctg, description=asm_rec.description)
                    SeqIO.write(out_rec,fh,'fasta')

    def write_result(self):
        with open(self.contig_file,'w') as ofh:
            count=1
            # write separated contigs
            for contig,psid,hpid in self.asm_list:
                asm_id=f'{contig}_{psid}_h{hpid}'
                asm_file=os.path.join(self.trimmed_dir,f'{asm_id}.fa')
                for rec in SeqIO.parse(asm_file,'fasta'):
                    rec.id=f'ctg_{count}'
                    rec.description=rec.description.split(maxsplit=1)[1]
                    SeqIO.write(rec,ofh,'fasta')
                    count+=1
            # write unphased contigs
            for rec in SeqIO.parse(self.unphased_file,'fasta'):
                rec.id=f'ctg_{count}'
                rec.description=rec.description.split(maxsplit=1)[1]
                SeqIO.write(rec,ofh,'fasta')
                count+=1

    def run(self):
        self.assemble_reads()
        self.split_asm_regions()
        self.trim_assemblies()
        self.write_result()

