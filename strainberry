#!/usr/bin/env python3

import sys,os,argparse,subprocess

from sberry.utils import *
from sberry.separate import LongshotReadSeparator
from sberry.phaseset import PhasesetCollection
from sberry.hratio import AverageHammingRatio
from sberry.assemble import WtdbgAssembler
from sberry.scaffolder import Scaffolder
from sberry.__version__ import __version__


sberry_root = os.path.dirname(os.path.abspath(os.path.realpath(__file__)))
sberry_scripts = os.path.join(sberry_root,'sberry')


def _version():
    sberry_version = f'{__version__}'
    git = subprocess.run(['git', '-C', sberry_root, 'describe', '--always'], capture_output=True, stderr=None, text=True)
    if git.returncode == 0:
        git_hash = git.stdout.strip().rsplit('-',1)[-1]
        sberry_version += f'-{git_hash}'
    return sberry_version


def main(argv=None):
    parser = argparse.ArgumentParser(prog='strainberry', description='Automated strain separation of low-complexity metagenomes', add_help=False)
    # MANDATORY ARGUMENTS
    required_args = parser.add_argument_group('Required arguments')
    required_args.add_argument('-r','--reference', dest='FASTA', metavar='PATH', required=True, help='Strain-oblivious assembly in FASTA format')
    required_args.add_argument('-b','--bam', dest='BAM', metavar='PATH', required=True, help='Read alignment in BAM format')
    required_args.add_argument('-o','--out-dir', dest='OUTDIR', metavar='PATH', required=True, help='Output directory of Strainberry assemblies')
    # OPTIONAL ARGUMENTS
    optional_args = parser.add_argument_group('Optional arguments')
    optional_args.add_argument('--nanopore', dest='ONT', action='store_true', help='Input consists of Oxford Nanopore raw reads')
    optional_args.add_argument('-n','--max-strains', dest='STRAINS', metavar='int', type=int, default=5, help='Attempt strain-separation at most for the provided strain multiplicity [%(default)s]')
    optional_args.add_argument('-s','--snv-density', dest='SNV_DENSITY', metavar='float', type=float, default=0.1, help='Minimum SNV percentage to consider haplotype blocks [%(default)s]')
    optional_args.add_argument('-c','--cpus', dest='CPUS', metavar='int', type=int, default=1, help='Maximum number of CPUs to be used [%(default)s]')
    # OTHER ARGUMENTS
    other_args = parser.add_argument_group('Other arguments')
    other_args.add_argument('-h','--help', action='help', help='Show this help message and exit')
    other_args.add_argument('-V','--version', action='version', version=f'Strainberry {_version()}', help='Show version number and exit')
    other_args.add_argument('-v','--verbose', dest='VERBOSE', action='count', default=0, help='Verbose output')
    other_args.add_argument('--debug', dest='DEBUG', action='store_true', help=argparse.SUPPRESS)
    opt = parser.parse_args()

    print_status(f'Starting Strainberry v{_version()}')
    if opt.VERBOSE > 0:
        print_status(f'parameters: {"--nanopore" if opt.ONT else ""} -n {opt.STRAINS} -s {opt.SNV_DENSITY} -c {opt.CPUS}')

    validParameters=True
    if not os.path.isfile(opt.FASTA):
        validParameters=False
        print_error(f'FASTA file "{opt.FASTA}" does not exist.')
    elif not os.path.isfile(f'{opt.FASTA}.fai'):
        validParameters=False
        print_error(f'FASTA index file "{opt.FASTA}.fai" does not exist. FASTA file must be indexed with samtools faidx.')
    if not os.path.isfile(opt.BAM):
        validParameters=False
        print_error(f'BAM file "{opt.BAM}" does not exist.')
    elif not os.path.isfile(f'{opt.BAM}.bai'):
        validParameters=False
        print_error(f'BAM file index "{opt.BAM}.bai" does not exist.')
    if os.path.isdir(opt.OUTDIR) and len(os.listdir(opt.OUTDIR))!=0:
        validParameters=False
        print_error(f'output directory "{opt.OUTDIR}" is not empty')
    if not validParameters:
        return 1

    if opt.CPUS < 0:
        opt.CPUS=1

    if opt.STRAINS < 2:
        print_warning(f'-n parameter has been set to 2')
        opt.STRAINS=2

    if not os.path.isdir(opt.OUTDIR):
        os.makedirs(opt.OUTDIR,exist_ok=True)

    nsep=2
    sep_hratio=1.0
    prev_ahr=None
    bamfile=opt.BAM
    fastafile=opt.FASTA
    scaffold_infofile=None
    while nsep <= opt.STRAINS:

        print_status(f'### performing {nsep}-strain separation')
        # Create output directory for current iteration
        out_dir = os.path.join(opt.OUTDIR,f'strainberry_n{nsep}')
        os.makedirs(out_dir,exist_ok=True)

        # SNP calling and phasing
        print_status(f'SNP calling and phasing')
        longshot = LongshotReadSeparator(fastafile,bamfile,out_dir,snv_dens=opt.SNV_DENSITY,nproc=opt.CPUS)
        if not longshot.phase_and_tag():
            print_error(f'haplotype phasing failed')
            return 2
        psc = longshot.psc

        new_ahr = AverageHammingRatio()
        new_hratio,num_phasesets,nreads=new_ahr.compute(bamfile,psc,nproc=opt.CPUS)

        with open(os.path.join(out_dir,'10-separation','reference.ahr.tsv'), 'w') as ofh:
            for ref_id in new_ahr.ahrdict:
                ref_hratio,ref_nreads = new_ahr.reference_average_hratio(ref_id)
                ofh.write(f'{ref_id}\t{ref_hratio:.4f}\t{ref_nreads}\n')

        print_status(f'hratio={new_hratio:.4f} |phasesets|={num_phasesets} nreads={nreads}')
        if num_phasesets == 0:
            print_status(f'no regions to separate')
            print_status(f'best separation available at: {fastafile}')
            break

        # check if current iteration improves average hamming ratio
        if nsep > 2 and (new_hratio > sep_hratio or (sep_hratio-new_hratio) < 0.01*sep_hratio):
            print_status(f'average Hamming ratio did not improve enough: {sep_hratio:.4f} -> {new_hratio:.4f}')
            print_status(f'best separation available at: {fastafile}')
            break
        print_status(f'average Hamming ratio improved: {sep_hratio:.4f} -> {new_hratio:.4f}')
        sep_hratio=new_hratio
        
        # remove phasesets that would not improve hamming ratio of a reference sequence
        if nsep > 2:
            scf_ahr={}
            with open(scaffold_infofile,'r') as ifh:
                for line in ifh:
                    rec=line.split('\t')
                    scf_id=rec[0]
                    ctg_info={ key:val for key,val in (x.split('=') for x in rec[3:]) }
                    if ctg_info['phased']=='true':
                        ref_hratio,ref_nreads=prev_ahr.phaseset_average_hratio(ctg_info['reference'],ctg_info['ps'])
                        scf_hratio,scf_nreads=scf_ahr[scf_id] if scf_id in scf_ahr else (1.0,0)
                        scf_ahr[scf_id] = ( (scf_hratio*scf_nreads+ref_hratio*ref_nreads)/(scf_nreads+ref_nreads), (scf_nreads+ref_nreads) )
            with open(os.path.join(out_dir,'scaffolds.retained.txt'),'w') as ret_fh, open(os.path.join(out_dir,'scaffolds.filtered.txt'),'w') as flt_fh:
                for scf_id in new_ahr.ahrdict:
                    new_scf_hratio,new_scf_nreads = new_ahr.reference_average_hratio(scf_id)
                    pre_scf_hratio,pre_scf_nreads = scf_ahr[scf_id] if scf_id in scf_ahr else (1.0,0)
                    if new_scf_hratio > pre_scf_hratio or (pre_scf_hratio-new_scf_hratio < 0.01*pre_scf_hratio):
                        psc.remove_reference(scf_id)
                        flt_fh.write(f'{scf_id}\t{pre_scf_hratio:.4f}\t{new_scf_hratio:.4f}\n')
                        continue
                    ret_fh.write(f'{scf_id}\t{pre_scf_hratio:.4f}\t{new_scf_hratio:.4f}\n')

        prev_ahr=new_ahr
        
        # read separation
        print_status(f'separating reads')
        if not longshot.separate_reads():
            print_error(f'read separation failed')
            return 2
        reads_dir = longshot.reads_dir
        hap_set = longshot.hap_set

        # print Phaseset statistics
        if opt.DEBUG:
            with open(os.path.join(out_dir,'10-separation','phaseset.stats.tsv'), 'w') as ofh:
                for ref_id in psc.references():
                    for ps in psc.phasesets(ref_id):
                        nreads_ratio = ps.nreads_tagged/ps.nreads_mapped if ps.nreads_mapped > 0 else 0.0
                        record = [ ps.reference, ps.psid, str(ps.start()), str(ps.end()), f'{100.0*ps.density():.2f}',
                                   str(ps.nreads_tagged), str(ps.nreads_mapped), f'{nreads_ratio:.4f}' ]
                        ofh.write('\t'.join(record) + '\n')

        # Haplotype assembly
        print_status(f'assembling strain haplotypes')
        assembly_dir = os.path.join(out_dir,'20-assembly')
        wtdbg = WtdbgAssembler(fastafile,psc,hap_set,reads_dir,assembly_dir,is_ont=opt.ONT,nproc=opt.CPUS)
        wtdbg.run()

        # Scaffolding
        contigfile=wtdbg.contig_file
        scaffolder = Scaffolder(contigfile,bamfile,out_dir,is_ont=opt.ONT,nproc=opt.CPUS)
        if not scaffolder.run():
            print_error(f'scaffolding failed')
            return 2
        print_status(f'scaffold file created at: {scaffolder.scaffold_file}')
        
        # Map input reads to current strain-aware scaffolds
        print_status('mapping reads to strain-separated scaffolds')
        samtools_fastq_cmd = [ 'samtools', 'fastq', bamfile ]
        minimap2_cmd = [ 'minimap2', '-ax', 'map-ont' if opt.ONT else 'map-pb', '-t', f'{opt.CPUS}', scaffolder.scaffold_file, '-' ]
        samtools_sort_cmd = [ 'samtools', 'sort', '--threads', f'{opt.CPUS}', '-o', os.path.join(out_dir,'assembly.scaffolds.bam') ]
        # Run mapping commands
        samtools_fastq = subprocess.Popen(samtools_fastq_cmd, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
        minimap2 = subprocess.Popen(minimap2_cmd, stdin=samtools_fastq.stdout, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
        samtools_sort = subprocess.Popen(samtools_sort_cmd, stdin=minimap2.stdout, stderr=subprocess.DEVNULL)
        samtools_sort.communicate()
        if samtools_sort.returncode != 0:
            print_error(f'read mapping to strain-separated scaffolds failed')
            return 2

        # update file names for next iteration
        bamfile=os.path.join(out_dir,'assembly.scaffolds.bam')
        fastafile=scaffolder.scaffold_file
        scaffold_infofile=scaffolder.scaffold_info_file

        # Index BAM and FASTA files for the next iteration
        if subprocess.run(['samtools','index',f'-@{min(4,opt.CPUS)}',bamfile]).returncode != 0:
            print_error(f'failed to create index for {bamfile}')
            return 2
        if subprocess.run(['samtools','faidx',fastafile]).returncode != 0:
            print_error(f'failed to create index for {fastafile}')
            return 2

        nsep+=1

    # no more iterations to perform: move last separation in the main output directory
    if fastafile != opt.FASTA:
        fastadir,fastaname=os.path.split(fastafile)
        os.rename(fastafile,os.path.join(opt.OUTDIR,fastaname))
        os.rename(f'{fastafile}.fai',os.path.join(opt.OUTDIR,f'{fastaname}.fai'))
        bamdir,bamname=os.path.split(bamfile)
        os.rename(bamfile,os.path.join(opt.OUTDIR,bamname))
        os.rename(f'{bamfile}.bai',os.path.join(opt.OUTDIR,f'{bamname}.bai'))

    print_status('Strainberry finished successfully')
    return 0


if __name__ == "__main__":
    sys.exit(main())
