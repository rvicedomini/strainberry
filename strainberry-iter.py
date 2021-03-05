#!/usr/bin/env python3

import sys,os,argparse,subprocess
import pysam

from sberry.utils import *
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
    parser = argparse.ArgumentParser(prog='Strainberry', description='Automated strain separation of low-complexity metagenomes', add_help=False)
    # MANDATORY ARGUMENTS
    required_args = parser.add_argument_group('Required arguments')
    required_args.add_argument('-r','--reference', dest='FASTA', metavar='PATH', required=True, help='Strain-oblivious assembly in FASTA format')
    required_args.add_argument('-b','--bam', dest='BAM', metavar='PATH', required=True, help='Read alignment in BAM format')
    required_args.add_argument('-o','--out-dir', dest='OUTDIR', metavar='PATH', required=True, help='Output directory of Strainberry assemblies')
    # OPTIONAL ARGUMENTS
    optional_args = parser.add_argument_group('Optional arguments')
    optional_args.add_argument('-n','--num-strains', dest='STRAINS', metavar='int', type=int, default=5, help='Attempt strain-separation at most for the provided strain multiplicity [%(default)s]')
    optional_args.add_argument('-q','--qual', dest='QUAL', metavar='int', type=int, default=50, help='Consider only variants with a minimum QUAL value [%(default)s]')
    optional_args.add_argument('-s','--snv-density', dest='SNV_DENSITY', metavar='float', type=float, default=0.1, help='Minimum SNV percentage to consider haplotype blocks [%(default)s]')
    optional_args.add_argument('-c','--cpus', dest='CPUS', metavar='int', type=int, default=1, help='Maximum number of CPUs to be used [%(default)s]')
    optional_args.add_argument('--ont','--nanopore', dest='ONT', action='store_true', help='Input consists of Oxford Nanopore raw reads')
    #optional_args.add_argument('--freebayes', dest='FREEBAYES', action='store_true', help='Use freebayes with default parameters instead Longshot for variant calling (SLOWER)')
    other_args = parser.add_argument_group('Other arguments')
    other_args.add_argument('-h','--help', action='help', help='Show this help message and exit')
    other_args.add_argument('-V','--version', action='version', version=f'Strainberry {_version()}', help='Show version number and exit')
    other_args.add_argument('-v','--verbose', dest='VERBOSE', action='count', default=0, help='Verbose output')
    opt = parser.parse_args()

    validParameters=True
    if not os.path.isfile(opt.FASTA):
        validParameters=False
        print_error(f'FASTA file "{opt.FASTA}" does not exist.')
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

    print_status(f'Starting Strainberry v{_version()}')
    if opt.VERBOSE > 0:
        print_status(f'parameters: -n {opt.STRAINS} -q {opt.QUAL} -s {opt.SNV_DENSITY} -c {opt.CPUS} {opt.ONT}')

    nsep=2
    bamfile=opt.BAM
    fastafile=opt.FASTA

    while nsep <= opt.STRAINS:

        print_status(f'### Strainberry {nsep}-strain separation')

        # Create directory for current iteration
        out_dir = os.path.join(opt.OUTDIR,f'strainberry_n{nsep}')
        os.makedirs(out_dir,exist_ok=True)

        # SNP calling, phasing, and read separation
        sberry_variants_cmd = [ os.path.join(sberry_root,'sberry-variants'), '-r', fastafile, '-b', bamfile, '-o', out_dir, '-q', str(opt.QUAL) ]
        sberry_variants = subprocess.run(sberry_variants_cmd)
        if sberry_variants.returncode != 0:
            print_error(f'sberry-variants command failed')
            return 2

        # Haplotype assembly
        sberry_assemble_cmd = [ os.path.join(sberry_root,'sberry-assemble'),
                '-r', os.path.join(out_dir,'00-preprocess','reference.fa'),
                '-b', os.path.join(out_dir,'20-separation','tagged'),
                '-v', os.path.join(out_dir,'20-separation','variants.phased.vcf.gz'),
                '-s', str(opt.SNV_DENSITY),
                '-o', out_dir,
                '-t', str(opt.CPUS),
                '--nanopore' if opt.ONT else '' ]
        sberry_assemble = subprocess.run(sberry_assemble_cmd)
        if sberry_assemble.returncode != 0:
            print_error(f'sberry-assemble command failed')
            return 2

        # Scaffolding
        sberry_sascf_cmd = [ os.path.join(sberry_root,'sberry-sascf'),
                '-r', os.path.join(out_dir,'assembly.contigs.fa'),
                '-b', os.path.join(out_dir,'00-preprocess','alignment.bam'),
                '-o', out_dir,
                '-t', str(opt.CPUS),
                '--nanopore' if opt.ONT else '' ]
        sberry_sascf = subprocess.run(sberry_sascf_cmd)
        if sberry_sascf.returncode != 0:
            print_error(f'sberry-sascf command failed')
            return 2
        
        # update file names for next iteration
        bamfile=os.path.join(out_dir,'assembly.scaffolds.bam')
        fastafile=os.path.join(out_dir,'assembly.scaffolds.fa')

        # Map input reads to current strain-aware scaffolds
        print_status('Mapping reads to Strainberry scaffolds')
        samtools_fastq_cmd = [ 'samtools', 'fastq', os.path.join(out_dir,'00-preprocess','alignment.bam') ]
        minimap2_cmd = [ 'minimap2', '-ax', 'map-ont' if opt.ONT else 'map-pb', '-t', f'{opt.CPUS}', os.path.join(out_dir,'assembly.scaffolds.fa'), '-' ]
        samtools_sort_cmd = [ 'samtools', 'sort', '--threads', f'{opt.CPUS}', '-o', os.path.join(out_dir,'assembly.scaffolds.bam') ]
        samtools_index_cmd = [ 'samtools', 'index', f'-@{min(2,opt.CPUS)}', os.path.join(out_dir,'assembly.scaffolds.bam') ]
        # Run mapping commands
        samtools_fastq = subprocess.Popen(samtools_fastq_cmd, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
        minimap2 = subprocess.Popen(minimap2_cmd, stdin=samtools_fastq.stdout, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
        samtools_sort = subprocess.Popen(samtools_sort_cmd, stdin=minimap2.stdout, stderr=subprocess.DEVNULL)
        samtools_sort.communicate()
        if samtools_sort.returncode != 0:
            print_error(f'read mapping to strain-separated scaffolds failed')
            return 2
        samtools_index = subprocess.run(samtools_index_cmd)
        if samtools_index.returncode != 0:
            print_error(f'failed to create index {bamfile}')
            return 2

        nsep+=1

    print_status('Strainberry finished successfully')
    return 0


if __name__ == "__main__":
    sys.exit(main())
