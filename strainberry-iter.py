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
    if git.returncode != 0:
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
    other_args.add_argument('-v','--verbose', dest='VERBOSE', action='count', default=0, help='Verbose output')
    other_args.add_argument('-V','--version', action='version', version=f'Strainberry {_version()}')
    opt = parser.parse_args()

    validParameters=True
    if not os.path.isfile(opt.FASTA):
        validParameters=False
        print_error(f'FASTA file "{opt.FASTA}" does not exist.')
    if not os.path.isfile(opt.BAM):
        validParameters=False
        print_error(f'BAM file "{opt.BAM}" does not exist.')
    if os.path.isdir(opt.OUTDIR) and len(os.listdir(opt.OUTDIR))!=0:
        validParameters=False
        print_error(f'Output directory "{opt.OUTDIR}" must be empty')
    if not validParameters:
        return 1

    if opt.CPUS < 0:
        opt.CPUS=1

    if not os.path.isdir(opt.OUTDIR):
        os.makedirs(opt.OUTDIR,exist_ok=True)

    print_status(f'Starting Strainberry v{_version()}')
    if opt.VERBOSE > 0:
        print_status(f'parameters: -n {opt.STRAINS} -q {opt.QUAL} -s {opt.SNV_DENSITY} -c {opt.CPUS} {opt.ONT}')
    
    # SNP calling, phasing, and read separation
    sberry_variants_cmd = [ os.path.join(sberry_root,'sberry-variants'), '-r', opt.FASTA, '-b', opt.BAM, '-o', opt.OUTDIR, '-q', str(opt.QUAL) ]
    sberry_variants = subprocess.run(sberry_variants_cmd)
    if sberry_variants.returncode != 0:
        print_error(f'sberry-variants command failed')
        return 2

    # Haplotype assembly
    sberry_assemble_cmd = [ os.path.join(sberry_root,'sberry-assemble'),
            '-r', os.path.join(opt.OUTDIR,'00-preprocess','reference.fa'),
            '-b', os.path.join(opt.OUTDIR,'20-separation','tagged'),
            '-v', os.path.join(opt.OUTDIR,'20-separation','variants.phased.vcf.gz'),
            '-s', str(opt.SNV_DENSITY),
            '-o', opt.OUTDIR,
            '-t', str(opt.CPUS),
            '--nanopore' if opt.ONT else '' ]
    sberry_assemble = subprocess.run(sberry_assemble_cmd)
    if sberry_assemble.returncode != 0:
        print_error(f'sberry-assemble command failed')
        return 2

    # Scaffolding
    sberry_sascf_cmd = [ os.path.join(sberry_root,'sberry-sascf'),
            '-r', os.path.join(opt.OUTDIR,'assembly.contigs.fa'),
            '-b', os.path.join(opt.OUTDIR,'00-preprocess','alignment.bam'),
            '-o', opt.OUTDIR,
            '-t', str(opt.CPUS),
            '--nanopore' if opt.ONT else '' ]
    sberry_sascf = subprocess.run(sberry_sascf_cmd)
    if sberry_sascf.returncode != 0:
        print_error(f'sberry-sascf command failed')
        return 2

    print_status('Strainberry finished successfully')
    return 0


if __name__ == "__main__":
    sys.exit(main())
