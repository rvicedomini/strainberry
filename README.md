![Strainberry logo](https://github.com/rvicedomini/strainberry/raw/master/images/sberry-logo_590x142.png)

# Automated strain separation of low-complexity metagenomes

Strainberry is a method that performs strain separation in low-complexity metagenomes using error-prone long-read technologies. 
It exploits state-of-the-art tools for variant calling, haplotype phasing, and genome assembly, in order to
achieve single-sample assembly of strains with higher quality than other state-of-the-art long-read assemblers.

# System requirements

+ Strainberry has been developed and tested under a Linux environment
+ The bash shell should be installed
+ To easily install Strainberry dependencies, [miniconda3](https://conda.io/en/latest/miniconda.html) is strongly recommended

### Installation

The recommended way to install Strainberry is through an isolated environment built with conda (*e.g.*, sberry):
```
git clone https://github.com/rvicedomini/strainberry.git
cd strainberry
conda env create -n sberry --file environment.yml
```

It is advised to include Strainberry directory in your PATH environment variable by adding the following line to your `~/.bashrc` file:
```
export PATH=/path/to/strainberry:${PATH}
```

# Usage

Activate the conda environment:
```
$ conda activate sberry
```

Running Strainberry:
```
$ strainberry [options] -r <FASTA> -b <BAM> -o <OUTPUT_DIR>
```

After Strainberry execution the conda environment can be deactivated:
```
$ conda deactivate sberry
```

### Command line options

```
$ strainberry --help

Strainberry 1.0
Automated strain separation of low-complexity metagenomes

  USAGE:
    strainberry [options] -r <FASTA> -b <BAM> -o <OUTPUT_DIR>

  MANDATORY OPTIONS:         
    -r, --reference <name>   Input (strain-oblivious) assembly in FASTA format
    -b, --bam <name>         Input read alignment in BAM format
                             (at least a 60X coverage is recommended)
    -o, --output-dir <name>  Output directory

  OTHER OPTIONS:            
    --nanopore              Input consists of Oxford Nanopore raw reads
    -n, --strains <num>     Strain multiplicity [2]
    -q, --qual <num>        Consider only variants with a minimum QUAL value [50]
    -s, --snv-dens <float>  Minimum SNV percentage to consider haplotype blocks [0.1]
                            It should be a value in [0,100]
    -t, --threads <num>     Number of threads [1]
    --freebayes             Use exclusively freebayes with default parameters instead
                            of Longshot for variant calling (SLOWER)

    -h, --help     Print this help message
    -V, --version  Print version
```

##

### conda packages
- Python 3.7+
- wtdbg2
- GNU parallel
- vcflib
- longshot
- freebayes
- whatshap polyploid (polyploid-haplotag branch)
- minimap2
- mummer4

### python modules
- biopython
- pysam
- networkx
- pygraphviz
