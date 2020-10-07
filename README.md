![Strainberry logo](https://github.com/rvicedomini/strainberry/raw/master/images/sberry-logo_590x142.png)

# Automated strain separation of low-complexity metagenomes

Strainberry is a method that performs strain separation in low-complexity metagenomes using error-prone long-read technologies. 
It exploits state-of-the-art tools for variant calling, haplotype phasing, and genome assembly, in order to
achieve single-sample assembly of strains with higher quality than other state-of-the-art long-read assemblers.

## System requirements

Strainberry has been developed and tested under a Linux environment.
It requires certain packages/tools in order to be installed/used: 
+ The bash shell
+ [miniconda3](https://conda.io/en/latest/miniconda.html)
+ Standard development packages:
    - Debian/Ubuntu: `build-essential` and `python3-dev`
    - RedHat/CentOS/Fedora: `gcc`, `gcc-c++`, `glibc-devel`, `make`, and `python3-devel`

## Installation

The simplest (and recommended) way to install Strainberry dependencies is by creating an isolated conda environment (*e.g.*, named sberry):
```
git clone https://github.com/rvicedomini/strainberry.git
cd strainberry
conda env create -n sberry --file environment.yml
```

To make the `strainberry` command available, it is advised to include the absolute path of Strainberry's directory in your PATH environment variable by adding the following line to your `~/.bashrc` file:
```
export PATH=/absolute/path/to/strainberry:${PATH}
```

### Updating to latest version

``` 
cd strainberry
git pull
conda env update -n sberry --file environment.yml
```
## Usage

Activate the conda environment:
```
$ conda activate sberry
```

Running Strainberry:
```
$ strainberry [options] -r <FASTA> -b <BAM> -o <OUTPUT_DIR>
```
where `<FASTA>` is a *strain-oblivious* metagenome assembly (*e.g.*, generated with metaFlye) 
and `<BAM>` is a *coordinate-sorted* alignment of long reads in BAM format. 
Strainberry output sequences are be stored in `<OUTPUT_DIR>`

After Strainberry execution the conda environment can be deactivated:
```
$ conda deactivate
```

## Output

The output directory of Strainberry has the following structure:
```
OUTPUT_DIR/
├── 00-preprocess/
├── 10-variants/
├── 20-separation/
├── 30-assembly/
├── 50-sascf/
├── assembly.contigs.fa
└── assembly.scaffolds.fa
```
Strainberry output assemblies are stored in the `assembly.contigs.fa` and `assembly.scaffolds.fa` files.
As usual, the contigs assembly is more conservative while the scaffolds assembly is more contiguous. 
In the contig assembly headers are of the form:
```
>sberry|[reference-name]_[phaseset-id]_h[haplotype]_[contig-index]
```
All the other sub-directories in Strainberry's output contain intermediate results and log files.
Therefore, after a successful run of Strainberry, they could be deleted.

### Test dataset

In order to verify that Strainberry works properly, it is possible to run it on a small dataset in the `example` sub-directory:
```
$ cd example
$ strainberry -r ecoli.fa -b ecoli.sorted.bam -o sberry_out -t 4
```
Strainberry should take around 5 minutes to finish. The input is a small fragment of a consensus *E. coli* sequence.
In the `sberry_out` output directory, both `assembly.contigs.fa` and `assembly.scaffolds.fa` files should contain two sequences
(one closer to strain K12, the other closer to strain W).


## Command line options

```
$ strainberry --help

Strainberry 1.0
Automated strain separation of low-complexity metagenomes

  USAGE:
    strainberry [options] -r <FASTA> -b <BAM> -o <OUTPUT_DIR>

  MANDATORY OPTIONS:         
    -r, --reference <name>   Strain-oblivious assembly in FASTA format
    -b, --bam <name>         Coordinate-sorted read alignment in BAM format
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

