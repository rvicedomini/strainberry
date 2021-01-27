![Strainberry logo](https://github.com/rvicedomini/strainberry/raw/master/images/sberry-logo_590x142.png)

# Automated strain separation of low-complexity metagenomes

Strainberry is a method that performs strain separation in low-complexity metagenomes using error-prone long-read technologies. 
It exploits state-of-the-art tools for variant calling, haplotype phasing, and genome assembly, in order to
achieve single-sample assembly of strains with higher quality than other state-of-the-art long-read assemblers.

+ [System requirements](#system-requirements)
+ [Installation](#installation)
+ [Usage](#usage)
+ [Output files](#output-files)
+ [Example](#example)
+ [Command-line options](#command-line-options)
+ [Results](#results)

## System requirements

Strainberry has been developed and tested under a Linux environment.
It requires certain packages/tools in order to be installed/used: 
+ GNU bash (version 4 or later recommended)
+ [miniconda3](https://conda.io/en/latest/miniconda.html)
+ Standard development packages with a GCC version supporting C++11:
    - Debian/Ubuntu: `build-essential` and `python3-dev`
    - RedHat/CentOS/Fedora: `gcc`, `gcc-c++`, `glibc-devel`, `make`, and `python3-devel`

## Installation

The simplest (and recommended) way to install Strainberry dependencies is by creating an isolated conda environment (*e.g.*, named sberry):
```bash
git clone https://github.com/rvicedomini/strainberry.git
cd strainberry
conda env create -n sberry --file environment.yml
```
The whole installation process should take about 5-10 minutes.

To make the `strainberry` command available, it is also advised to include the absolute path of Strainberry's directory in your PATH environment variable by adding the following line to your `~/.bashrc` file:

```
export PATH=/absolute/path/to/strainberry:${PATH}
```

### Updating to the latest version

```bash
cd strainberry
git pull
conda env update -n sberry --file environment.yml
```

## Usage

Activate the conda environment:
```
conda activate sberry
```

Running Strainberry:
```
strainberry [options] -r <FASTA> -b <BAM> -o <OUTPUT_DIR>
```

where `<FASTA>` is a *strain-oblivious* metagenome assembly (*e.g.*, generated with metaFlye) 
and `<BAM>` is a *coordinate-sorted* long-read alignment in BAM format. 
Strainberry output sequences are stored in `<OUTPUT_DIR>`

After Strainberry execution the conda environment can be deactivated with the command:
```
conda deactivate
```

## Output files

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
In the contig assembly headers are of the form
```
>sberry|[reference-id]_[phaseset-id]_h[haplotype-id]_[contig-id]
```
for the separated haplotypes/strains, while are of the form:
```
>[reference-id]_[reference-position]
```
for the non-separated regions of the input reference assembly.

All the sub-directories in Strainberry's output contain intermediate results and log files.
Therefore, after a successful run of Strainberry, they could be deleted.

## Example

In order to verify that Strainberry works properly, it is possible to run it on a small dataset in the `example` sub-directory.

### Generating the input from the reads (optional)
In order to generate a strain-oblivious assembly and a read alignment, we recommend to use metaFlye, minimap2, and samtools.
Assuming these tools are available, it is possible to run the following commands, using 12 threads:

```
cd example
flye --meta --pacbio-raw reads.fq.gz --out-dir flye_out --genome-size 300k --threads 12
minimap2 -ax map-pb -t 12 ./flye_out/assembly.fasta reads.fq.gz | samtools sort >./flye_out/alignment.sorted.bam
```

where `--genome-size` provides an estimate of the metagenome to metaFlye (not required from version 2.8).
The assembly and read alignment are then available in the `flye_out` directory as `assembly.fasta` and `alignment.sorted.bam` respectively.

### Running Strainberry
Given a strain-oblivious assembly (file `assembly.fasta`) and a long-read mapping (file `alignment.sorted.bam`), it is possible to run Strainberry using 4 threads as follows:

```
$ cd example
$ strainberry -r assembly.fasta -b alignment.sorted.bam -o sberry_out -t 4
```

Strainberry should take around 5 minutes to finish. The file `assembly.fasta` contains a single sequence which is a consensus of a small region of *E. coli* strains K12 and W.
After a successful run of Strainberry, in the `sberry_out` directory, the file `assembly.scaffolds.fa` should contain two scaffolds (one closer to strain K12, the other closer to strain W).

## Command line options

```
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

## Results

A manuscript presenting Strainberry has been submitted and is currently under review.
A pre-print will be made available soon.

The analysis workflow and all the scripts necessary to reproduce the main results of Strainberry are available at [https://github.com/rvicedomini/strainberry-analyses](https://github.com/rvicedomini/strainberry-analyses).
