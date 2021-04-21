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
+ [Reference](#reference)

## System requirements

Strainberry has been developed and tested under a Linux environment.
It requires certain packages/tools in order to be installed/used: 
+ GNU bash (version 4 or later recommended)
+ [miniconda3](https://conda.io/en/latest/miniconda.html)
+ Standard development packages with a GCC version supporting C++11:
    - Debian/Ubuntu: `build-essential` and `python3-dev`
    - RedHat/CentOS/Fedora: `gcc`, `gcc-c++`, `glibc-devel`, `make`, and `python3-devel`

## Installation

The simplest (and recommended) way to install Strainberry dependencies is by creating an isolated conda environment (*e.g.*, named `sberry`):
```bash
git clone https://github.com/rvicedomini/strainberry.git
cd strainberry
conda env create -n sberry --file environment.yml
```
The whole installation process should take about 5-10 minutes.

To make the `strainberry` command available, it is advised to include the absolute path of Strainberry's directory in your PATH environment variable by adding the following line to your `~/.bashrc` file:

```
export PATH=/absolute/path/to/strainberry:${PATH}
```

### Updating to the latest version

Assuming Strainberry conda environment has `sberry` name, the following commands allow to update Strainberry to the latest (development) version.
The `git pull` command is needed only if Strainberry has been cloned from this git repository.
```bash
cd strainberry
git pull
conda env update -n sberry --file environment.yml --prune
```

## Usage

Activate Strainberry conda environment:
```
conda activate sberry
```

Running Strainberry:
```
strainberry [options] -r FASTA -b BAM -o OUTPUT_DIR
```

where `FASTA` is a *strain-oblivious* metagenome assembly (*e.g.*, generated with metaFlye) and `BAM` is a *coordinate-sorted* long-read alignment in BAM format.
Both `FASTA` and `BAM` files are expected to be indexed with the `samtools faidx` and `samtools index` commands, respectively.
Strainberry's output is stored in `OUTPUT_DIR`.

After Strainberry execution the conda environment can be deactivated with the command:
```
conda deactivate
```

## Output files

The output directory of Strainberry has the following structure:

```
OUTPUT_DIR/
├── strainberry_n2/
├── strainberry_n3/
├── ...
├── strainberry_nK/
├── assembly.scaffolds.bam
├── assembly.scaffolds.bam.bai
├── assembly.scaffolds.fa
└── assembly.scaffolds.fa.fai
```

Strainberry output assembly is stored in the `assembly.scaffolds.fa` file.
A minimap2-based alignment of input reads on the output assembly is also available in the `assembly.scaffolds.bam` file.

All sub-directories with prefix `strainberry_n{k}` contain intermediate results of Strainberry iterations and log files.
After a successful run of Strainberry, they could be deleted.

## Example

In order to verify that Strainberry has been correctly installed, it is possible to test it on a small dataset in the `example` sub-directory.

### Generating the input from the reads (optional)
In order to generate a strain-oblivious assembly and a read alignment, we recommend to use metaFlye, minimap2, and samtools.
Assuming these tools are available, it is possible to run the following commands, using 12 threads:
```
cd example
flye --meta --pacbio-raw reads.fq.gz --out-dir flye_out --genome-size 300k --threads 12
minimap2 -ax map-pb -t 12 ./flye_out/assembly.fasta reads.fq.gz | samtools sort >./flye_out/alignment.sorted.bam
samtools faidx ./flye_out/assembly.fasta
samtools index ./flye_out/alignment.sorted.bam
```
where `--genome-size` provides an estimate of the metagenome to metaFlye (not required from version 2.8).
The assembly and read alignment are then available in the `flye_out` directory as `assembly.fasta` and `alignment.sorted.bam` respectively.

### Running Strainberry
Given a strain-oblivious assembly (file `assembly.fasta`) and a long-read mapping (file `alignment.sorted.bam`), it is possible to run Strainberry using 4 threads as follows:

```
$ cd example
$ strainberry -r assembly.fasta -b alignment.sorted.bam -o sberry_out -c 4
```

Strainberry should take around 5 minutes to finish. The file `assembly.fasta` contains a single sequence which is a consensus of a small region of *E. coli* strains K12 and W.
After a successful run of Strainberry, in the `sberry_out` directory, the file `assembly.scaffolds.fa` should contain two scaffolds (one closer to strain K12, the other closer to strain W).

## Command line options

```
usage: strainberry -r PATH -b PATH -o PATH [--nanopore] [-n int] [-s float]
                   [-c int] [-h] [-V] [-v]

Automated strain separation of low-complexity metagenomes

Required arguments:
  -r PATH, --reference PATH
                        Strain-oblivious assembly in FASTA format
  -b PATH, --bam PATH   Read alignment in BAM format
  -o PATH, --out-dir PATH
                        Output directory of Strainberry assemblies

Optional arguments:
  --nanopore            Input consists of Oxford Nanopore raw reads
  -n int, --max-strains int
                        Attempt strain-separation at most for the provided
                        strain multiplicity [5]
  -s float, --snv-density float
                        Minimum SNV percentage to consider haplotype blocks
                        [0.1]
  -c int, --cpus int    Maximum number of CPUs to be used [1]

Other arguments:
  -h, --help            Show this help message and exit
  -V, --version         Show version number and exit
  -v, --verbose         Verbose output
```

## Reference

A manuscript presenting Strainberry has been submitted and is currently under review.
A pre-print is available on bioRxiv:

R. Vicedomini, C. Quince, A. E. Darling, R. Chikhi, 
"Automated strain separation in low-complexity metagenomes using long reads."
*bioRxiv*, 2021. [https://doi.org/10.1101/2021.02.24.429166](https://doi.org/10.1101/2021.02.24.429166)

