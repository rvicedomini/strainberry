![Strainberry logo](https://github.com/rvicedomini/strainberry/raw/master/images/sberry-logo_590x142.png)

# Automated strain separation of low-complexity metagenomes

WORK IN PROGRESS

## requirements
- Linux environment with Bash
- Conda (to install following software)

### software
- Python 3.7+
- Flye 2.7
- GNU parallel
- tabix, vcftools, vcflib (verify precisely)
- longshot
- freebayes
- whatshap polyploid (polyploid_haplotag branch)
- minimap2
- mummer4 (this one can be optional if I manage to integrate minimap2 in the main pipeline and leave mummer only for evaluation)

### python modules
- biopython
- pysam
- networkx
- pygraphviz
