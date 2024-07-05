# OPRD1
<!-- badges: start -->
<!-- badges: end -->

This repository contains all code and scripts used to perform the analyses described in the paper "Mirror effect of delta opioid receptor mutations on human metabolic homeostasis".

The repository includes both [bioinformatics](bioinformatics) and [biostatistics](biostatistics) analyses.

**Note: the code has not been tested outside of its original environment and depends on data not publicly available.**

## Bioinformatics

The bioinformatics analysis was performed using the scripts in the [bioinformatics](bioinformatics) directory.

### RNAseq
The different analysis steps are listed in the Nextflow script [bioinformatics/RnaSeqAnalysis.nf](bioinformatics/RnaSeqAnalysis.nf).

Nextflow version used is 20.04.1.5335, R version used is 3.6.3 (February, 2020) on Debian GNU/Linux 9 (stretch, x86_64).

## Biostatistics

The statistical analyses were performed using the scripts in the [biostatistics](biostatistics) directory.  
The directory [biostatistics/utils](biostatistics/utils) contains additional scripts needed in the analyses.

The R version used is the 4.0.2 (June, 2020) on Debian GNU/Linux 9 (stretch, x86_64) available as a Docker image ([umr1283/stat:R402](https://hub.docker.com/r/umr1283/stat)).

## Contact

For questions and other discussion, please contact the authors.
