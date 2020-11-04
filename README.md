# OPRD1
<!-- badges: start -->
<!-- badges: end -->

This repository contains all code and scripts used to perform the analyses described in the paper "... OPRD1 ...".

The repository includes both [bioinformatics](bioinformatics) and [biostatistics](biostatistics) analyses.

**Note: the code has not been tested outside of its original environment and depends on data not publicly available.**

## Bioinformatics

The bioinformatics analysis was performed using the scripts in the [bioinformatics](bioinformatics) directory.

# DNAseq
The different analysis steps are listed in the plain text file [bioinformatics/template_analysis](bioinformatics/template_analysis), which is provided to the [bioinformatics/start_multiserveur.pl](bioinformatics/start_multiserveur.pl) script.  
Then the [bioinformatics/start_sample.pl](bioinformatics/start_sample.pl) script is launched with all the different steps of the analysis.

All the different scripts to perform each steps can be found in the [bioinformatics/Include](bioinformatics/Include) directory.  
These scripts use local databases:  
- MySQL databases with informations about the different samples (*e.g.*, sample ID, sex, demultiplexing index and lane, raw data path) and the dbSNP annotations (version 135).  
- A MongoDB database with the dbNSFP annotations (version 3.0).  
- A Redis database with mutation counting informations.

Tools, captures and sequencing adaptors paths can be passed by the files in the [bioinformatics/Config](bioinformatics/Config) directory.

Perl version used is Perl 5, version 16, subversion 3 (v5.16.3) on CentOS Linux 7 (Core, x86_64).

# RNAseq
The different analysis steps are listed in the Nextflow script [bioinformatics/RnaSeqAnalysis.nf](bioinformatics/RnaSeqAnalysis.nf). A supplementary script [bioinformatics/rnaseq_add_genename.R](bioinformatics/rnaseq_add_genename.R) called during the analysis is also included.

Nextflow version used is 20.04.1.5335, R version used is 3.6.3 (February, 2020) on Debian GNU/Linux 9 (stretch, x86_64).

## Biostatistics

The statistical analyses were performed using the scripts in the [biostatistics](biostatistics) directory.  
The directory [biostatistics/utils](biostatistics/utils) contains additional scripts needed in the analyses.

The R version used is the 4.0.2 (June, 2020) on Debian GNU/Linux 9 (stretch, x86_64) available as a Docker image ([umr1283/stat:R402](https://hub.docker.com/r/umr1283/stat)).

## Contact

For questions and other discussion, please contact the authors.
