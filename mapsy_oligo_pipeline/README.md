# Neanderthal Updated Library

Updated: 2022-10-18. Contact: Stephen Rong.


### Introduction

This repository contains the scripts for generating the oligo sequences for the Neanderthal MaPSy experiments.


### Contents

**scripts/**: Contains the scripts for generating oligo sequences.
  * master_script.sh: Runs individual scripts to generate all results.
  * mapsy_helper_v1.R: Helper functions for automating oligo sequence design.
  * mapsy_create_v1.R: Command line script which takes as input a fasta reference genome, a gff transcript annotation, a vcf variant file, a param file, and a list of gene names to filter, and generates oligos.
  * mapsy_library_2_params.R: Param file for generating variable part of oligos.
  * mapsy_library_2_scripts_archaic.sh: Script file for generating variable part of oligos.
  * mapsy_library_2_finish_archaic.R: Add primer and common sequences to create final oligos.
  * for_analysis_archaic.R: Create input and output reference genomes for alignment and variant to sequence association tables.

**results/**: 
  * mapsy_libraries/: Contains intermediate files from converting variants to oligo sequences.
  * mapsy_orders/: Contains final files after adding primer and common sequences to oligos.
  * for_analysis/: Contains input and output reference genomes for alignment and variant to sequence association tables.

**data/**: Contains the input data needed for running the scripts.
  * commons/: Contains common half exon sequences.
  * fasta/: Contains human reference genome sequence.
  * gff/: Contains gene annotations files.
  * merge_variants/: Contains VCF of variants for Neanderthal updated library.
  * primers/: Contains primer pair sequences.


### Dependencies

R (version 4.0.2), with packages: optparse (1.6.6), tidyverse (1.3.0), data.table (1.13.4), plyranges (4.0.0), rtracklayer (4.0.2), BSgenome (1.56.0), GenomicFeatures (1.40.1), VariantAnnotation (1.34.0), stringdist (0.9.6.3), stringi (1.5.3).
