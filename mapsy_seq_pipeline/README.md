# Neanderthal MaPSy analysis

Updated: 2022-10-18. Contact: Stephen Rong.


### Introduction

This repository contains scripts for analyzing NGS data from the Neanderthal MaPSy experiments.


### Contents

**scripts/**: 
  * star_index_Neanderthal_exonic_updated.sh: 
  * sequencing_helper_files/: 
    * star_all_standard_neanderthal.sh: STAR alignment for Neanderthal experiments.
    * star_all_standard_process_neanderthal.sh: Postprocess STAR alignment for Neanderthal experiments.
  * sequencing_Neanderthal_updated/: 
    * master_Neanderthal_updated_exonic.sh: Runs all scripts in this subdirectory.
    * star_all_standard_Neanderthal_updated_exonic.sh: Runs star_all_standard_neanderthal.sh on each paired-end sequencing dataset.
    * star_all_standard_process_Neanderthal_updated_exonic.sh: Runs star_all_standard_process_neanderthal.sh on each paired-end sequencing dataset.
  * mpralm_splicing_scores/: 
    * mpralm_analyze_exp.R: Computes mpralm scores from MaPSy experiment.
  * postprocess_Neanderthal_updated/: 
    * postprocess_Neanderthal_updated.R: MaPSy analysis on Neanderthal updated library.

**results/**: 
  * postprocess_Neanderthal_updated/: Outputs of postprocess_Neanderthal_updated.R.

**data/**:
  * custom_reference/: Custom reference genomes for STAR alignment.
  * custom_star_index/: STAR index files for STAR alignment.
  * known_canonical/: GRanges files for canonical transcript exons.
  * Neanderthal-HEK_293T-final_expanded_exonic/: 
    * fastqc/: Folder for FASTQC results.
      * multiqc/: MultiQC results for FASTQ files.
    * star_align/: Folder for STAR alignment BAM files.
      * multiqc/: MultiQC results for BAM files.
    * star_idxstats/: Counts of uniquely mapping reads.
    * star_sort/: Temp folder for sorting BAM files.
    * star_unique/: Temp folder for BAM uniquely mapping reads.
  * premapsy_variants/: VCF files for Neanderthal updated variants.


### Dependencies

R (version 4.0.2), with packages: tidyverse (1.3.0), ggpubr (0.4.0), GGally (2.0.0), DescTools (0.99.39), mpra (1.10.0). Other Software: bedtools (2.29.1), samtools (1.9), STAR (2.7.3a), fastqc (0.11.8), multiqc (1.8).
