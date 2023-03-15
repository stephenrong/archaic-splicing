#!/bin/sh

# Get VEP annotations from vcf files
Rscript preprocess_annotations_SNPs.R
sh get_vep_annotations_SNPs.sh

# Join VEP annotations to hub files
Rscript join_vep_annotations_SNPs_v2.R
