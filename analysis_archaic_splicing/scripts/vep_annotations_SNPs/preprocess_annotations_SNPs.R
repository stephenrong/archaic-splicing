#!/bin/R

# Preprocess hub files to vcf files

library(tidyverse)
library(data.table)
source("../get_helper.R")

# # output to VCF
final_v2_variants_B_stat_mask_1KGP_archaic_gnomAD_hub <- as_tibble(fread(
	"../../results/preprocess_1KGP_SNPs/final_v2_variants_B_stat_mask_1KGP_archaic_gnomAD_hub.txt.gz"))
hub_to_vcf(final_v2_variants_B_stat_mask_1KGP_archaic_gnomAD_hub, 
	"final_v2_variants_B_stat_mask_1KGP_archaic_gnomAD_hub", 
	"../../results/vep_annotations_SNPs/final_v2_variants_B_stat_mask_1KGP_archaic_gnomAD_hub.vcf", 
	BSgenome.Hsapiens.UCSC.hg19, partial=FALSE, compress=TRUE)

ALL_1KGP_phase3_MAF0.01_hub <- as_tibble(fread(
	"../../results/preprocess_1KGP_SNPs/ALL_1KGP_phase3_MAF0.01_hub.txt.gz"))
hub_to_vcf(ALL_1KGP_phase3_MAF0.01_hub, 
	"ALL_1KGP_phase3_MAF0.01_hub", 
	"../../results/vep_annotations_SNPs/ALL_1KGP_phase3_MAF0.01_hub.vcf", 
	BSgenome.Hsapiens.UCSC.hg19, partial=FALSE, compress=TRUE)
