#!/bin/R

# Join VEP annotations to hub files

library(tidyverse)
library(data.table)
library(vcfR)
# source("../get_helper.R")

final_v2_variants_B_stat_mask_1KGP_archaic_gnomAD_hub_vep <- as_tibble(fread("../../../final-analysis-archaic-EndToEnd2/results/vep_annotations_SNPs/final_v2_variants_B_stat_mask_1KGP_archaic_gnomAD_hub_vep.txt.gz"))
final_v2_variants_B_stat_mask_1KGP_archaic_gnomAD_hub_final <- left_join(as_tibble(fread("../../results/preprocess_1KGP_SNPs/final_v2_variants_B_stat_mask_1KGP_archaic_gnomAD_hub.txt.gz")), final_v2_variants_B_stat_mask_1KGP_archaic_gnomAD_hub_vep)
write_tsv(final_v2_variants_B_stat_mask_1KGP_archaic_gnomAD_hub_final, "../../results/vep_annotations_SNPs/final_v2_variants_B_stat_mask_1KGP_archaic_gnomAD_hub_vep.txt.gz")

final_v2_variants_B_stat_mask_1KGP_archaic_gnomAD_hapR2_hub_vep <- as_tibble(fread("../../../final-analysis-archaic-EndToEnd2/results/vep_annotations_SNPs/final_v2_variants_B_stat_mask_1KGP_archaic_gnomAD_hub_vep.txt.gz"))
final_v2_variants_B_stat_mask_1KGP_archaic_gnomAD_hapR2_hub_final <- left_join(as_tibble(fread("../../results/preprocess_1KGP_SNPs/final_v2_variants_B_stat_mask_1KGP_archaic_gnomAD_hapR2_hub.txt.gz")), final_v2_variants_B_stat_mask_1KGP_archaic_gnomAD_hapR2_hub_vep)
write_tsv(final_v2_variants_B_stat_mask_1KGP_archaic_gnomAD_hapR2_hub_final, "../../results/vep_annotations_SNPs/final_v2_variants_B_stat_mask_1KGP_archaic_gnomAD_hapR2_hub_vep.txt.gz")
