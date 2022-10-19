#!/bin/R

library(tidyverse)
library(data.table)

source("../get_helper.R")

# load table
final_v2_variants_B_stat_mask_1KGP_archaic_gnomAD_introgr_hub <- as_tibble(fread("../../results/preprocess_1KGP_SNPs/final_v2_variants_B_stat_mask_1KGP_archaic_gnomAD_introgr_hub.txt.gz"))

# output bed
hub_to_bed(final_v2_variants_B_stat_mask_1KGP_archaic_gnomAD_introgr_hub, "../../results/annotate_splice_prediction/final_v2_variants_B_stat_mask_1KGP_archaic_gnomAD_introgr_hub.bed")
