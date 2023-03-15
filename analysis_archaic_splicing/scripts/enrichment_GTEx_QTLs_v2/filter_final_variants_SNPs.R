#!/bin/R

# Filter final library SNP info, filter AI, filter to MAF >= 1% in EUR

library(tidyverse)
library(data.table)

# fix files
final_v2_variants_B_stat_mask_1KGP_archaic_gnomAD_hapR2_hub <- as_tibble(fread("../../results/preprocess_1KGP_SNPs/final_v2_variants_B_stat_mask_1KGP_archaic_gnomAD_hapR2_hub.txt.gz")) %>% 
	mutate(hub_variant_CHROM=as.character(hub_variant_CHROM)) %>% 
	mutate(CHROM=as.character(CHROM))

final_v2_variants_B_stat_mask_1KGP_archaic_gnomAD_introgr_hub <- as_tibble(fread("../../results/preprocess_1KGP_SNPs/final_v2_variants_B_stat_mask_1KGP_archaic_gnomAD_introgr_hub.txt.gz")) %>% 
	mutate(hub_variant_CHROM=as.character(hub_variant_CHROM)) %>% 
	mutate(CHROM=as.character(CHROM))

final_v2_variants_B_stat_mask_1KGP_archaic_gnomAD_introgr_hapR2_hub <- bind_cols(final_v2_variants_B_stat_mask_1KGP_archaic_gnomAD_introgr_hub, final_v2_variants_B_stat_mask_1KGP_archaic_gnomAD_hapR2_hub %>% 
	dplyr::select(hub_variant_ID, EUR_AC_bin, EUR_LDtagN, EUR_LDtagN_bin))
write_tsv(final_v2_variants_B_stat_mask_1KGP_archaic_gnomAD_introgr_hapR2_hub, gzfile("../../results/preprocess_1KGP_SNPs/final_v2_variants_B_stat_mask_1KGP_archaic_gnomAD_introgr_hapR2_hub.txt.gz"))

# load files
final_v2_variants_B_stat_mask_1KGP_archaic_gnomAD_introgr_hapR2_hub <- as_tibble(fread("../../results/preprocess_1KGP_SNPs/final_v2_variants_B_stat_mask_1KGP_archaic_gnomAD_introgr_hapR2_hub.txt.gz")) %>% 
	mutate(hub_variant_CHROM=as.character(hub_variant_CHROM)) %>% mutate(CHROM=as.character(CHROM))

# filter variants
filter_final_v2_variants_B_stat_mask_1KGP_archaic_gnomAD_introgr_hapR2_hub <- final_v2_variants_B_stat_mask_1KGP_archaic_gnomAD_introgr_hapR2_hub %>% 
	filter(hub_variant_REF %in% c("A", "C", "T", "G")) %>%  # SNP only
	filter(hub_variant_ALT %in% c("A", "C", "T", "G")) %>%  # SNP only
	filter(!(hub_variant_CHROM %in% c("X", "Y")))  # autosomes only

# filter by European MAF
filter_final_v2_variants_B_stat_mask_1KGP_archaic_gnomAD_introgr_hapR2_hub <- filter_final_v2_variants_B_stat_mask_1KGP_archaic_gnomAD_introgr_hapR2_hub %>% 
	filter(EUR_AF >= 0.01, EUR_AF <= 0.99)
write_tsv(filter_final_v2_variants_B_stat_mask_1KGP_archaic_gnomAD_introgr_hapR2_hub, gzfile("../../results/enrichment_GTEx_QTLs_v2/filter_final_v2_variants_B_stat_mask_1KGP_archaic_gnomAD_introgr_hapR2_hub.txt.gz"))

# save ldvernot akey 2016
filter_final_v2_variants_B_stat_mask_1KGP_archaic_gnomAD_introgr_hapR2_ldvernot_akey_2016_hub <- filter_final_v2_variants_B_stat_mask_1KGP_archaic_gnomAD_introgr_hapR2_hub %>% 
	filter(hub_in_ldvernot_akey_2016_EUR)
write_tsv(filter_final_v2_variants_B_stat_mask_1KGP_archaic_gnomAD_introgr_hapR2_ldvernot_akey_2016_hub, gzfile("../../results/enrichment_GTEx_QTLs_v2/filter_final_v2_variants_B_stat_mask_1KGP_archaic_gnomAD_introgr_hapR2_ldvernot_akey_2016_hub.txt.gz"))

# save browning 2016
filter_final_v2_variants_B_stat_mask_1KGP_archaic_gnomAD_introgr_hapR2_browning_2018_hub <- filter_final_v2_variants_B_stat_mask_1KGP_archaic_gnomAD_introgr_hapR2_hub %>% 
	filter(hub_in_browning_2018_CEU|hub_in_browning_2018_IBS|hub_in_browning_2018_GBR|hub_in_browning_2018_FIN|hub_in_browning_2018_TSI)
write_tsv(filter_final_v2_variants_B_stat_mask_1KGP_archaic_gnomAD_introgr_hapR2_browning_2018_hub, gzfile("../../results/enrichment_GTEx_QTLs_v2/filter_final_v2_variants_B_stat_mask_1KGP_archaic_gnomAD_introgr_hapR2_browning_2018_hub.txt.gz"))

# save gittelman 2016
filter_final_v2_variants_B_stat_mask_1KGP_archaic_gnomAD_introgr_hapR2_gittelman_2016_hub <- filter_final_v2_variants_B_stat_mask_1KGP_archaic_gnomAD_introgr_hapR2_hub %>% 
	filter(hub_in_gittelman_2016_EUR)
write_tsv(filter_final_v2_variants_B_stat_mask_1KGP_archaic_gnomAD_introgr_hapR2_gittelman_2016_hub, gzfile("../../results/enrichment_GTEx_QTLs_v2/filter_final_v2_variants_B_stat_mask_1KGP_archaic_gnomAD_introgr_hapR2_gittelman_2016_hub.txt.gz"))

# save archaic introgressed all
filter_final_v2_variants_B_stat_mask_1KGP_archaic_gnomAD_introgr_hapR2_introgressed_hub <- filter_final_v2_variants_B_stat_mask_1KGP_archaic_gnomAD_introgr_hapR2_hub %>% 
	filter(hub_in_ldvernot_akey_2016_EUR|hub_in_browning_2018_CEU|hub_in_browning_2018_IBS|hub_in_browning_2018_GBR|hub_in_browning_2018_FIN|hub_in_browning_2018_TSI, hub_in_final_study_introgressed)
write_tsv(filter_final_v2_variants_B_stat_mask_1KGP_archaic_gnomAD_introgr_hapR2_introgressed_hub, gzfile("../../results/enrichment_GTEx_QTLs_v2/filter_final_v2_variants_B_stat_mask_1KGP_archaic_gnomAD_introgr_hapR2_introgressed_hub.txt.gz"))

# save adaptively introgressed all
filter_final_v2_variants_B_stat_mask_1KGP_archaic_gnomAD_introgr_hapR2_adaptive_hub <- filter_final_v2_variants_B_stat_mask_1KGP_archaic_gnomAD_introgr_hapR2_hub %>% 
	filter(hub_in_gittelman_2016_EUR, hub_in_final_study_adaptive)
write_tsv(filter_final_v2_variants_B_stat_mask_1KGP_archaic_gnomAD_introgr_hapR2_adaptive_hub, gzfile("../../results/enrichment_GTEx_QTLs_v2/filter_final_v2_variants_B_stat_mask_1KGP_archaic_gnomAD_introgr_hapR2_adaptive_hub.txt.gz"))
