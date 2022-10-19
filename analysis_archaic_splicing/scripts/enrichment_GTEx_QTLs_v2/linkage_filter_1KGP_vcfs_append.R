#!/bin/R

# Add EUR_AC_bin, EUR_hapR2tag_bin

# load packages
library(tidyverse)
library(data.table)

# load EUR_AC
ALL_1KGP_phase3_hub <- as_tibble(fread("../../results/preprocess_1KGP_SNPs/ALL_1KGP_phase3_hub.txt.gz"))
ALL_1KGP_phase3_hub_temp <- ALL_1KGP_phase3_hub %>% 
	mutate(temp_key = paste(hub_variant_CHROM, hub_variant_POS, sep="_")) %>% 
	dplyr::select(hub_variant_ID, temp_key, EUR_AC)

# create EUR_AC_bin
ALL_1KGP_phase3_hub_temp$EUR_AC_bin <- 
	as.numeric(cut(ALL_1KGP_phase3_hub_temp$EUR_AC, c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 
	quantile(filter(ALL_1KGP_phase3_hub_temp, EUR_AC>=10)$EUR_AC, c(1:40)/40))))

# load hapR2 info
linkage_hapR2_SNP_info <- as_tibble(fread("../../results/linkage_hapR2_SNPs/ALL.chr.phase3_shapeit2_mvncall_integrated_v5a.20130502.hapR2.collated.gz"))
linkage_hapR2_SNP_info <- linkage_hapR2_SNP_info %>% 
	mutate(temp_key = paste(CHR, POS1, sep="_")) %>% 
	mutate(EUR_hapR2tag = n) %>% dplyr::select(temp_key, EUR_hapR2tag)
linkage_hapR2_SNP_info$EUR_hapR2tag_bin <- 
	as.numeric(cut_number(linkage_hapR2_SNP_info$EUR_hapR2tag, 20))

# add EUR_hapR2tag
final_v2_variants_B_stat_mask_1KGP_archaic_gnomAD_hub <- as_tibble(fread("../../results/preprocess_1KGP_SNPs/final_v2_variants_B_stat_mask_1KGP_archaic_gnomAD_hub.txt.gz"))
final_v2_variants_B_stat_mask_1KGP_archaic_gnomAD_hapR2_hub <- final_v2_variants_B_stat_mask_1KGP_archaic_gnomAD_hub %>% 
	# dplyr::select(-starts_with("B_statistic")) %>% 
	left_join(ALL_1KGP_phase3_hub_temp) %>% left_join(linkage_hapR2_SNP_info) %>% dplyr::select(-temp_key)
write_tsv(final_v2_variants_B_stat_mask_1KGP_archaic_gnomAD_hapR2_hub, gzfile("../../results/preprocess_1KGP_SNPs/final_v2_variants_B_stat_mask_1KGP_archaic_gnomAD_hapR2_hub.txt.gz"))

ALL_1KGP_phase3_greater0.9_hub <- as_tibble(fread("../../results/preprocess_1KGP_SNPs/ALL_1KGP_phase3_greater0.9_hub.txt.gz"))
ALL_1KGP_phase3_greater0.9_hapR2_hub <- ALL_1KGP_phase3_greater0.9_hub %>% 
	# dplyr::select(-starts_with("B_statistic")) %>% 
	left_join(ALL_1KGP_phase3_hub_temp) %>% left_join(linkage_hapR2_SNP_info) %>% dplyr::select(-temp_key)
write_tsv(ALL_1KGP_phase3_greater0.9_hapR2_hub, gzfile("../../results/preprocess_1KGP_SNPs/ALL_1KGP_phase3_greater0.9_hapR2_hub.txt.gz"))

ALL_1KGP_phase3_MAF0.01_hub <- as_tibble(fread("../../results/preprocess_1KGP_SNPs/ALL_1KGP_phase3_MAF0.01_hub.txt.gz"))
ALL_1KGP_phase3_MAF0.01_hapR2_hub <- ALL_1KGP_phase3_MAF0.01_hub %>% 
	# dplyr::select(-starts_with("B_statistic")) %>% 
	left_join(ALL_1KGP_phase3_hub_temp) %>% left_join(linkage_hapR2_SNP_info) %>% dplyr::select(-temp_key)
write_tsv(ALL_1KGP_phase3_MAF0.01_hapR2_hub, gzfile("../../results/preprocess_1KGP_SNPs/ALL_1KGP_phase3_MAF0.01_hapR2_hub.txt.gz"))

archaics_all_1KGP_archaic_gnomAD_hub <- as_tibble(fread("../../results/preprocess_1KGP_SNPs/archaics_all_1KGP_archaic_gnomAD_hub.txt.gz"))
archaics_all_1KGP_archaic_gnomAD_hapR2_hub <- archaics_all_1KGP_archaic_gnomAD_hub %>% 
	# dplyr::select(-starts_with("B_statistic")) %>% 
	left_join(ALL_1KGP_phase3_hub_temp) %>% left_join(linkage_hapR2_SNP_info) %>% dplyr::select(-temp_key)
write_tsv(archaics_all_1KGP_archaic_gnomAD_hapR2_hub, gzfile("../../results/preprocess_1KGP_SNPs/archaics_all_1KGP_archaic_gnomAD_hapR2_hub.txt.gz"))

hg19_REFisDER_1KGP_archaic_gnomAD_hub <- as_tibble(fread("../../results/preprocess_1KGP_SNPs/hg19_REFisDER_1KGP_archaic_gnomAD_hub.txt.gz"))
hg19_REFisDER_1KGP_archaic_gnomAD_hapR2_hub <- hg19_REFisDER_1KGP_archaic_gnomAD_hub %>% 
	# dplyr::select(-starts_with("B_statistic")) %>% 
	left_join(ALL_1KGP_phase3_hub_temp) %>% left_join(linkage_hapR2_SNP_info) %>% dplyr::select(-temp_key)
write_tsv(hg19_REFisDER_1KGP_archaic_gnomAD_hapR2_hub, gzfile("../../results/preprocess_1KGP_SNPs/hg19_REFisDER_1KGP_archaic_gnomAD_hapR2_hub.txt.gz"))

ALL_1KGP_phase3_hub <- as_tibble(fread("../../results/preprocess_1KGP_SNPs/ALL_1KGP_phase3_hub.txt.gz"))
ALL_1KGP_phase3_hapR2_hub <- ALL_1KGP_phase3_hub %>% 
	# dplyr::select(-starts_with("B_statistic")) %>% 
	left_join(ALL_1KGP_phase3_hub_temp) %>% left_join(linkage_hapR2_SNP_info) %>% dplyr::select(-temp_key)
write_tsv(ALL_1KGP_phase3_hapR2_hub, gzfile("../../results/preprocess_1KGP_SNPs/ALL_1KGP_phase3_hapR2_hub.txt.gz"))
