#!/bin/R

library(tidyverse)
library(data.table)

final_v2_variants_B_stat_mask_1KGP_archaic_gnomAD_hub <- as_tibble(fread("../../results/preprocess_1KGP_SNPs/final_v2_variants_B_stat_mask_1KGP_archaic_gnomAD_hub.txt.gz"))

# Partition based on AFR allele frequency
final_v2_variants_B_stat_mask_1KGP_archaic_gnomAD_introgr_hub <- final_v2_variants_B_stat_mask_1KGP_archaic_gnomAD_hub %>% 
	mutate(hub_variant_introgr_filter = case_when(  # identify introgressed allele based on ALT or REF frequency in AFR
		!(hub_in_ldvernot_akey_2016|hub_in_browning_2018) ~ "INTROGR_FAIL",  # not introgressed
		# is.na(AFR_AF)&(is.na(gnomAD_AF)) ~ "AFR_ALT_NA",  # if absent in 1KGP, and also absent in gnomAD, assume AFR AF is 0
		# is.na(AFR_AF)&(gnomAD_AF <= 0.05) ~ "AFR_ALT_NA",  # if absent in 1KGP, and gnomAD AF is low, assume AFR AF is low
		# is.na(AFR_AF) ~ "AFR_NA_FAIL",  # otherwise, exclude (unresolved absent in 1KGP)
		is.na(AFR_AF) ~ "AFR_ALT",  # if absent in 1KGP, assume AFR AF is 0, very few exceptions
		AFR_AF <= 0.01 ~ "AFR_ALT",  # low in 1KGP
		AFR_AF >= 0.99 ~ "AFR_REF",  # high in 1KGP
		TRUE ~ "AFR_AF_FAIL"  # otherwise, exclude (intermediate AF)
	))

sort(table(final_v2_variants_B_stat_mask_1KGP_archaic_gnomAD_introgr_hub$hub_variant_introgr_filter))

# Must match one of the archaic alleles
final_v2_variants_B_stat_mask_1KGP_archaic_gnomAD_introgr_hub <- final_v2_variants_B_stat_mask_1KGP_archaic_gnomAD_introgr_hub %>% 
	mutate(hub_variant_introgr_filter = 
		ifelse(hub_in_final_study_introgressed&(!archaic_mask_union), 
			"ARC_MATCH_FAIL", hub_variant_introgr_filter)) %>% 
	mutate(hub_variant_introgr_filter = 
		ifelse((hub_variant_introgr_filter == "AFR_ALT"), 
			ifelse(
				((altai_denisovan_AC %in% c(1, 2)) | (altai_neanderthal_AC %in% c(1, 2)) | (vindija_neanderthal_AC %in% c(1, 2)) | (chagyrskaya_neanderthal_AC %in% c(1, 2))), 
				hub_variant_introgr_filter, "ARC_MATCH_FAIL"),
			hub_variant_introgr_filter)) %>% 
	mutate(hub_variant_introgr_filter = 
		ifelse((hub_variant_introgr_filter == "AFR_REF"), 
			ifelse(
				((altai_denisovan_AC %in% c(1, NA)) | (altai_neanderthal_AC %in% c(1, NA)) | (vindija_neanderthal_AC %in% c(1, NA)) | (chagyrskaya_neanderthal_AC %in% c(1, NA))), 
				hub_variant_introgr_filter, "ARC_MATCH_FAIL"),
			hub_variant_introgr_filter))
sort(table(final_v2_variants_B_stat_mask_1KGP_archaic_gnomAD_introgr_hub$hub_variant_introgr_filter))

# Get introgressed allele
final_v2_variants_B_stat_mask_1KGP_archaic_gnomAD_introgr_hub <- final_v2_variants_B_stat_mask_1KGP_archaic_gnomAD_introgr_hub %>% 
	mutate(hub_variant_introgr_allele = case_when(
		hub_variant_introgr_filter == "AFR_ALT" ~ hub_variant_ALT, 
		hub_variant_introgr_filter == "AFR_REF" ~ hub_variant_REF, 
		TRUE ~ "NA"
	)) %>% 
	mutate(hub_variant_introgr_allele = ifelse(hub_variant_introgr_allele == "NA", NA, hub_variant_introgr_allele))

# Get introgressed derived vs ancestral
final_v2_variants_B_stat_mask_1KGP_archaic_gnomAD_introgr_hub <- final_v2_variants_B_stat_mask_1KGP_archaic_gnomAD_introgr_hub %>% 
	mutate(hub_variant_introgr_class = case_when(
		hub_variant_introgr_allele == hub_variant_DER ~ "introgr_DER",
		hub_variant_introgr_allele == hub_variant_ANC ~ "introgr_ANC",
		hub_variant_introgr_allele %in% c("A", "C", "T", "G") ~ "introgr_ANC_FAIL", 
		hub_variant_introgr_filter == "AFR_AF_FAIL" ~ "introgr_AF_FAIL", 
		hub_variant_introgr_filter == "ARC_MATCH_FAIL" ~ "introgr_MATCH_FAIL",
		TRUE ~ "NA"
	)) %>% 
	mutate(hub_variant_introgr_class = ifelse(hub_variant_introgr_class == "NA", NA, hub_variant_introgr_class))

final_v2_variants_B_stat_mask_1KGP_archaic_gnomAD_introgr_hub <- final_v2_variants_B_stat_mask_1KGP_archaic_gnomAD_introgr_hub %>% 
	mutate(hub_in_final_study_introgressed_DER = ((hub_variant_introgr_class=="introgr_DER")&(!is.na(hub_variant_introgr_class)))) %>% 
	mutate(hub_in_final_study_introgressed_ANC = ((hub_variant_introgr_class=="introgr_ANC")&(!is.na(hub_variant_introgr_class)))) %>% 
	mutate(hub_in_final_study_introgressed_AF_FAIL = ((hub_variant_introgr_class=="introgr_AF_FAIL")&(!is.na(hub_variant_introgr_class)))) %>% 
	mutate(hub_in_final_study_introgressed_MATCH_FAIL = ((hub_variant_introgr_class=="introgr_MATCH_FAIL")&(!is.na(hub_variant_introgr_class))))

final_v2_variants_B_stat_mask_1KGP_archaic_gnomAD_introgr_hub <- final_v2_variants_B_stat_mask_1KGP_archaic_gnomAD_introgr_hub %>% 
	mutate(hub_in_final_study_adaptive_DER = (hub_in_final_study_adaptive & hub_in_final_study_introgressed_DER)) %>% 
	mutate(hub_in_final_study_adaptive_ANC = (hub_in_final_study_adaptive & hub_in_final_study_introgressed_ANC)) %>% 
	mutate(hub_in_final_study_adaptive_AF_FAIL = (hub_in_final_study_adaptive & hub_in_final_study_introgressed_AF_FAIL)) %>% 
	mutate(hub_in_final_study_adaptive_MATCH_FAIL = (hub_in_final_study_adaptive & hub_in_final_study_introgressed_MATCH_FAIL))

# check
sort(table(final_v2_variants_B_stat_mask_1KGP_archaic_gnomAD_introgr_hub$hub_variant_introgr_filter))
sort(table(final_v2_variants_B_stat_mask_1KGP_archaic_gnomAD_introgr_hub$hub_variant_introgr_class))

# save
tab_hub_variant_introgr_filter <- as_tibble(data.frame(sort(table(final_v2_variants_B_stat_mask_1KGP_archaic_gnomAD_introgr_hub$hub_variant_introgr_filter))))
tab_hub_variant_introgr_class <- as_tibble(data.frame(sort(table(final_v2_variants_B_stat_mask_1KGP_archaic_gnomAD_introgr_hub$hub_variant_introgr_class))))
names(tab_hub_variant_introgr_filter) <- c("hub_variant_introgr_filter", "n")
names(tab_hub_variant_introgr_class) <- c("hub_variant_introgr_class", "n")
write_tsv(tab_hub_variant_introgr_filter, "../../results/preprocess_1KGP_SNPs/tab_hub_variant_introgr_filter.txt")
write_tsv(tab_hub_variant_introgr_class, "../../results/preprocess_1KGP_SNPs/tab_hub_variant_introgr_class.txt")

# # Get neanderthal and denisovan desert variants
final_v2_variants_B_stat_mask_1KGP_archaic_gnomAD_introgr_hub <- final_v2_variants_B_stat_mask_1KGP_archaic_gnomAD_introgr_hub %>% 
	mutate(hub_in_final_study_archaic_introgressed = (hub_in_final_study_archaic & (hub_in_final_study_introgressed))) %>% 
	mutate(hub_in_final_study_nean_introgressed = (hub_in_final_study_nean & (hub_in_final_study_introgressed))) %>% 
	mutate(hub_in_final_study_deni_introgressed = (hub_in_final_study_deni & (hub_in_final_study_introgressed))) %>% 
	mutate(hub_in_final_study_archaic_nonintrogressed = (hub_in_final_study_archaic & !(hub_in_final_study_introgressed))) %>% 
	mutate(hub_in_final_study_nean_nonintrogressed = (hub_in_final_study_nean & !(hub_in_final_study_introgressed))) %>% 
	mutate(hub_in_final_study_deni_nonintrogressed = (hub_in_final_study_deni & !(hub_in_final_study_introgressed)))

# Save
write_tsv(final_v2_variants_B_stat_mask_1KGP_archaic_gnomAD_introgr_hub, gzfile("../../results/preprocess_1KGP_SNPs/final_v2_variants_B_stat_mask_1KGP_archaic_gnomAD_introgr_hub.txt.gz"))

# # 	test
# a <- final_v2_variants_B_stat_mask_1KGP_archaic_gnomAD_introgr_hub %>% filter(hub_in_final_study_adaptive_DER) %>% arrange(AFR_AF)
# write_tsv(a, "test_afa.txt")
# b <- final_v2_variants_B_stat_mask_1KGP_archaic_gnomAD_introgr_hub %>% filter(hub_in_final_study_adaptive_ANC) %>% arrange(AFR_AF)
# write_tsv(b, "test_afb.txt")
