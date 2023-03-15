#!/bin/sh

# Reprocess to get updated lineage-specific variants and combine with old library data

library(tidyverse)
library(data.table)

# old version
merge_variants_B_stat_mask_1KGP_archaic_gnomAD_hub <- as_tibble(fread("../../results/preprocess_1KGP_SNPs/merge_variants_B_stat_mask_1KGP_archaic_gnomAD_hub.txt.gz"))

check_within <- function(query_tb, subject_tb) {
	return(table(query_tb$hub_variant_ID %in% subject_tb$hub_variant_ID))
}

# load helpers
ALL_1KGP_phase3_greater0.9_gnomAD_hub <- as_tibble(fread("../../results/preprocess_1KGP_SNPs/ALL_1KGP_phase3_greater0.9_gnomAD_hub.txt.gz"))
ALL_1KGP_phase3_greater0.5_gnomAD_hub <- as_tibble(fread("../../results/preprocess_1KGP_SNPs/ALL_1KGP_phase3_greater0.5_gnomAD_hub.txt.gz"))
archaics_all_1KGP_archaic_gnomAD_hub <- as_tibble(fread("../../results/preprocess_1KGP_SNPs/archaics_all_1KGP_archaic_gnomAD_hub.txt.gz"))
hg19_REFisDER_1KGP_archaic_gnomAD_hub <- as_tibble(fread("../../results/preprocess_1KGP_SNPs/hg19_REFisDER_1KGP_archaic_gnomAD_hub.txt.gz"))

# modern nearly:
# 	REFisANC, modern > 0.9, !(archaic >= 1), archaic_mask_inter
modern_nearly_REFisANC <- 
 	# add annotations, is a nearly_0.99 subset, so works for all SNPs, in particular gnomAD annotations
	ALL_1KGP_phase3_greater0.9_gnomAD_hub %>% 
	left_join(dplyr::select(merge_variants_B_stat_mask_1KGP_archaic_gnomAD_hub, !starts_with("hub_in"))) %>% 
	filter(hub_variant_REF == hub_variant_ANC) %>% 
	filter(archaic_mask_inter) %>% 
	filter((AF > 0.9) & (gnomAD_AF > 0.9)) %>% 
	filter(!(hub_variant_ID %in% filter(archaics_all_1KGP_archaic_gnomAD_hub, archaic_AC >= 1)$hub_variant_ID))
check_within(modern_nearly_REFisANC, merge_variants_B_stat_mask_1KGP_archaic_gnomAD_hub)

# 	REFisDER, modern < 0.1, archaic >= 8, archaic_mask_inter
modern_nearly_REFisDER <- 
	hg19_REFisDER_1KGP_archaic_gnomAD_hub %>% 
	filter(archaic_AC >= 8) %>% 
	filter(archaic_mask_inter) %>% 
	filter(((AF < 0.1) | (is.na(AF))) & ((gnomAD_AF < 0.1) | (is.na(gnomAD_AF))))
check_within(modern_nearly_REFisDER, merge_variants_B_stat_mask_1KGP_archaic_gnomAD_hub)
# a variant missing, was previously eliminated due to having high gnomAD frequency, now not eliminated due to not passing gnomAD filter


# modern nearly_0.5:
# 	REFisANC, modern > 0.5, !(archaic >= 1), archaic_mask_inter
modern_nearly_0.5_REFisANC <- 
 	# add annotations, is a nearly_0.5 subset, so works for all SNPs, in particular gnomAD annotations
	ALL_1KGP_phase3_greater0.5_gnomAD_hub %>% 
	left_join(dplyr::select(merge_variants_B_stat_mask_1KGP_archaic_gnomAD_hub, !starts_with("hub_in"))) %>% 
	filter(hub_variant_REF == hub_variant_ANC) %>% 
	filter(archaic_mask_inter) %>% 
	filter((AF > 0.5) & (gnomAD_AF > 0.5)) %>% 
	filter(!(hub_variant_ID %in% filter(archaics_all_1KGP_archaic_gnomAD_hub, archaic_AC >= 1)$hub_variant_ID))
check_within(modern_nearly_0.5_REFisANC, merge_variants_B_stat_mask_1KGP_archaic_gnomAD_hub)

# 	REFisDER, modern < 0.5, archaic >= 8, archaic_mask_inter
modern_nearly_0.5_REFisDER <- 
	hg19_REFisDER_1KGP_archaic_gnomAD_hub %>% 
	filter(archaic_AC >= 8) %>% 
	filter(archaic_mask_inter) %>% 
	filter(((AF < 0.5) | (is.na(AF))) & ((gnomAD_AF < 0.5) | (is.na(gnomAD_AF))))
check_within(modern_nearly_0.5_REFisDER, merge_variants_B_stat_mask_1KGP_archaic_gnomAD_hub)


# modern nearly_0.6:
# 	REFisANC, modern > 0.6, !(archaic >= 1), archaic_mask_inter
modern_nearly_0.6_REFisANC <- 
 	# add annotations, is a nearly_0.6 subset, so works for all SNPs, in particular gnomAD annotations
	ALL_1KGP_phase3_greater0.5_gnomAD_hub %>% 
	left_join(dplyr::select(merge_variants_B_stat_mask_1KGP_archaic_gnomAD_hub, !starts_with("hub_in"))) %>% 
	filter(hub_variant_REF == hub_variant_ANC) %>% 
	filter(archaic_mask_inter) %>% 
	filter((AF > 0.6) & (gnomAD_AF > 0.6)) %>% 
	filter(!(hub_variant_ID %in% filter(archaics_all_1KGP_archaic_gnomAD_hub, archaic_AC >= 1)$hub_variant_ID))
check_within(modern_nearly_0.6_REFisANC, merge_variants_B_stat_mask_1KGP_archaic_gnomAD_hub)

# 	REFisDER, modern < 0.6, archaic >= 8, archaic_mask_inter
modern_nearly_0.6_REFisDER <- 
	hg19_REFisDER_1KGP_archaic_gnomAD_hub %>% 
	filter(archaic_AC >= 8) %>% 
	filter(archaic_mask_inter) %>% 
	filter(((AF < 0.4) | (is.na(AF))) & ((gnomAD_AF < 0.4) | (is.na(gnomAD_AF))))
check_within(modern_nearly_0.6_REFisDER, merge_variants_B_stat_mask_1KGP_archaic_gnomAD_hub)


# modern nearly_0.7:
# 	REFisANC, modern > 0.7, !(archaic >= 1), archaic_mask_inter
modern_nearly_0.7_REFisANC <- 
 	# add annotations, is a nearly_0.7 subset, so works for all SNPs, in particular gnomAD annotations
	ALL_1KGP_phase3_greater0.5_gnomAD_hub %>% 
	left_join(dplyr::select(merge_variants_B_stat_mask_1KGP_archaic_gnomAD_hub, !starts_with("hub_in"))) %>% 
	filter(hub_variant_REF == hub_variant_ANC) %>% 
	filter(archaic_mask_inter) %>% 
	filter((AF > 0.7) & (gnomAD_AF > 0.7)) %>% 
	filter(!(hub_variant_ID %in% filter(archaics_all_1KGP_archaic_gnomAD_hub, archaic_AC >= 1)$hub_variant_ID))
check_within(modern_nearly_0.7_REFisANC, merge_variants_B_stat_mask_1KGP_archaic_gnomAD_hub)

# 	REFisDER, modern < 0.7, archaic >= 8, archaic_mask_inter
modern_nearly_0.7_REFisDER <- 
	hg19_REFisDER_1KGP_archaic_gnomAD_hub %>% 
	filter(archaic_AC >= 8) %>% 
	filter(archaic_mask_inter) %>% 
	filter(((AF < 0.3) | (is.na(AF))) & ((gnomAD_AF < 0.3) | (is.na(gnomAD_AF))))
check_within(modern_nearly_0.7_REFisDER, merge_variants_B_stat_mask_1KGP_archaic_gnomAD_hub)


# modern nearly_0.8:
# 	REFisANC, modern > 0.8, !(archaic >= 1), archaic_mask_inter
modern_nearly_0.8_REFisANC <- 
 	# add annotations, is a nearly_0.8 subset, so works for all SNPs, in particular gnomAD annotations
	ALL_1KGP_phase3_greater0.5_gnomAD_hub %>% 
	left_join(dplyr::select(merge_variants_B_stat_mask_1KGP_archaic_gnomAD_hub, !starts_with("hub_in"))) %>% 
	filter(hub_variant_REF == hub_variant_ANC) %>% 
	filter(archaic_mask_inter) %>% 
	filter((AF > 0.8) & (gnomAD_AF > 0.8)) %>% 
	filter(!(hub_variant_ID %in% filter(archaics_all_1KGP_archaic_gnomAD_hub, archaic_AC >= 1)$hub_variant_ID))
check_within(modern_nearly_0.8_REFisANC, merge_variants_B_stat_mask_1KGP_archaic_gnomAD_hub)

# 	REFisDER, modern < 0.8, archaic >= 8, archaic_mask_inter
modern_nearly_0.8_REFisDER <- 
	hg19_REFisDER_1KGP_archaic_gnomAD_hub %>% 
	filter(archaic_AC >= 8) %>% 
	filter(archaic_mask_inter) %>% 
	filter(((AF < 0.2) | (is.na(AF))) & ((gnomAD_AF < 0.2) | (is.na(gnomAD_AF))))
check_within(modern_nearly_0.8_REFisDER, merge_variants_B_stat_mask_1KGP_archaic_gnomAD_hub)


# modern nearly_0.9:
# 	REFisANC, modern > 0.9, !(archaic >= 1), archaic_mask_inter
modern_nearly_0.9_REFisANC <- 
 	# add annotations, is a nearly_0.9 subset, so works for all SNPs, in particular gnomAD annotations
	ALL_1KGP_phase3_greater0.9_gnomAD_hub %>% 
	left_join(dplyr::select(merge_variants_B_stat_mask_1KGP_archaic_gnomAD_hub, !starts_with("hub_in"))) %>% 
	filter(hub_variant_REF == hub_variant_ANC) %>% 
	filter(archaic_mask_inter) %>% 
	filter((AF > 0.9) & (gnomAD_AF > 0.9)) %>% 
	filter(!(hub_variant_ID %in% filter(archaics_all_1KGP_archaic_gnomAD_hub, archaic_AC >= 1)$hub_variant_ID))
check_within(modern_nearly_0.9_REFisANC, merge_variants_B_stat_mask_1KGP_archaic_gnomAD_hub)

# 	REFisDER, modern < 0.9, archaic >= 8, archaic_mask_inter
modern_nearly_0.9_REFisDER <- 
	hg19_REFisDER_1KGP_archaic_gnomAD_hub %>% 
	filter(archaic_AC >= 8) %>% 
	filter(archaic_mask_inter) %>% 
	filter(((AF < 0.1) | (is.na(AF))) & ((gnomAD_AF < 0.1) | (is.na(gnomAD_AF))))
check_within(modern_nearly_0.9_REFisDER, merge_variants_B_stat_mask_1KGP_archaic_gnomAD_hub)


# modern nearly_0.99:
# 	REFisANC, modern > 0.99, !(archaic >= 1), archaic_mask_inter
modern_nearly_0.99_REFisANC <- 
 	# add annotations, is a nearly_0.99 subset, so works for all SNPs, in particular gnomAD annotations
	ALL_1KGP_phase3_greater0.9_gnomAD_hub %>% 
	left_join(dplyr::select(merge_variants_B_stat_mask_1KGP_archaic_gnomAD_hub, !starts_with("hub_in"))) %>% 
	filter(hub_variant_REF == hub_variant_ANC) %>% 
	filter(archaic_mask_inter) %>% 
	filter((AF > 0.99) & (gnomAD_AF > 0.99)) %>% 
	filter(!(hub_variant_ID %in% filter(archaics_all_1KGP_archaic_gnomAD_hub, archaic_AC >= 1)$hub_variant_ID))
check_within(modern_nearly_0.99_REFisANC, merge_variants_B_stat_mask_1KGP_archaic_gnomAD_hub)

# 	REFisDER, modern < 0.99, archaic >= 8, archaic_mask_inter
modern_nearly_0.99_REFisDER <- 
	hg19_REFisDER_1KGP_archaic_gnomAD_hub %>% 
	filter(archaic_AC >= 8) %>% 
	filter(archaic_mask_inter) %>% 
	filter(((AF < 0.01) | (is.na(AF))) & ((gnomAD_AF < 0.01) | (is.na(gnomAD_AF))))
check_within(modern_nearly_0.99_REFisDER, merge_variants_B_stat_mask_1KGP_archaic_gnomAD_hub)


# modern nearly_0.999:
# 	REFisANC, modern > 0.999, !(archaic >= 1), archaic_mask_inter
modern_nearly_0.999_REFisANC <- 
 	# add annotations, is a nearly_0.999 subset, so works for all SNPs, in particular gnomAD annotations
	ALL_1KGP_phase3_greater0.9_gnomAD_hub %>% 
	left_join(dplyr::select(merge_variants_B_stat_mask_1KGP_archaic_gnomAD_hub, !starts_with("hub_in"))) %>% 
	filter(hub_variant_REF == hub_variant_ANC) %>% 
	filter(archaic_mask_inter) %>% 
	filter((AF > 0.999) & (gnomAD_AF > 0.999)) %>% 
	filter(!(hub_variant_ID %in% filter(archaics_all_1KGP_archaic_gnomAD_hub, archaic_AC >= 1)$hub_variant_ID))
check_within(modern_nearly_0.999_REFisANC, merge_variants_B_stat_mask_1KGP_archaic_gnomAD_hub)

# 	REFisDER, modern < 0.999, archaic >= 8, archaic_mask_inter
modern_nearly_0.999_REFisDER <- 
	hg19_REFisDER_1KGP_archaic_gnomAD_hub %>% 
	filter(archaic_AC >= 8) %>% 
	filter(archaic_mask_inter) %>% 
	filter(((AF < 0.001) | (is.na(AF))) & ((gnomAD_AF < 0.001) | (is.na(gnomAD_AF))))
check_within(modern_nearly_0.999_REFisDER, merge_variants_B_stat_mask_1KGP_archaic_gnomAD_hub)


# modern fixed:
# 	REFisANC, modern = 1, !(archaic >= 1), archaic_mask_inter
modern_fixed_REFisANC <- 
 	# add annotations, is a fixed subset, so works for all SNPs, in particular gnomAD annotations
	ALL_1KGP_phase3_greater0.9_gnomAD_hub %>% 
	left_join(dplyr::select(merge_variants_B_stat_mask_1KGP_archaic_gnomAD_hub, !starts_with("hub_in"))) %>% 
	filter(hub_variant_REF == hub_variant_ANC) %>% 
	filter(archaic_mask_inter) %>% 
	filter((AF == 1) & (gnomAD_AF == 1)) %>% 
	filter(!(hub_variant_ID %in% filter(archaics_all_1KGP_archaic_gnomAD_hub, archaic_AC >= 1)$hub_variant_ID))
check_within(modern_fixed_REFisANC, merge_variants_B_stat_mask_1KGP_archaic_gnomAD_hub)

# 	REFisDER, modern = 0, archaic >= 8, archaic_mask_inter
modern_fixed_REFisDER <- 
	hg19_REFisDER_1KGP_archaic_gnomAD_hub %>% 
	filter(archaic_AC >= 8) %>% 
	filter(archaic_mask_inter) %>% 
	filter(((AF == 0) | (is.na(AF))) & ((gnomAD_AF == 0) | (is.na(gnomAD_AF))))
check_within(modern_fixed_REFisDER, merge_variants_B_stat_mask_1KGP_archaic_gnomAD_hub)


# archaic nearly:
# 	REFisANC, archaic >= 7, modern < 0.01, archaic_mask_inter
archaic_nearly_REFisANC <- 
	archaics_all_1KGP_archaic_gnomAD_hub %>% 
	filter(hub_variant_REF == hub_variant_ANC) %>% 
	filter(archaic_mask_inter) %>% 
	filter(archaic_AC >= 7) %>% 
	filter(((AFR_AF < 0.01) | (is.na(AFR_AF))) & ((gnomAD_AF_afr < 0.01) | (is.na(gnomAD_AF_afr))))
check_within(archaic_nearly_REFisANC, merge_variants_B_stat_mask_1KGP_archaic_gnomAD_hub)

# 	REFisDER, modern > 0.99, archaic <= 1, archaic_mask_inter
archaic_nearly_REFisDER <- 
	hg19_REFisDER_1KGP_archaic_gnomAD_hub %>% 
	filter(((AFR_AF > 0.99) & (gnomAD_AF_afr > 0.99)) | ((is.na(AFR_AF) & (gnomAD_AF_afr > 0.99))) | ((AFR_AF > 0.99) & (is.na(gnomAD_AF_afr)))) %>% 
	filter((archaic_AC <= 1) | (is.na(archaic_AC))) %>% 
	filter(archaic_mask_inter)
check_within(archaic_nearly_REFisDER, merge_variants_B_stat_mask_1KGP_archaic_gnomAD_hub)
# small number of new variants

# note 02-06-2023:
# issue with 277 archaic specific variants
# mistakenly retained during earlier revision
# these have AF > 0.9 but AFR_AF < 0.99 
# or gnomAD_AFR_AF < 0.99

# archaic fixed:
# 	REFisANC, archaic >= 8, modern = 0, archaic_mask_inter
archaic_fixed_REFisANC <- 
	archaics_all_1KGP_archaic_gnomAD_hub %>% 
	filter(hub_variant_REF == hub_variant_ANC) %>% 
	filter(archaic_mask_inter) %>% 
	filter(archaic_AC >= 8) %>% 
	filter(((AFR_AF == 0) | (is.na(AFR_AF))) & ((gnomAD_AF_afr == 0) | (is.na(gnomAD_AF_afr))))
check_within(archaic_fixed_REFisANC, merge_variants_B_stat_mask_1KGP_archaic_gnomAD_hub)

# 	REFisDER, modern = 1, archaic <= 0, archaic_mask_inter
archaic_fixed_REFisDER <- 
	hg19_REFisDER_1KGP_archaic_gnomAD_hub %>% 
	filter(((AFR_AF == 1) & (gnomAD_AF_afr == 1)) | ((is.na(AFR_AF) & (gnomAD_AF_afr == 1))) | ((AFR_AF > 0.99) & (is.na(gnomAD_AF_afr)))) %>% 
	filter((archaic_AC <= 0) | (is.na(archaic_AC))) %>% 
	filter(archaic_mask_inter)
check_within(archaic_fixed_REFisDER, merge_variants_B_stat_mask_1KGP_archaic_gnomAD_hub)
# small number of new variants


# neanderthal nearly:
# 	REFisANC, neanderthal >= 5, modern < 0.1, archaic_mask_inter
neanderthal_nearly_REFisANC <- 
	archaics_all_1KGP_archaic_gnomAD_hub %>% 
	filter(hub_variant_REF == hub_variant_ANC) %>% 
	filter(archaic_mask_inter) %>% 
	filter(neanderthal_AC >= 5, denisovan_AC == 0) %>% 
	filter(((AFR_AF < 0.01) | (is.na(AFR_AF))) & ((gnomAD_AF_afr < 0.01) | (is.na(gnomAD_AF_afr))))
check_within(neanderthal_nearly_REFisANC, merge_variants_B_stat_mask_1KGP_archaic_gnomAD_hub)

# 	REFisDER, modern > 0.9, neanderthal <= 1, denisovan == 2, archaic_mask_inter
neanderthal_nearly_REFisDER <- 
	hg19_REFisDER_1KGP_archaic_gnomAD_hub %>% 
	filter(((AFR_AF > 0.99) & (gnomAD_AF_afr > 0.99)) | ((is.na(AFR_AF) & (gnomAD_AF_afr > 0.99))) | ((AFR_AF > 0.99) & (is.na(gnomAD_AF_afr)))) %>% 
	filter(neanderthal_AC <= 1, denisovan_AC == 2) %>% 
	filter(archaic_mask_inter)
check_within(neanderthal_nearly_REFisDER, merge_variants_B_stat_mask_1KGP_archaic_gnomAD_hub)
# tiny number of new variants 


# neanderthal fixed:
# 	REFisANC, neanderthal >= 6, modern = 0, archaic_mask_inter
neanderthal_fixed_REFisANC <- 
	archaics_all_1KGP_archaic_gnomAD_hub %>% 
	filter(hub_variant_REF == hub_variant_ANC) %>% 
	filter(archaic_mask_inter) %>% 
	filter(neanderthal_AC >= 6, denisovan_AC == 0) %>% 
	filter(((AFR_AF == 0) | (is.na(AFR_AF))) & ((gnomAD_AF_afr == 0) | (is.na(gnomAD_AF_afr))))
check_within(neanderthal_fixed_REFisANC, merge_variants_B_stat_mask_1KGP_archaic_gnomAD_hub)

# 	REFisDER, modern = 1, neanderthal <= 0, denisovan == 2, archaic_mask_inter
neanderthal_fixed_REFisDER <- 
	hg19_REFisDER_1KGP_archaic_gnomAD_hub %>% 
	filter(((AFR_AF == 1) & (gnomAD_AF_afr == 1)) | ((is.na(AFR_AF) & (gnomAD_AF_afr == 1))) | ((AFR_AF == 1) & (is.na(gnomAD_AF_afr)))) %>% 
	filter(neanderthal_AC <= 0, denisovan_AC == 2) %>% 
	filter(archaic_mask_inter)
check_within(neanderthal_fixed_REFisDER, merge_variants_B_stat_mask_1KGP_archaic_gnomAD_hub)
# tiny number of new variants 


# denisovan nearly:
# 	REFisANC, denisovan == 2, modern < 0.1, archaic_mask_inter
denisovan_nearly_REFisANC <- 
	archaics_all_1KGP_archaic_gnomAD_hub %>% 
	filter(hub_variant_REF == hub_variant_ANC) %>% 
	filter(archaic_mask_inter) %>% 
	filter(denisovan_AC == 2, neanderthal_AC == 0) %>% 
	filter(((AFR_AF < 0.01) | (is.na(AFR_AF))) & ((gnomAD_AF_afr < 0.01) | (is.na(gnomAD_AF_afr))))
check_within(denisovan_nearly_REFisANC, merge_variants_B_stat_mask_1KGP_archaic_gnomAD_hub)

# 	REFisDER, modern > 0.9, neandertal == 6, denisovan == 0, archaic_mask_inter
denisovan_nearly_REFisDER <- 
	hg19_REFisDER_1KGP_archaic_gnomAD_hub %>% 
	filter(((AFR_AF > 0.99) & (gnomAD_AF_afr > 0.99)) | ((is.na(AFR_AF) & (gnomAD_AF_afr > 0.99))) | ((AFR_AF > 0.99) & (is.na(gnomAD_AF_afr)))) %>% 
	filter(denisovan_AC == 0, neanderthal_AC == 6) %>% 
	filter(archaic_mask_inter)
check_within(denisovan_nearly_REFisDER, merge_variants_B_stat_mask_1KGP_archaic_gnomAD_hub)
# tiny number of new variants


# denisovan fixed:
# 	REFisANC, denisovan == 2, modern = 0, archaic_mask_inter
denisovan_fixed_REFisANC <- 
	archaics_all_1KGP_archaic_gnomAD_hub %>% 
	filter(hub_variant_REF == hub_variant_ANC) %>% 
	filter(archaic_mask_inter) %>% 
	filter(denisovan_AC == 2, neanderthal_AC == 0) %>% 
	filter(((AFR_AF == 0) | (is.na(AFR_AF))) & ((gnomAD_AF_afr == 0) | (is.na(gnomAD_AF_afr))))
check_within(denisovan_fixed_REFisANC, merge_variants_B_stat_mask_1KGP_archaic_gnomAD_hub)

# 	REFisDER, modern = 1, neandertal == 6, denisovan == 0, archaic_mask_inter
denisovan_fixed_REFisDER <- 
	hg19_REFisDER_1KGP_archaic_gnomAD_hub %>% 
	filter(((AFR_AF == 1) & (gnomAD_AF_afr == 1)) | ((is.na(AFR_AF) & (gnomAD_AF_afr == 1))) | ((AFR_AF == 1) & (is.na(gnomAD_AF_afr)))) %>% 
	filter(denisovan_AC == 0, neanderthal_AC == 6) %>% 
	filter(archaic_mask_inter)
check_within(denisovan_fixed_REFisDER, merge_variants_B_stat_mask_1KGP_archaic_gnomAD_hub)
# tiny number of new variants


# clean up and combine
modern_nearly_REFisANC_final <- modern_nearly_REFisANC %>% 
	dplyr::select(names(modern_nearly_REFisANC)) %>% 
	mutate(names_temp = "hub_in_final_study_modern", values_temp = TRUE)
modern_nearly_REFisDER_final <- modern_nearly_REFisDER %>% 
	dplyr::select(names(modern_nearly_REFisANC)) %>% 
	mutate(names_temp = "hub_in_final_study_modern", values_temp = TRUE)
archaic_nearly_REFisANC_final <- archaic_nearly_REFisANC %>% 
	dplyr::select(names(modern_nearly_REFisANC)) %>% 
	mutate(names_temp = "hub_in_final_study_archaic", values_temp = TRUE)
archaic_nearly_REFisDER_final <- archaic_nearly_REFisDER %>% 
	dplyr::select(names(modern_nearly_REFisANC)) %>% 
	mutate(names_temp = "hub_in_final_study_archaic", values_temp = TRUE)
neanderthal_nearly_REFisANC_final <- neanderthal_nearly_REFisANC %>% 
	dplyr::select(names(modern_nearly_REFisANC)) %>% 
	mutate(names_temp = "hub_in_final_study_nean", values_temp = TRUE)
neanderthal_nearly_REFisDER_final <- neanderthal_nearly_REFisDER %>% 
	dplyr::select(names(modern_nearly_REFisANC)) %>% 
	mutate(names_temp = "hub_in_final_study_nean", values_temp = TRUE)
denisovan_nearly_REFisANC_final <- denisovan_nearly_REFisANC %>% 
	dplyr::select(names(modern_nearly_REFisANC)) %>% 
	mutate(names_temp = "hub_in_final_study_deni", values_temp = TRUE)
denisovan_nearly_REFisDER_final <- denisovan_nearly_REFisDER %>% 
	dplyr::select(names(modern_nearly_REFisANC)) %>% 
	mutate(names_temp = "hub_in_final_study_deni", values_temp = TRUE)


modern_nearly_0.5_REFisANC_final <- modern_nearly_0.5_REFisANC %>% 
	dplyr::select(names(modern_nearly_REFisANC)) %>% 
	mutate(names_temp = "hub_in_final_study_modern_nearly_0.5", values_temp = TRUE)
modern_nearly_0.5_REFisDER_final <- modern_nearly_0.5_REFisDER %>% 
	dplyr::select(names(modern_nearly_REFisANC)) %>% 
	mutate(names_temp = "hub_in_final_study_modern_nearly_0.5", values_temp = TRUE)
modern_nearly_0.6_REFisANC_final <- modern_nearly_0.6_REFisANC %>% 
	dplyr::select(names(modern_nearly_REFisANC)) %>% 
	mutate(names_temp = "hub_in_final_study_modern_nearly_0.6", values_temp = TRUE)
modern_nearly_0.6_REFisDER_final <- modern_nearly_0.6_REFisDER %>% 
	dplyr::select(names(modern_nearly_REFisANC)) %>% 
	mutate(names_temp = "hub_in_final_study_modern_nearly_0.6", values_temp = TRUE)
modern_nearly_0.7_REFisANC_final <- modern_nearly_0.7_REFisANC %>% 
	dplyr::select(names(modern_nearly_REFisANC)) %>% 
	mutate(names_temp = "hub_in_final_study_modern_nearly_0.7", values_temp = TRUE)
modern_nearly_0.7_REFisDER_final <- modern_nearly_0.7_REFisDER %>% 
	dplyr::select(names(modern_nearly_REFisANC)) %>% 
	mutate(names_temp = "hub_in_final_study_modern_nearly_0.7", values_temp = TRUE)
modern_nearly_0.8_REFisANC_final <- modern_nearly_0.8_REFisANC %>% 
	dplyr::select(names(modern_nearly_REFisANC)) %>% 
	mutate(names_temp = "hub_in_final_study_modern_nearly_0.8", values_temp = TRUE)
modern_nearly_0.8_REFisDER_final <- modern_nearly_0.8_REFisDER %>% 
	dplyr::select(names(modern_nearly_REFisANC)) %>% 
	mutate(names_temp = "hub_in_final_study_modern_nearly_0.8", values_temp = TRUE)
modern_nearly_0.9_REFisANC_final <- modern_nearly_0.9_REFisANC %>% 
	dplyr::select(names(modern_nearly_REFisANC)) %>% 
	mutate(names_temp = "hub_in_final_study_modern_nearly_0.9", values_temp = TRUE)
modern_nearly_0.9_REFisDER_final <- modern_nearly_0.9_REFisDER %>% 
	dplyr::select(names(modern_nearly_REFisANC)) %>% 
	mutate(names_temp = "hub_in_final_study_modern_nearly_0.9", values_temp = TRUE)

modern_nearly_0.99_REFisANC_final <- modern_nearly_0.99_REFisANC %>% 
	dplyr::select(names(modern_nearly_REFisANC)) %>% 
	mutate(names_temp = "hub_in_final_study_modern_nearly_0.99", values_temp = TRUE)
modern_nearly_0.99_REFisDER_final <- modern_nearly_0.99_REFisDER %>% 
	dplyr::select(names(modern_nearly_REFisANC)) %>% 
	mutate(names_temp = "hub_in_final_study_modern_nearly_0.99", values_temp = TRUE)
modern_nearly_0.999_REFisANC_final <- modern_nearly_0.999_REFisANC %>% 
	dplyr::select(names(modern_nearly_REFisANC)) %>% 
	mutate(names_temp = "hub_in_final_study_modern_nearly_0.999", values_temp = TRUE)
modern_nearly_0.999_REFisDER_final <- modern_nearly_0.999_REFisDER %>% 
	dplyr::select(names(modern_nearly_REFisANC)) %>% 
	mutate(names_temp = "hub_in_final_study_modern_nearly_0.999", values_temp = TRUE)

modern_fixed_REFisANC_final <- modern_fixed_REFisANC %>% 
	dplyr::select(names(modern_nearly_REFisANC)) %>% 
	mutate(names_temp = "hub_in_final_study_modern_fixed", values_temp = TRUE)
modern_fixed_REFisDER_final <- modern_fixed_REFisDER %>% 
	dplyr::select(names(modern_nearly_REFisANC)) %>% 
	mutate(names_temp = "hub_in_final_study_modern_fixed", values_temp = TRUE)
archaic_fixed_REFisANC_final <- archaic_fixed_REFisANC %>% 
	dplyr::select(names(modern_nearly_REFisANC)) %>% 
	mutate(names_temp = "hub_in_final_study_archaic_fixed", values_temp = TRUE)
archaic_fixed_REFisDER_final <- archaic_fixed_REFisDER %>% 
	dplyr::select(names(modern_nearly_REFisANC)) %>% 
	mutate(names_temp = "hub_in_final_study_archaic_fixed", values_temp = TRUE)
neanderthal_fixed_REFisANC_final <- neanderthal_fixed_REFisANC %>% 
	dplyr::select(names(modern_nearly_REFisANC)) %>% 
	mutate(names_temp = "hub_in_final_study_nean_fixed", values_temp = TRUE)
neanderthal_fixed_REFisDER_final <- neanderthal_fixed_REFisDER %>% 
	dplyr::select(names(modern_nearly_REFisANC)) %>% 
	mutate(names_temp = "hub_in_final_study_nean_fixed", values_temp = TRUE)
denisovan_fixed_REFisANC_final <- denisovan_fixed_REFisANC %>% 
	dplyr::select(names(modern_nearly_REFisANC)) %>% 
	mutate(names_temp = "hub_in_final_study_deni_fixed", values_temp = TRUE)
denisovan_fixed_REFisDER_final <- denisovan_fixed_REFisDER %>% 
	dplyr::select(names(modern_nearly_REFisANC)) %>% 
	mutate(names_temp = "hub_in_final_study_deni_fixed", values_temp = TRUE)

# save temp files
write_tsv(modern_nearly_REFisANC_final, gzfile("../../results/preprocess_1KGP_SNPs/modern_nearly_REFisANC_final.txt.gz"))
write_tsv(modern_nearly_REFisDER_final, gzfile("../../results/preprocess_1KGP_SNPs/modern_nearly_REFisDER_final.txt.gz"))
write_tsv(archaic_nearly_REFisANC_final, gzfile("../../results/preprocess_1KGP_SNPs/archaic_nearly_REFisANC_final.txt.gz"))
write_tsv(archaic_nearly_REFisDER_final, gzfile("../../results/preprocess_1KGP_SNPs/archaic_nearly_REFisDER_final.txt.gz"))
write_tsv(neanderthal_nearly_REFisANC_final, gzfile("../../results/preprocess_1KGP_SNPs/neanderthal_nearly_REFisANC_final.txt.gz"))
write_tsv(neanderthal_nearly_REFisDER_final, gzfile("../../results/preprocess_1KGP_SNPs/neanderthal_nearly_REFisDER_final.txt.gz"))
write_tsv(denisovan_nearly_REFisANC_final, gzfile("../../results/preprocess_1KGP_SNPs/denisovan_nearly_REFisANC_final.txt.gz"))
write_tsv(denisovan_nearly_REFisDER_final, gzfile("../../results/preprocess_1KGP_SNPs/denisovan_nearly_REFisDER_final.txt.gz"))

write_tsv(modern_nearly_0.5_REFisANC_final, gzfile("../../results/preprocess_1KGP_SNPs/modern_nearly_0.5_REFisANC_final.txt.gz"))
write_tsv(modern_nearly_0.5_REFisDER_final, gzfile("../../results/preprocess_1KGP_SNPs/modern_nearly_0.5_REFisDER_final.txt.gz"))
write_tsv(modern_nearly_0.6_REFisANC_final, gzfile("../../results/preprocess_1KGP_SNPs/modern_nearly_0.6_REFisANC_final.txt.gz"))
write_tsv(modern_nearly_0.6_REFisDER_final, gzfile("../../results/preprocess_1KGP_SNPs/modern_nearly_0.6_REFisDER_final.txt.gz"))
write_tsv(modern_nearly_0.7_REFisANC_final, gzfile("../../results/preprocess_1KGP_SNPs/modern_nearly_0.7_REFisANC_final.txt.gz"))
write_tsv(modern_nearly_0.7_REFisDER_final, gzfile("../../results/preprocess_1KGP_SNPs/modern_nearly_0.7_REFisDER_final.txt.gz"))
write_tsv(modern_nearly_0.8_REFisANC_final, gzfile("../../results/preprocess_1KGP_SNPs/modern_nearly_0.8_REFisANC_final.txt.gz"))
write_tsv(modern_nearly_0.8_REFisDER_final, gzfile("../../results/preprocess_1KGP_SNPs/modern_nearly_0.8_REFisDER_final.txt.gz"))
write_tsv(modern_nearly_0.9_REFisANC_final, gzfile("../../results/preprocess_1KGP_SNPs/modern_nearly_0.9_REFisANC_final.txt.gz"))
write_tsv(modern_nearly_0.9_REFisDER_final, gzfile("../../results/preprocess_1KGP_SNPs/modern_nearly_0.9_REFisDER_final.txt.gz"))
write_tsv(modern_nearly_0.99_REFisANC_final, gzfile("../../results/preprocess_1KGP_SNPs/modern_nearly_0.99_REFisANC_final.txt.gz"))
write_tsv(modern_nearly_0.99_REFisDER_final, gzfile("../../results/preprocess_1KGP_SNPs/modern_nearly_0.99_REFisDER_final.txt.gz"))
write_tsv(modern_nearly_0.999_REFisANC_final, gzfile("../../results/preprocess_1KGP_SNPs/modern_nearly_0.999_REFisANC_final.txt.gz"))
write_tsv(modern_nearly_0.999_REFisDER_final, gzfile("../../results/preprocess_1KGP_SNPs/modern_nearly_0.999_REFisDER_final.txt.gz"))

write_tsv(modern_fixed_REFisANC_final, gzfile("../../results/preprocess_1KGP_SNPs/modern_fixed_REFisANC_final.txt.gz"))
write_tsv(modern_fixed_REFisDER_final, gzfile("../../results/preprocess_1KGP_SNPs/modern_fixed_REFisDER_final.txt.gz"))
write_tsv(archaic_fixed_REFisANC_final, gzfile("../../results/preprocess_1KGP_SNPs/archaic_fixed_REFisANC_final.txt.gz"))
write_tsv(archaic_fixed_REFisDER_final, gzfile("../../results/preprocess_1KGP_SNPs/archaic_fixed_REFisDER_final.txt.gz"))
write_tsv(neanderthal_fixed_REFisANC_final, gzfile("../../results/preprocess_1KGP_SNPs/neanderthal_fixed_REFisANC_final.txt.gz"))
write_tsv(neanderthal_fixed_REFisDER_final, gzfile("../../results/preprocess_1KGP_SNPs/neanderthal_fixed_REFisDER_final.txt.gz"))
write_tsv(denisovan_fixed_REFisANC_final, gzfile("../../results/preprocess_1KGP_SNPs/denisovan_fixed_REFisANC_final.txt.gz"))
write_tsv(denisovan_fixed_REFisDER_final, gzfile("../../results/preprocess_1KGP_SNPs/denisovan_fixed_REFisDER_final.txt.gz"))
