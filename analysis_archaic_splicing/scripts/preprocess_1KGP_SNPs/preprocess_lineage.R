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
archaics_all_1KGP_archaic_gnomAD_hub <- as_tibble(fread("../../results/preprocess_1KGP_SNPs/archaics_all_1KGP_archaic_gnomAD_hub.txt.gz"))
hg19_REFisDER_1KGP_archaic_gnomAD_hub <- as_tibble(fread("../../results/preprocess_1KGP_SNPs/hg19_REFisDER_1KGP_archaic_gnomAD_hub.txt.gz"))

# modern nearly:
# 	REFisANC, modern > 0.95, !(archaic >= 1), archaic_mask_inter
modern_nearly_REFisANC <- 
 	# add annotations, is a strict subset, so works for all SNPs, in particular gnomAD annotations
	ALL_1KGP_phase3_greater0.9_gnomAD_hub %>% 
	left_join(dplyr::select(merge_variants_B_stat_mask_1KGP_archaic_gnomAD_hub, !starts_with("hub_in"))) %>% 
	filter(hub_variant_REF == hub_variant_ANC) %>% 
	filter(archaic_mask_inter) %>% 
	filter((AF > 0.9) & (gnomAD_AF > 0.9)) %>% 
	filter(!(hub_variant_ID %in% filter(archaics_all_1KGP_archaic_gnomAD_hub, archaic_AC >= 1)$hub_variant_ID))
check_within(modern_nearly_REFisANC, merge_variants_B_stat_mask_1KGP_archaic_gnomAD_hub)

# 	REFisDER, modern < 0.05, archaic >= 8, archaic_mask_inter
modern_nearly_REFisDER <- 
	hg19_REFisDER_1KGP_archaic_gnomAD_hub %>% 
	filter(archaic_mask_inter) %>% 
	filter(archaic_AC >= 8) %>% 
	filter(((AF < 0.1) | (is.na(AF))) & ((gnomAD_AF < 0.1) | (is.na(gnomAD_AF))))
check_within(modern_nearly_REFisDER, merge_variants_B_stat_mask_1KGP_archaic_gnomAD_hub)
# a variant missing, was previously eliminated due to having high gnomAD frequency, now not eliminated due to not passing gnomAD filter

# archaic nearly:
# 	REFisANC, archaic >= 7, modern < 0.05, archaic_mask_inter
archaic_nearly_REFisANC <- 
	archaics_all_1KGP_archaic_gnomAD_hub %>% 
	filter(hub_variant_REF == hub_variant_ANC) %>% 
	filter(archaic_mask_inter) %>% 
	filter(archaic_AC >= 7) %>% 
	filter(((AFR_AF < 0.01) | (is.na(AFR_AF))) & ((gnomAD_AF_afr < 0.01) | (is.na(gnomAD_AF_afr))))
check_within(archaic_nearly_REFisANC, merge_variants_B_stat_mask_1KGP_archaic_gnomAD_hub)

# 	REFisDER, modern > 0.9, archaic <= 1, archaic_mask_inter
archaic_nearly_REFisDER <- 
	hg19_REFisDER_1KGP_archaic_gnomAD_hub %>% 
	filter(archaic_mask_inter) %>% 
	filter((archaic_AC <= 1) | (is.na(archaic_AC))) %>% 
	filter(((AFR_AF > 0.99) & (gnomAD_AF_afr > 0.99)) | ((is.na(AFR_AF) & (gnomAD_AF_afr > 0.99))) | ((AFR_AF > 0.99) & (is.na(gnomAD_AF_afr))))
check_within(archaic_nearly_REFisDER, merge_variants_B_stat_mask_1KGP_archaic_gnomAD_hub)
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
	filter(archaic_mask_inter) %>% 
	filter(neanderthal_AC <= 1, denisovan_AC == 2) %>% 
	filter(((AFR_AF > 0.99) & (gnomAD_AF_afr > 0.99)) | ((is.na(AFR_AF) & (gnomAD_AF_afr > 0.99))) | ((AFR_AF > 0.99) & (is.na(gnomAD_AF_afr))))
check_within(neanderthal_nearly_REFisDER, merge_variants_B_stat_mask_1KGP_archaic_gnomAD_hub)
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
	filter(archaic_mask_inter) %>% 
	filter(denisovan_AC == 0, neanderthal_AC == 6) %>% 
	filter(((AFR_AF > 0.99) & (gnomAD_AF_afr > 0.99)) | ((is.na(AFR_AF) & (gnomAD_AF_afr > 0.99))) | ((AFR_AF > 0.99) & (is.na(gnomAD_AF_afr))))
check_within(denisovan_nearly_REFisDER, merge_variants_B_stat_mask_1KGP_archaic_gnomAD_hub)
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

join_variants_B_stat_mask_1KGP_archaic_gnomAD_hub <- bind_rows(
		modern_nearly_REFisANC_final, modern_nearly_REFisDER_final, 
		archaic_nearly_REFisANC_final, archaic_nearly_REFisDER_final, 
		neanderthal_nearly_REFisANC_final, neanderthal_nearly_REFisDER_final, 
		denisovan_nearly_REFisANC_final, denisovan_nearly_REFisDER_final
	)

join_variants_B_stat_mask_1KGP_archaic_gnomAD_hub <- 
	join_variants_B_stat_mask_1KGP_archaic_gnomAD_hub %>% 
	pivot_wider(names_from = names_temp, values_from = values_temp) %>% 
	dplyr::select(starts_with("hub"), !starts_with("hub"))

# rejoin to variant table and check
preprocess_join <- function(tb) { 
	# recalculate gnomAD_AF before join
	tb <- tb %>% 
		mutate(gnomAD_AF = gnomAD_AC/gnomAD_AN)
}

final_variants_B_stat_mask_1KGP_archaic_gnomAD_hub <- 
	full_join(
		preprocess_join(merge_variants_B_stat_mask_1KGP_archaic_gnomAD_hub), 
		preprocess_join(join_variants_B_stat_mask_1KGP_archaic_gnomAD_hub)
	) %>% 
	mutate(hub_in_final_study_introgressed = hub_in_this_study_introgressed) %>% 
	mutate(hub_in_final_study_adaptive = hub_in_this_study_adaptive) %>% 
	dplyr::select(starts_with("hub"), !starts_with("hub")) %>% 
	arrange(hub_variant_CHROM, hub_variant_POS, hub_variant_REF, hub_variant_ALT)

# save final table
write_tsv(final_variants_B_stat_mask_1KGP_archaic_gnomAD_hub, gzfile("../../results/preprocess_1KGP_SNPs/final_variants_B_stat_mask_1KGP_archaic_gnomAD_hub.txt.gz"))

# save dup rows
final_variants_B_stat_mask_1KGP_archaic_gnomAD_dup_hub <- 
	final_variants_B_stat_mask_1KGP_archaic_gnomAD_hub %>% 
	filter(hub_variant_ID %in% 
		(final_variants_B_stat_mask_1KGP_archaic_gnomAD_hub %>% 
			mutate(temp_dup = duplicated(hub_variant_ID)) %>% 
			filter(temp_dup) %>% .$hub_variant_ID)) %>% 
	arrange(hub_variant_CHROM, hub_variant_POS, hub_variant_REF, hub_variant_ALT)
print(final_variants_B_stat_mask_1KGP_archaic_gnomAD_dup_hub)
write_tsv(final_variants_B_stat_mask_1KGP_archaic_gnomAD_dup_hub, gzfile("../../results/preprocess_1KGP_SNPs/final_variants_B_stat_mask_1KGP_archaic_gnomAD_dup_hub.txt.gz"))

# save new rows
final_variants_B_stat_mask_1KGP_archaic_gnomAD_new_hub <- 
	final_variants_B_stat_mask_1KGP_archaic_gnomAD_hub %>% 
	filter(
		xor(hub_in_this_study_modern, hub_in_final_study_modern)|
		xor(hub_in_this_study_archaic, hub_in_final_study_archaic)|
		xor(hub_in_this_study_nean, hub_in_final_study_nean)|
		xor(hub_in_this_study_deni, hub_in_final_study_deni)
	)
write_tsv(final_variants_B_stat_mask_1KGP_archaic_gnomAD_new_hub, gzfile("../../results/preprocess_1KGP_SNPs/final_variants_B_stat_mask_1KGP_archaic_gnomAD_new_hub.txt.gz"))
