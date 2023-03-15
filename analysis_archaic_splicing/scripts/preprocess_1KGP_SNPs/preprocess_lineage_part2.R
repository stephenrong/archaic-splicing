#!/bin/sh

# Reprocess to get updated lineage-specific variants and combine with old library data

library(tidyverse)
library(data.table)

merge_variants_B_stat_mask_1KGP_archaic_gnomAD_hub <- as_tibble(fread("../../results/preprocess_1KGP_SNPs/merge_variants_B_stat_mask_1KGP_archaic_gnomAD_hub.txt.gz")) %>% 
	mutate(across(everything(), as.character))

# load files
modern_nearly_REFisANC_final <- as_tibble(fread("../../results/preprocess_1KGP_SNPs/modern_nearly_REFisANC_final.txt.gz")) %>% 
	mutate(across(everything(), as.character))
modern_nearly_REFisDER_final <- as_tibble(fread("../../results/preprocess_1KGP_SNPs/modern_nearly_REFisDER_final.txt.gz")) %>% 
	mutate(across(everything(), as.character))
archaic_nearly_REFisANC_final <- as_tibble(fread("../../results/preprocess_1KGP_SNPs/archaic_nearly_REFisANC_final.txt.gz")) %>% 
	mutate(across(everything(), as.character))
archaic_nearly_REFisDER_final <- as_tibble(fread("../../results/preprocess_1KGP_SNPs/archaic_nearly_REFisDER_final.txt.gz")) %>% 
	mutate(across(everything(), as.character))
neanderthal_nearly_REFisANC_final <- as_tibble(fread("../../results/preprocess_1KGP_SNPs/neanderthal_nearly_REFisANC_final.txt.gz")) %>% 
	mutate(across(everything(), as.character))
neanderthal_nearly_REFisDER_final <- as_tibble(fread("../../results/preprocess_1KGP_SNPs/neanderthal_nearly_REFisDER_final.txt.gz")) %>% 
	mutate(across(everything(), as.character))
denisovan_nearly_REFisANC_final <- as_tibble(fread("../../results/preprocess_1KGP_SNPs/denisovan_nearly_REFisANC_final.txt.gz")) %>% 
	mutate(across(everything(), as.character))
denisovan_nearly_REFisDER_final <- as_tibble(fread("../../results/preprocess_1KGP_SNPs/denisovan_nearly_REFisDER_final.txt.gz")) %>% 
	mutate(across(everything(), as.character))
modern_nearly_0.5_REFisANC_final <- as_tibble(fread("../../results/preprocess_1KGP_SNPs/modern_nearly_0.5_REFisANC_final.txt.gz")) %>% 
	mutate(across(everything(), as.character))
modern_nearly_0.5_REFisDER_final <- as_tibble(fread("../../results/preprocess_1KGP_SNPs/modern_nearly_0.5_REFisDER_final.txt.gz")) %>% 
	mutate(across(everything(), as.character))
modern_nearly_0.6_REFisANC_final <- as_tibble(fread("../../results/preprocess_1KGP_SNPs/modern_nearly_0.6_REFisANC_final.txt.gz")) %>% 
	mutate(across(everything(), as.character))
modern_nearly_0.6_REFisDER_final <- as_tibble(fread("../../results/preprocess_1KGP_SNPs/modern_nearly_0.6_REFisDER_final.txt.gz")) %>% 
	mutate(across(everything(), as.character))
modern_nearly_0.7_REFisANC_final <- as_tibble(fread("../../results/preprocess_1KGP_SNPs/modern_nearly_0.7_REFisANC_final.txt.gz")) %>% 
	mutate(across(everything(), as.character))
modern_nearly_0.7_REFisDER_final <- as_tibble(fread("../../results/preprocess_1KGP_SNPs/modern_nearly_0.7_REFisDER_final.txt.gz")) %>% 
	mutate(across(everything(), as.character))
modern_nearly_0.8_REFisANC_final <- as_tibble(fread("../../results/preprocess_1KGP_SNPs/modern_nearly_0.8_REFisANC_final.txt.gz")) %>% 
	mutate(across(everything(), as.character))
modern_nearly_0.8_REFisDER_final <- as_tibble(fread("../../results/preprocess_1KGP_SNPs/modern_nearly_0.8_REFisDER_final.txt.gz")) %>% 
	mutate(across(everything(), as.character))
modern_nearly_0.9_REFisANC_final <- as_tibble(fread("../../results/preprocess_1KGP_SNPs/modern_nearly_0.9_REFisANC_final.txt.gz")) %>% 
	mutate(across(everything(), as.character))
modern_nearly_0.9_REFisDER_final <- as_tibble(fread("../../results/preprocess_1KGP_SNPs/modern_nearly_0.9_REFisDER_final.txt.gz")) %>% 
	mutate(across(everything(), as.character))
modern_nearly_0.99_REFisANC_final <- as_tibble(fread("../../results/preprocess_1KGP_SNPs/modern_nearly_0.99_REFisANC_final.txt.gz")) %>% 
	mutate(across(everything(), as.character))
modern_nearly_0.99_REFisDER_final <- as_tibble(fread("../../results/preprocess_1KGP_SNPs/modern_nearly_0.99_REFisDER_final.txt.gz")) %>% 
	mutate(across(everything(), as.character))
modern_nearly_0.999_REFisANC_final <- as_tibble(fread("../../results/preprocess_1KGP_SNPs/modern_nearly_0.999_REFisANC_final.txt.gz")) %>% 
	mutate(across(everything(), as.character))
modern_nearly_0.999_REFisDER_final <- as_tibble(fread("../../results/preprocess_1KGP_SNPs/modern_nearly_0.999_REFisDER_final.txt.gz")) %>% 
	mutate(across(everything(), as.character))
modern_fixed_REFisANC_final <- as_tibble(fread("../../results/preprocess_1KGP_SNPs/modern_fixed_REFisANC_final.txt.gz")) %>% 
	mutate(across(everything(), as.character))
modern_fixed_REFisDER_final <- as_tibble(fread("../../results/preprocess_1KGP_SNPs/modern_fixed_REFisDER_final.txt.gz")) %>% 
	mutate(across(everything(), as.character))
archaic_fixed_REFisANC_final <- as_tibble(fread("../../results/preprocess_1KGP_SNPs/archaic_fixed_REFisANC_final.txt.gz")) %>% 
	mutate(across(everything(), as.character))
archaic_fixed_REFisDER_final <- as_tibble(fread("../../results/preprocess_1KGP_SNPs/archaic_fixed_REFisDER_final.txt.gz")) %>% 
	mutate(across(everything(), as.character))
neanderthal_fixed_REFisANC_final <- as_tibble(fread("../../results/preprocess_1KGP_SNPs/neanderthal_fixed_REFisANC_final.txt.gz")) %>% 
	mutate(across(everything(), as.character))
neanderthal_fixed_REFisDER_final <- as_tibble(fread("../../results/preprocess_1KGP_SNPs/neanderthal_fixed_REFisDER_final.txt.gz")) %>% 
	mutate(across(everything(), as.character))
denisovan_fixed_REFisANC_final <- as_tibble(fread("../../results/preprocess_1KGP_SNPs/denisovan_fixed_REFisANC_final.txt.gz")) %>% 
	mutate(across(everything(), as.character))
denisovan_fixed_REFisDER_final <- as_tibble(fread("../../results/preprocess_1KGP_SNPs/denisovan_fixed_REFisDER_final.txt.gz")) %>% 
	mutate(across(everything(), as.character))

# join files
join_variants_B_stat_mask_1KGP_archaic_gnomAD_hub <- bind_rows(
	modern_nearly_REFisANC_final, modern_nearly_REFisDER_final, 
	archaic_nearly_REFisANC_final, archaic_nearly_REFisDER_final, 
	neanderthal_nearly_REFisANC_final, neanderthal_nearly_REFisDER_final, 
	denisovan_nearly_REFisANC_final, denisovan_nearly_REFisDER_final, 
	modern_nearly_0.5_REFisANC_final, modern_nearly_0.5_REFisDER_final, 
	modern_nearly_0.6_REFisANC_final, modern_nearly_0.6_REFisDER_final, 
	modern_nearly_0.7_REFisANC_final, modern_nearly_0.7_REFisDER_final, 
	modern_nearly_0.8_REFisANC_final, modern_nearly_0.8_REFisDER_final, 
	modern_nearly_0.9_REFisANC_final, modern_nearly_0.9_REFisDER_final, 
	modern_nearly_0.99_REFisANC_final, modern_nearly_0.99_REFisDER_final, 
	modern_nearly_0.999_REFisANC_final, modern_nearly_0.999_REFisDER_final, 
	modern_fixed_REFisANC_final, modern_fixed_REFisDER_final, 
	archaic_fixed_REFisANC_final, archaic_fixed_REFisDER_final, 
	neanderthal_fixed_REFisANC_final, neanderthal_fixed_REFisDER_final, 
	denisovan_fixed_REFisANC_final, denisovan_fixed_REFisDER_final
)

join_variants_B_stat_mask_1KGP_archaic_gnomAD_hub <- 
	join_variants_B_stat_mask_1KGP_archaic_gnomAD_hub %>% 
	pivot_wider(names_from = names_temp, values_from = values_temp) %>% 
	dplyr::select(starts_with("hub"), !starts_with("hub"))
write_tsv(join_variants_B_stat_mask_1KGP_archaic_gnomAD_hub, gzfile("../../results/preprocess_1KGP_SNPs/join_variants_B_stat_mask_1KGP_archaic_gnomAD_hub.txt.gz"))
join_variants_B_stat_mask_1KGP_archaic_gnomAD_hub <- as_tibble(fread("../../results/preprocess_1KGP_SNPs/join_variants_B_stat_mask_1KGP_archaic_gnomAD_hub.txt.gz"))%>% 
	mutate(across(everything(), as.character))

# rejoin to variant table and check
preprocess_join <- function(tb) { 
	# recalculate gnomAD_AF before join
	tb <- tb %>% 
		mutate(gnomAD_AF = as.character(as.numeric(gnomAD_AC)/as.numeric(gnomAD_AN)))
}

# save final table
final_variants_B_stat_mask_1KGP_archaic_gnomAD_hub <- 
	full_join(
		preprocess_join(merge_variants_B_stat_mask_1KGP_archaic_gnomAD_hub), 
		preprocess_join(join_variants_B_stat_mask_1KGP_archaic_gnomAD_hub)
	) %>% 
	mutate(hub_in_final_study_introgressed = hub_in_this_study_introgressed) %>% 
	mutate(hub_in_final_study_adaptive = hub_in_this_study_adaptive) %>% 
	dplyr::select(starts_with("hub"), !starts_with("hub")) %>% 
	mutate(hub_variant_CHROM = factor(hub_variant_CHROM, levels=c(as.character(c(1:22)), "X"))) %>% 
	arrange(hub_variant_CHROM, as.numeric(hub_variant_POS), hub_variant_REF, hub_variant_ALT)
write_tsv(final_variants_B_stat_mask_1KGP_archaic_gnomAD_hub, gzfile("../../results/preprocess_1KGP_SNPs/final_variants_B_stat_mask_1KGP_archaic_gnomAD_hub.txt.gz"))

# save dup rows
final_variants_B_stat_mask_1KGP_archaic_gnomAD_dup_hub <- 
	final_variants_B_stat_mask_1KGP_archaic_gnomAD_hub %>% 
	filter(hub_variant_ID %in% 
		(final_variants_B_stat_mask_1KGP_archaic_gnomAD_hub %>% 
			mutate(temp_dup = duplicated(hub_variant_ID)) %>% 
			filter(temp_dup) %>% .$hub_variant_ID)) %>% 
	mutate(hub_variant_CHROM = factor(hub_variant_CHROM, levels=c(as.character(c(1:22)), "X"))) %>% 
	arrange(hub_variant_CHROM, as.numeric(hub_variant_POS), hub_variant_REF, hub_variant_ALT)
print(final_variants_B_stat_mask_1KGP_archaic_gnomAD_dup_hub)
write_tsv(final_variants_B_stat_mask_1KGP_archaic_gnomAD_dup_hub, gzfile("../../results/preprocess_1KGP_SNPs/final_variants_B_stat_mask_1KGP_archaic_gnomAD_dup_hub.txt.gz"))

# save new rows
final_variants_B_stat_mask_1KGP_archaic_gnomAD_new_hub <- 
	final_variants_B_stat_mask_1KGP_archaic_gnomAD_hub %>% 
	filter(
		xor(as.logical(hub_in_this_study_modern), as.logical(hub_in_final_study_modern))|
		xor(as.logical(hub_in_this_study_archaic), as.logical(hub_in_final_study_archaic))|
		xor(as.logical(hub_in_this_study_nean), as.logical(hub_in_final_study_nean))|
		xor(as.logical(hub_in_this_study_deni), as.logical(hub_in_final_study_deni))
	)
write_tsv(final_variants_B_stat_mask_1KGP_archaic_gnomAD_new_hub, gzfile("../../results/preprocess_1KGP_SNPs/final_variants_B_stat_mask_1KGP_archaic_gnomAD_new_hub.txt.gz"))
