#!/bin/R

library(tidyverse)
library(data.table)

# Cache
ukb_fine_mapping_final <- as_tibble(fread("../../data/finemap_overlap/ukb/release1.1/UKBB_94traits_release1.bed.gz"))
names(ukb_fine_mapping_final) <- c("chromosome", "start", "end", "variant", "rsid", "allele1", "allele2", "minorallele", "cohort", "model_marginal", "method", "trait", "region", "maf", "beta_marginal", "se_marginal", "chisq_marginal", "pip", "cs_id", "beta_posterior", "sd_posterior", "LD_HWE", "LD_SV")
ukb_fine_mapping_traits <- as_tibble(fread("../../data/finemap_overlap/ukb/release1.1/UKBB_94traits_release1.traits"))
ukb_fine_mapping_final <- ukb_fine_mapping_final %>% left_join(ukb_fine_mapping_traits)
ukb_fine_mapping_final <- ukb_fine_mapping_final %>% 
	mutate(hub_variant_ID = paste(gsub("chr", "", chromosome), "_", start+1, "_", allele1, "/", allele2, sep=""))
ukb_fine_mapping_final_pip <- ukb_fine_mapping_final %>% 
	filter(pip >= 0.05)
write_tsv(ukb_fine_mapping_final_pip, gzfile("../../results/additional_analyses_fine_map_overlap/ukb_fine_mapping_final_pip.txt.gz"))

# Overlap
ukb_fine_mapping_final_pip <- read_tsv("../../results/additional_analyses_fine_map_overlap/ukb_fine_mapping_final_pip.txt.gz") %>% 
	filter(pip >= 0.1)
mapsy_variant_table <- as_tibble(fread("../../results/mapsy_to_variant_table_updated/Neanderthal_updated_mapsy_to_variant_table.txt.gz"))
final_v2_variants_table <- as_tibble(fread("../../results/preprocess_1KGP_SNPs/final_v2_variants_B_stat_mask_1KGP_archaic_gnomAD_introgr_hub_spliceai_gs_raw_dedup.txt.gz"))

mapsy_variant_table_finemap_ukb <- mapsy_variant_table %>% 
	filter(mpralm.sigvar) %>% 
	filter(hub_variant_ID %in% ukb_fine_mapping_final_pip$hub_variant_ID) %>% 
	left_join(ukb_fine_mapping_final_pip)
mapsy_variant_table_finemap_ukb_summ <- enframe(colSums(mapsy_variant_table_finemap_ukb[c("hub_in_final_study_modern", "hub_in_final_study_archaic", "hub_in_final_study_nean", "hub_in_final_study_deni", "hub_in_final_study_introgressed", "hub_in_final_study_adaptive")]))

final_v2_variants_table_finemap_ukb <- final_v2_variants_table %>% 
	filter(spliceai_max >= 0.1) %>% 
	filter(hub_variant_ID %in% ukb_fine_mapping_final_pip$hub_variant_ID) %>% 
	left_join(ukb_fine_mapping_final_pip)
final_v2_variants_table_finemap_ukb_summ <- enframe(colSums(final_v2_variants_table_finemap_ukb[c("hub_in_final_study_modern", "hub_in_final_study_archaic", "hub_in_final_study_nean", "hub_in_final_study_deni", "hub_in_final_study_introgressed", "hub_in_final_study_adaptive")]))

write_tsv(mapsy_variant_table_finemap_ukb, gzfile("../../results/additional_analyses_fine_map_overlap/mapsy_variant_table_finemap_ukb.txt.gz"))
write_tsv(final_v2_variants_table_finemap_ukb, gzfile("../../results/additional_analyses_fine_map_overlap/final_v2_variants_table_finemap_ukb.txt.gz"))

write_tsv(mapsy_variant_table_finemap_ukb_summ, "../../results/additional_analyses_fine_map_overlap/mapsy_variant_table_finemap_ukb_summ.txt")
write_tsv(final_v2_variants_table_finemap_ukb_summ, "../../results/additional_analyses_fine_map_overlap/final_v2_variants_table_finemap_ukb_summ.txt")
