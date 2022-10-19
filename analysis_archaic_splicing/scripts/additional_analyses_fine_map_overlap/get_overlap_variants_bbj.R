#!/bin/R

library(tidyverse)
library(data.table)

# Cache
bbj_fine_mapping_FINEMAP_final <- as_tibble(fread("../../results/additional_analyses_fine_map_overlap/bbj_fine_mapping_FINEMAP_final.txt.gz", sep="\t"))
bbj_fine_mapping_FINEMAP_traits <- as_tibble(fread("../../data/finemap_overlap/bbj/ST2_overview_traits_extract.txt"))
bbj_fine_mapping_FINEMAP_final <- bbj_fine_mapping_FINEMAP_final %>% left_join(bbj_fine_mapping_FINEMAP_traits)
bbj_fine_mapping_FINEMAP_final <- bbj_fine_mapping_FINEMAP_final %>% 
	mutate(hub_variant_ID = paste(gsub("chr", "", chromosome), "_", position, "_", allele1, "/", allele2, sep=""))
bbj_fine_mapping_FINEMAP_final_pip <- bbj_fine_mapping_FINEMAP_final %>% 
	filter(pip >= 0.05)
write_tsv(bbj_fine_mapping_FINEMAP_final_pip, gzfile("../../results/additional_analyses_fine_map_overlap/bbj_fine_mapping_FINEMAP_final_pip.txt.gz"))

bbj_fine_mapping_SUSIE_final <- as_tibble(fread("../../results/additional_analyses_fine_map_overlap/bbj_fine_mapping_SUSIE_final.txt.gz", sep="\t"))
bbj_fine_mapping_SUSIE_traits <- as_tibble(fread("../../data/finemap_overlap/bbj/ST2_overview_traits_extract.txt"))
bbj_fine_mapping_SUSIE_final <- bbj_fine_mapping_SUSIE_final %>% left_join(bbj_fine_mapping_SUSIE_traits)
bbj_fine_mapping_SUSIE_final <- bbj_fine_mapping_SUSIE_final %>% 
	mutate(hub_variant_ID = paste(gsub("chr", "", chromosome), "_", position, "_", allele1, "/", allele2, sep=""))
bbj_fine_mapping_SUSIE_final_pip <- bbj_fine_mapping_SUSIE_final %>% 
	filter(pip >= 0.05)
write_tsv(bbj_fine_mapping_SUSIE_final_pip, gzfile("../../results/additional_analyses_fine_map_overlap/bbj_fine_mapping_SUSIE_final_pip.txt.gz"))

# Overlap
bbj_fine_mapping_FINEMAP_final_pip <- read_tsv("../../results/additional_analyses_fine_map_overlap/bbj_fine_mapping_FINEMAP_final_pip.txt.gz") %>% 
	filter(pip >= 0.1)
mapsy_variant_table <- as_tibble(fread("../../results/mapsy_to_variant_table_updated/Neanderthal_updated_mapsy_to_variant_table.txt.gz"))
final_v2_variants_table <- as_tibble(fread("../../results/preprocess_1KGP_SNPs/final_v2_variants_B_stat_mask_1KGP_archaic_gnomAD_introgr_hub_spliceai_gs_raw_dedup.txt.gz"))

mapsy_variant_table_finemap_bbj_FINEMAP <- mapsy_variant_table %>% 
	filter(mpralm.sigvar) %>% 
	filter(hub_variant_ID %in% bbj_fine_mapping_FINEMAP_final_pip$hub_variant_ID) %>% 
	left_join(bbj_fine_mapping_FINEMAP_final_pip)
mapsy_variant_table_finemap_bbj_FINEMAP_summ <- enframe(colSums(mapsy_variant_table_finemap_bbj_FINEMAP[c("hub_in_final_study_modern", "hub_in_final_study_archaic", "hub_in_final_study_nean", "hub_in_final_study_deni", "hub_in_final_study_introgressed", "hub_in_final_study_adaptive")]))

final_v2_variants_table_finemap_bbj_FINEMAP <- final_v2_variants_table %>% 
	filter(spliceai_max >= 0.1) %>% 
	filter(hub_variant_ID %in% bbj_fine_mapping_FINEMAP_final_pip$hub_variant_ID) %>% 
	left_join(bbj_fine_mapping_FINEMAP_final_pip)
final_v2_variants_table_finemap_bbj_FINEMAP_summ <- enframe(colSums(final_v2_variants_table_finemap_bbj_FINEMAP[c("hub_in_final_study_modern", "hub_in_final_study_archaic", "hub_in_final_study_nean", "hub_in_final_study_deni", "hub_in_final_study_introgressed", "hub_in_final_study_adaptive")]))

write_tsv(mapsy_variant_table_finemap_bbj_FINEMAP, gzfile("../../results/additional_analyses_fine_map_overlap/mapsy_variant_table_finemap_bbj_FINEMAP.txt.gz"))
write_tsv(final_v2_variants_table_finemap_bbj_FINEMAP, gzfile("../../results/additional_analyses_fine_map_overlap/final_v2_variants_table_finemap_bbj_FINEMAP.txt.gz"))

write_tsv(mapsy_variant_table_finemap_bbj_FINEMAP_summ, "../../results/additional_analyses_fine_map_overlap/mapsy_variant_table_finemap_bbj_FINEMAP_summ.txt")
write_tsv(final_v2_variants_table_finemap_bbj_FINEMAP_summ, "../../results/additional_analyses_fine_map_overlap/final_v2_variants_table_finemap_bbj_FINEMAP_summ.txt")


bbj_fine_mapping_SUSIE_final_pip <- read_tsv("../../results/additional_analyses_fine_map_overlap/bbj_fine_mapping_SUSIE_final_pip.txt.gz") %>% 
	filter(pip >= 0.1)
mapsy_variant_table <- as_tibble(fread("../../results/mapsy_to_variant_table_updated/Neanderthal_updated_mapsy_to_variant_table.txt.gz"))
final_v2_variants_table <- as_tibble(fread("../../results/preprocess_1KGP_SNPs/final_v2_variants_B_stat_mask_1KGP_archaic_gnomAD_introgr_hub_spliceai_gs_raw_dedup.txt.gz"))

mapsy_variant_table_finemap_bbj_SUSIE <- mapsy_variant_table %>% 
	filter(mpralm.sigvar) %>% 
	filter(hub_variant_ID %in% bbj_fine_mapping_SUSIE_final_pip$hub_variant_ID) %>% 
	left_join(bbj_fine_mapping_SUSIE_final_pip)
mapsy_variant_table_finemap_bbj_SUSIE_summ <- enframe(colSums(mapsy_variant_table_finemap_bbj_SUSIE[c("hub_in_final_study_modern", "hub_in_final_study_archaic", "hub_in_final_study_nean", "hub_in_final_study_deni", "hub_in_final_study_introgressed", "hub_in_final_study_adaptive")]))

final_v2_variants_table_finemap_bbj_SUSIE <- final_v2_variants_table %>% 
	filter(spliceai_max >= 0.1) %>% 
	filter(hub_variant_ID %in% bbj_fine_mapping_SUSIE_final_pip$hub_variant_ID) %>% 
	left_join(bbj_fine_mapping_SUSIE_final_pip)
final_v2_variants_table_finemap_bbj_SUSIE_summ <- enframe(colSums(final_v2_variants_table_finemap_bbj_SUSIE[c("hub_in_final_study_modern", "hub_in_final_study_archaic", "hub_in_final_study_nean", "hub_in_final_study_deni", "hub_in_final_study_introgressed", "hub_in_final_study_adaptive")]))

write_tsv(mapsy_variant_table_finemap_bbj_SUSIE, gzfile("../../results/additional_analyses_fine_map_overlap/mapsy_variant_table_finemap_bbj_SUSIE.txt.gz"))
write_tsv(final_v2_variants_table_finemap_bbj_SUSIE, gzfile("../../results/additional_analyses_fine_map_overlap/final_v2_variants_table_finemap_bbj_SUSIE.txt.gz"))

write_tsv(mapsy_variant_table_finemap_bbj_SUSIE_summ, "../../results/additional_analyses_fine_map_overlap/mapsy_variant_table_finemap_bbj_SUSIE_summ.txt")
write_tsv(final_v2_variants_table_finemap_bbj_SUSIE_summ, "../../results/additional_analyses_fine_map_overlap/final_v2_variants_table_finemap_bbj_SUSIE_summ.txt")
