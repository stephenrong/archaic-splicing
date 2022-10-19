#!/bin/R

#!/bin/R
library(tidyverse)
library(data.table)
library(UpSetR)

# variant data
variant_table <- as_tibble(fread("../../results/preprocess_1KGP_SNPs/final_v2_variants_B_stat_mask_1KGP_archaic_gnomAD_introgr_hub.txt.gz")) %>% 
	mutate(hub_in_final_study_introgressed_ANC_FAIL = ((hub_variant_introgr_class == "introgr_ANC_FAIL")&hub_in_final_study_adaptive)) %>%   # add col, not used elsewhere
	mutate(hub_in_final_study_adaptive_ANC_FAIL = ((hub_variant_introgr_class == "introgr_ANC_FAIL")&hub_in_final_study_adaptive))  # add col, not used elsewhere

# variant upset
variant_table_short_overview <- variant_table %>% dplyr::select(c("hub_in_final_study_modern", "hub_in_final_study_nean", "hub_in_final_study_deni", "hub_in_final_study_archaic", "hub_in_final_study_introgressed", "hub_in_final_study_adaptive")) %>% mutate_all(function(x) as.numeric(x))
names(variant_table_short_overview) <- c("Modern specific", "Neanderthal specific", "Denisovan specific", "Archaic specific", "Introgressed", "Adaptively introgressed")
pdf(file="../../results/upset_plot_variants/variant_table_short_overview.pdf", width=6, height=4, onefile=F)
upset(as.data.frame(variant_table_short_overview), nsets=6)
dev.off()

# mapsy data
mapsy_table <- as_tibble(fread("../../results/mapsy_to_variant_table_updated/Neanderthal_updated_mapsy_to_variant_table.txt.gz")) %>% 
	mutate(hub_in_final_study_introgressed_ANC_FAIL = ((hub_variant_introgr_class == "introgr_ANC_FAIL")&hub_in_final_study_adaptive)) %>%   # add col, not used elsewhere
	mutate(hub_in_final_study_adaptive_ANC_FAIL = ((hub_variant_introgr_class == "introgr_ANC_FAIL")&hub_in_final_study_adaptive))  # add col, not used elsewhere

# mapsy upset
mapsy_table_short_overview <- mapsy_table %>% dplyr::select(c("hub_in_final_study_modern", "hub_in_final_study_nean", "hub_in_final_study_deni", "hub_in_final_study_archaic", "hub_in_final_study_introgressed", "hub_in_final_study_adaptive")) %>% mutate_all(function(x) as.numeric(x))
names(mapsy_table_short_overview) <- c("Modern specific", "Neanderthal specific", "Denisovan specific", "Archaic specific", "Introgressed", "Adaptively introgressed")
pdf(file="../../results/upset_plot_variants/mapsy_table_short_overview.pdf", width=6, height=4, onefile=F)
upset(as.data.frame(mapsy_table_short_overview), nsets=6)
dev.off()
