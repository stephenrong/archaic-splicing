#!/bin/R

library(tidyverse)
library(data.table)
library(rstatix)
library(ggpubr)
library(Hmisc)

# load tables
mapsy_variant_table <- as_tibble(fread("../../results/mapsy_to_variant_table_updated/Neanderthal_updated_mapsy_to_variant_table.txt.gz"))
final_v2_variants_B_stat_mask_1KGP_archaic_gnomAD_introgr_hub_spliceai_gs_raw <- as_tibble(fread("../../results/preprocess_1KGP_SNPs/final_v2_variants_B_stat_mask_1KGP_archaic_gnomAD_introgr_hub_spliceai_gs_raw_dedup.txt.gz"))

# spliceai effect bin
final_v2_variants_B_stat_mask_1KGP_archaic_gnomAD_introgr_hub_spliceai_gs_raw <- final_v2_variants_B_stat_mask_1KGP_archaic_gnomAD_introgr_hub_spliceai_gs_raw %>% 
	mutate(spliceai_max = ifelse(is.na(spliceai_max), 0, spliceai_max))

# join tables
mapsy_variant_table_spliceai <- left_join(mapsy_variant_table, final_v2_variants_B_stat_mask_1KGP_archaic_gnomAD_introgr_hub_spliceai_gs_raw)

# visualize in bins
mapsy_variant_table_spliceai <- mapsy_variant_table_spliceai %>% 
	mutate(mpralm.ANCDER.logFC.bin = cut_width(mpralm.ANCDER.logFC, 2, boundary=1))
write_tsv(mapsy_variant_table_spliceai, gzfile("../../results/distribution_splice_prediction/mapsy_variant_table_spliceai_gs_raw.txt.gz"))

ggplot(mapsy_variant_table_spliceai %>% filter(!is.na(mpralm.ANCDER.logFC.bin))) + 
	geom_boxplot(aes(x=mpralm.ANCDER.logFC.bin, y=spliceai_max), fill="#beaed4") + 
	theme_pubr() + xlab(expression(paste("MaPSy functional score (log"[2], " FC) bin"))) + ylab("SpliceAI max score") + theme(aspect.ratio=0.3) + 
	theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
ggsave("../../results/distribution_splice_prediction/distribution_spliceai_gs_raw.pdf", scale=0.55)
