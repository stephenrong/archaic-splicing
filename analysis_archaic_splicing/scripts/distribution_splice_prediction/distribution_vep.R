#!/bin/R

library(tidyverse)
library(data.table)
library(rstatix)
library(ggpubr)
library(Hmisc)
library(MetBrewer)

# load tables
mapsy_variant_table <- as_tibble(fread("../../results/mapsy_to_variant_table_updated/Neanderthal_updated_mapsy_to_variant_table.txt.gz"))
final_v2_variants_B_stat_mask_1KGP_archaic_gnomAD_hub_vep <- as_tibble(fread("../../results/vep_annotations_SNPs/final_v2_variants_B_stat_mask_1KGP_archaic_gnomAD_hub_vep.txt.gz"))

# join tables
mapsy_variant_table_vep <- left_join(mapsy_variant_table, final_v2_variants_B_stat_mask_1KGP_archaic_gnomAD_hub_vep)

# visualize in bins
mapsy_variant_table_vep <- mapsy_variant_table_vep %>% 
	mutate(mpralm.ANCDER.logFC.bin = cut_width(mpralm.ANCDER.logFC, 2, boundary=1))
write_tsv(mapsy_variant_table_vep, gzfile("../../results/distribution_splice_prediction/mapsy_variant_table_vep.txt.gz"))

# count VEP
mapsy_variant_table_vep_long <- mapsy_variant_table_vep %>% 
	pivot_longer(cols=c("stop_gained", "missense_variant", "synonymous_variant", "splice_region_variant", "5_prime_UTR_variant", "3_prime_UTR_variant"), names_to="VEP", values_to="VEP_TRUE") %>% 
	mutate(VEP=recode(VEP, stop_gained="Stop gain", missense_variant="Missense", synonymous_variant="Synonymous", splice_region_variant="Splice region", `5_prime_UTR_variant`="5' UTR", `3_prime_UTR_variant`="3' UTR")) %>% 
	filter(VEP_TRUE)

# deduplicate
mapsy_variant_table_vep_long <- mapsy_variant_table_vep_long %>% 
	mutate(VEP=factor(VEP, levels=c("Splice region", "Stop gain", "Missense", "Synonymous", "5' UTR", "3' UTR"))) %>% 
	arrange(VEP) %>% 
	filter(!duplicated(hub_variant_ID))
write_tsv(mapsy_variant_table_vep_long, gzfile("../../results/distribution_splice_prediction/mapsy_variant_table_vep_long.txt.gz"))

# get counts and proportions
mapsy_variant_table_vep_long_summ <- mapsy_variant_table_vep_long %>% 
	group_by(mpralm.ANCDER.logFC.bin, VEP) %>% count()
mapsy_variant_table_vep_long_summ <- mapsy_variant_table_vep_long_summ %>% 
	group_by(mpralm.ANCDER.logFC.bin) %>% 
	mutate(p = n/sum(n)) %>% 
	ungroup()
write_tsv(mapsy_variant_table_vep_long_summ, gzfile("../../results/distribution_splice_prediction/mapsy_variant_table_vep_long_summ.txt.gz"))

# visualize
ggplot(mapsy_variant_table_vep_long_summ %>% filter(!is.na(mpralm.ANCDER.logFC.bin))) + 
	geom_bar(aes(x=mpralm.ANCDER.logFC.bin, y=n, fill=VEP), stat="identity") + 
	xlab(expression(paste("MaPSy functional score (log"[2], " FC) bin"))) + ylab("VEP count") + scale_fill_met_d("Homer2") + theme(aspect.ratio=0.3) + 
	theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
ggsave("../../results/distribution_splice_prediction/distribution_vep_count.pdf", scale=0.7)

ggplot(mapsy_variant_table_vep_long_summ %>% filter(!is.na(mpralm.ANCDER.logFC.bin))) + 
	geom_bar(aes(x=mpralm.ANCDER.logFC.bin, y=n, fill=VEP), position="fill", stat="identity") + 
	xlab(expression(paste("MaPSy functional score (log"[2], " FC) bin"))) + ylab("VEP proportion") + scale_fill_met_d("Homer2") + theme(aspect.ratio=0.3) + 
	theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
ggsave("../../results/distribution_splice_prediction/distribution_vep_proportion.pdf", scale=0.7)
