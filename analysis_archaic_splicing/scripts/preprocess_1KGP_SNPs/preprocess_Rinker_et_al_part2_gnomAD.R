#!/bin/R

# Join gnomAD SNP info to Rinker data

library(tidyverse)
library(data.table)
library(plyranges)
library(wrapr)
library(vcfR)
library(rtracklayer)
source("../get_helper.R")

# gnomad sites
# 	temp, remove
Rinker_et_al_join_B_stat_mask_1KGP_archaic_hub <- as_tibble(fread("../../results/preprocess_1KGP_SNPs/Rinker_et_al_join_B_stat_mask_1KGP_archaic_hub.txt.gz"))

# 	vcf to table
print("gnomAD vcf to table")
col_list <- c("AC", "AF", "AN", "AC_afr", "AF_afr", "AN_afr")
Rinker_et_al_join_B_stat_mask_1KGP_archaic_gnomAD_file <- "../../results/preprocess_1KGP_SNPs/Rinker_et_al_join_B_stat_mask_1KGP_archaic_gnomAD_file.vcf.gz"
Rinker_et_al_join_B_stat_mask_1KGP_archaic_gnomAD_input <- read.vcfR(Rinker_et_al_join_B_stat_mask_1KGP_archaic_gnomAD_file) %>% 
	.@fix %>% as_tibble() %>% unique() %>% 
	# filter(  # unnecessary, because of later left join
	# 	REF %in% c("A", "C", "T", "G"),
	# 	ALT %in% c("A", "C", "T", "G")
	# ) %>% 
	separate(col=INFO, into=col_list, sep=";")
for (col in col_list) {
	let(c(VAR=col), 
		Rinker_et_al_join_B_stat_mask_1KGP_archaic_gnomAD_input <- Rinker_et_al_join_B_stat_mask_1KGP_archaic_gnomAD_input %>% 
			mutate(VAR = as.numeric(gsub(".*=", "", VAR))))
}
names(Rinker_et_al_join_B_stat_mask_1KGP_archaic_gnomAD_input) <- 
	c("hub_variant_CHROM", "hub_variant_POS", "gnomAD_ID", 
		"hub_variant_REF", "hub_variant_ALT", "gnomAD_QUAL", 
		"gnomAD_FILTER", "gnomAD_AC", "gnomAD_AN", "gnomAD_AF", 
		"gnomAD_AC_afr", "gnomAD_AN_afr", "gnomAD_AF_afr")

Rinker_et_al_join_B_stat_mask_1KGP_archaic_gnomAD_input <- Rinker_et_al_join_B_stat_mask_1KGP_archaic_gnomAD_input %>% 
	mutate(hub_variant_POS = as.numeric(hub_variant_POS)) %>% 
	mutate_hub_variant_ID()  # order doesn't matter, don't dplyr::select

# 	join to table
print("gnomAD join to table")
Rinker_et_al_join_B_stat_mask_1KGP_archaic_gnomAD_hub <- 
	left_join(
		Rinker_et_al_join_B_stat_mask_1KGP_archaic_hub %>% mutate(hub_variant_CHROM = as.character(hub_variant_CHROM)), 
		Rinker_et_al_join_B_stat_mask_1KGP_archaic_gnomAD_input %>% mutate(hub_variant_CHROM = as.character(hub_variant_CHROM))
	)

# 	save file
print("gnomAD save file")
write_tsv(Rinker_et_al_join_B_stat_mask_1KGP_archaic_gnomAD_hub, gzfile("../../results/preprocess_1KGP_SNPs/Rinker_et_al_join_B_stat_mask_1KGP_archaic_gnomAD_hub.txt.gz"))

# 	check duplicates
print(table(table(Rinker_et_al_join_B_stat_mask_1KGP_archaic_gnomAD_hub$hub_variant_ID)))
