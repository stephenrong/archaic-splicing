#!/bin/R

# Append Rinker et al NDA, RAA, RAH SNPs to final

# load packages
library(tidyverse)
library(data.table)

# load data
final_variants_B_stat_mask_1KGP_archaic_gnomAD_hub <- as_tibble(fread("../../results/preprocess_1KGP_SNPs/final_variants_B_stat_mask_1KGP_archaic_gnomAD_hub.txt.gz"))
Rinker_et_al_join_B_stat_mask_1KGP_archaic_gnomAD_hub <- as_tibble(fread("../../results/preprocess_1KGP_SNPs/Rinker_et_al_join_B_stat_mask_1KGP_archaic_gnomAD_hub.txt.gz"))

# append
final_v2_variants_B_stat_mask_1KGP_archaic_gnomAD_hub <- 
	full_join(
		final_variants_B_stat_mask_1KGP_archaic_gnomAD_hub %>% 
			mutate(hub_variant_CHROM = as.character(hub_variant_CHROM)) %>% 
			mutate(CHROM = as.character(CHROM)) %>% 
			mutate(gnomAD_AF = gnomAD_AC/gnomAD_AN),
		Rinker_et_al_join_B_stat_mask_1KGP_archaic_gnomAD_hub %>% 
			mutate(hub_variant_CHROM = as.character(hub_variant_CHROM)) %>% 
			mutate(CHROM = as.character(CHROM)) %>% 
			mutate(gnomAD_AF = gnomAD_AC/gnomAD_AN)
	) %>% 
	dplyr::select(starts_with("hub_"), !starts_with("hub_"))

# arrange
final_v2_variants_B_stat_mask_1KGP_archaic_gnomAD_hub <- final_v2_variants_B_stat_mask_1KGP_archaic_gnomAD_hub %>% 
	arrange(hub_variant_CHROM, hub_variant_POS, hub_variant_REF, hub_variant_ALT)

# fill na
temp <- final_v2_variants_B_stat_mask_1KGP_archaic_gnomAD_hub[,10:79]
temp[is.na(temp)] <- FALSE
final_v2_variants_B_stat_mask_1KGP_archaic_gnomAD_hub[,10:79] <- temp

# fix error, vernot not in ldvernot
final_v2_variants_B_stat_mask_1KGP_archaic_gnomAD_hub <- final_v2_variants_B_stat_mask_1KGP_archaic_gnomAD_hub %>% 
	mutate(hub_in_ldvernot_akey_2016 = (hub_in_ldvernot_akey_2016|hub_in_vernot_akey_2016)) %>% 
	mutate(hub_in_vernot_akey_2016_EAS = (hub_in_vernot_akey_2016_EAS|hub_in_vernot_akey_2016_EAS)) %>% 
	mutate(hub_in_vernot_akey_2016_EUR = (hub_in_vernot_akey_2016_EUR|hub_in_vernot_akey_2016_EUR)) %>% 
	mutate(hub_in_vernot_akey_2016_SAS = (hub_in_vernot_akey_2016_SAS|hub_in_vernot_akey_2016_SAS)) %>% 
	mutate(hub_in_vernot_akey_2016_MEL = (hub_in_vernot_akey_2016_MEL|hub_in_vernot_akey_2016_MEL))

# check overlaps
table(final_v2_variants_B_stat_mask_1KGP_archaic_gnomAD_hub[c("hub_in_ldvernot_akey_2016", "hub_in_vernot_akey_2016")])

table(final_v2_variants_B_stat_mask_1KGP_archaic_gnomAD_hub[c("hub_in_ldvernot_akey_2016", "hub_in_rinker_2020_NDA")])
table(final_v2_variants_B_stat_mask_1KGP_archaic_gnomAD_hub[c("hub_in_ldvernot_akey_2016", "hub_in_rinker_2020_RAA")])
table(final_v2_variants_B_stat_mask_1KGP_archaic_gnomAD_hub[c("hub_in_ldvernot_akey_2016", "hub_in_rinker_2020_RHA")])

table(final_v2_variants_B_stat_mask_1KGP_archaic_gnomAD_hub[c("hub_in_browning_2018", "hub_in_rinker_2020_NDA")])
table(final_v2_variants_B_stat_mask_1KGP_archaic_gnomAD_hub[c("hub_in_browning_2018", "hub_in_rinker_2020_RAA")])
table(final_v2_variants_B_stat_mask_1KGP_archaic_gnomAD_hub[c("hub_in_browning_2018", "hub_in_rinker_2020_RHA")])

# add another column for RA
final_v2_variants_B_stat_mask_1KGP_archaic_gnomAD_hub <- final_v2_variants_B_stat_mask_1KGP_archaic_gnomAD_hub %>% 
	mutate(hub_in_rinker_2020_RA = (hub_in_rinker_2020_RAA | hub_in_rinker_2020_RHA)) %>%  # all reintroduced
	mutate(hub_in_rinker_2020_ALL = (hub_in_rinker_2020_NDA | hub_in_rinker_2020_RAA | hub_in_rinker_2020_RHA))

# adjustment, introgressed should include rinker
final_v2_variants_B_stat_mask_1KGP_archaic_gnomAD_hub <- final_v2_variants_B_stat_mask_1KGP_archaic_gnomAD_hub %>% 
	mutate(hub_in_final_study_introgressed = (hub_in_ldvernot_akey_2016|hub_in_browning_2018))  # no Rinker et al, because they are derived from Vernot et al

# adjustment, introgressed not include racimo
final_v2_variants_B_stat_mask_1KGP_archaic_gnomAD_hub <- final_v2_variants_B_stat_mask_1KGP_archaic_gnomAD_hub %>% 
	mutate(hub_in_final_study_adaptive = (hub_in_gittelman_2016|hub_in_racimo_2017_AI))  # no Racimo et al, because ... 

# rearrange columns
final_v2_variants_B_stat_mask_1KGP_archaic_gnomAD_hub <- final_v2_variants_B_stat_mask_1KGP_archaic_gnomAD_hub %>% 
	dplyr::select(starts_with("hub"), !starts_with("hub")) 

# save
write_tsv(final_v2_variants_B_stat_mask_1KGP_archaic_gnomAD_hub, gzfile("../../results/preprocess_1KGP_SNPs/final_v2_variants_B_stat_mask_1KGP_archaic_gnomAD_hub.txt.gz"))
