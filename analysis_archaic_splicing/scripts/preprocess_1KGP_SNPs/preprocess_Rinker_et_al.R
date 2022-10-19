#!/bin/R

# Add Rinker et al NDA, RAA, RAH SNPs (use 1KGP SNP data to get AC, get B, get mask, convert to hub format, add archaic SNP info, add gnomAD SNP info [convert to BED, tabix gnomAD, join gnomAD])

# load packages
library(tidyverse)
library(data.table)
source("../get_helper.R")

# load Rinker et al
Rinker_et_al_EUR <- as_tibble(fread("../../data/annotate_rinker_et_al/Rinker_et_al_EUR.txt")) %>% mutate(POP = "EUR")
Rinker_et_al_EAS <- as_tibble(fread("../../data/annotate_rinker_et_al/Rinker_et_al_EAS.txt")) %>% mutate(POP = "EAS")
Rinker_et_al_SAS <- as_tibble(fread("../../data/annotate_rinker_et_al/Rinker_et_al_SAS.txt")) %>% mutate(POP = "SAS")

# process Rinker et al EUR
Rinker_et_al_EUR <- Rinker_et_al_EUR %>% 
	mutate(
		hub_reference_genome = "hg19", 
		hub_variant_CHROM = gsub("chr", "", CHROM), 
		hub_variant_POS = POS
	)
Rinker_et_al_EUR <- Rinker_et_al_EUR %>% 
	mutate(hub_in_rinker_2020_NDA_EUR = ifelse(CLASS=="NDA", TRUE, FALSE)) %>% 
	mutate(hub_in_rinker_2020_RAA_EUR = ifelse(CLASS=="RAA", TRUE, FALSE)) %>% 
	mutate(hub_in_rinker_2020_RHA_EUR = ifelse(CLASS=="RHA", TRUE, FALSE)) %>% 
	mutate(hub_in_rinker_2020_RA_EUR = (hub_in_rinker_2020_RAA_EUR | hub_in_rinker_2020_RHA_EUR))
Rinker_et_al_EUR <- Rinker_et_al_EUR %>% 
	dplyr::select(starts_with("hub_"), !starts_with("hub_")) %>% 
	dplyr::rename(EUR_INTROG_ALLELE=INTROG_ALLELE, EUR_INTROG_AF=INTROG_AF, EUR_INTROG_V16_LD_haplotype=V16_LD_haplotype) %>% 
	dplyr::select(-c("CHROM", "POS", "CLASS", "POP"))

# process Rinker et al EAS
Rinker_et_al_EAS <- Rinker_et_al_EAS %>% 
	mutate(
		hub_reference_genome = "hg19", 
		hub_variant_CHROM = gsub("chr", "", CHROM), 
		hub_variant_POS = POS
	)
Rinker_et_al_EAS <- Rinker_et_al_EAS %>% 
	mutate(hub_in_rinker_2020_NDA_EAS = ifelse(CLASS=="NDA", TRUE, FALSE)) %>% 
	mutate(hub_in_rinker_2020_RAA_EAS = ifelse(CLASS=="RAA", TRUE, FALSE)) %>% 
	mutate(hub_in_rinker_2020_RHA_EAS = ifelse(CLASS=="RHA", TRUE, FALSE)) %>% 
	mutate(hub_in_rinker_2020_RA_EAS = (hub_in_rinker_2020_RAA_EAS | hub_in_rinker_2020_RHA_EAS))
Rinker_et_al_EAS <- Rinker_et_al_EAS %>% 
	dplyr::select(starts_with("hub_"), !starts_with("hub_")) %>% 
	dplyr::rename(EAS_INTROG_ALLELE=INTROG_ALLELE, EAS_INTROG_AF=INTROG_AF, EAS_INTROG_V16_LD_haplotype=V16_LD_haplotype) %>% 
	dplyr::select(-c("CHROM", "POS", "CLASS", "POP"))

# process Rinker et al SAS
Rinker_et_al_SAS <- Rinker_et_al_SAS %>% 
	mutate(
		hub_reference_genome = "hg19", 
		hub_variant_CHROM = gsub("chr", "", CHROM), 
		hub_variant_POS = POS
	)
Rinker_et_al_SAS <- Rinker_et_al_SAS %>% 
	mutate(hub_in_rinker_2020_NDA_SAS = ifelse(CLASS=="NDA", TRUE, FALSE)) %>% 
	mutate(hub_in_rinker_2020_RAA_SAS = ifelse(CLASS=="RAA", TRUE, FALSE)) %>% 
	mutate(hub_in_rinker_2020_RHA_SAS = ifelse(CLASS=="RHA", TRUE, FALSE)) %>% 
	mutate(hub_in_rinker_2020_RA_SAS = (hub_in_rinker_2020_RAA_SAS | hub_in_rinker_2020_RHA_SAS))
Rinker_et_al_SAS <- Rinker_et_al_SAS %>% 
	dplyr::select(starts_with("hub_"), !starts_with("hub_")) %>% 
	dplyr::rename(SAS_INTROG_ALLELE=INTROG_ALLELE, SAS_INTROG_AF=INTROG_AF, SAS_INTROG_V16_LD_haplotype=V16_LD_haplotype) %>% 
	dplyr::select(-c("CHROM", "POS", "CLASS", "POP"))

# join Rinker et al
Rinker_et_al_join <- full_join(Rinker_et_al_EUR, full_join(Rinker_et_al_EAS, Rinker_et_al_SAS))
Rinker_et_al_join <- Rinker_et_al_join %>% 
	mutate(hub_in_rinker_2020_NDA = (hub_in_rinker_2020_NDA_EUR|hub_in_rinker_2020_NDA_EAS|hub_in_rinker_2020_NDA_SAS)) %>% 
	mutate(hub_in_rinker_2020_RAA = (hub_in_rinker_2020_RAA_EUR|hub_in_rinker_2020_RAA_EAS|hub_in_rinker_2020_RAA_SAS)) %>% 
	mutate(hub_in_rinker_2020_RHA = (hub_in_rinker_2020_RHA_EUR|hub_in_rinker_2020_RHA_EAS|hub_in_rinker_2020_RHA_SAS)) %>% 
	mutate(hub_in_rinker_2020_RA = (hub_in_rinker_2020_RA_EUR|hub_in_rinker_2020_RA_EAS|hub_in_rinker_2020_RA_SAS))
Rinker_et_al_join <- Rinker_et_al_join %>% 
	dplyr::select(starts_with("hub_"), !starts_with("hub_"))

# save
write_tsv(Rinker_et_al_EUR, gzfile("../../results/preprocess_1KGP_SNPs/Rinker_et_al_EUR_hub.txt.gz"))
write_tsv(Rinker_et_al_EAS, gzfile("../../results/preprocess_1KGP_SNPs/Rinker_et_al_EAS_hub.txt.gz"))
write_tsv(Rinker_et_al_SAS, gzfile("../../results/preprocess_1KGP_SNPs/Rinker_et_al_SAS_hub.txt.gz"))
write_tsv(Rinker_et_al_join, gzfile("../../results/preprocess_1KGP_SNPs/Rinker_et_al_join_hub.txt.gz"))

# load 1KGP SNPs
ALL_1KGP_phase3_hub <- as_tibble(fread("../../results/preprocess_1KGP_SNPs/ALL_1KGP_phase3_hub.txt.gz"))

# load archaic SNPs
archaics_all_hub <- as_tibble(fread("../../results/preprocess_1KGP_SNPs/archaics_all_hub.txt.gz"))

# join 1KGP SNP info
Rinker_et_al_EUR_anno <- ALL_1KGP_phase3_hub %>% right_join(Rinker_et_al_EUR) %>% 
	dplyr::select(!contains("_INTROG_")) %>% 
	dplyr::select(starts_with("hub_"), !starts_with("hub_"))
Rinker_et_al_EAS_anno <- ALL_1KGP_phase3_hub %>% right_join(Rinker_et_al_EAS) %>% 
	dplyr::select(!contains("_INTROG_")) %>% 
	dplyr::select(starts_with("hub_"), !starts_with("hub_"))
Rinker_et_al_SAS_anno <- ALL_1KGP_phase3_hub %>% right_join(Rinker_et_al_SAS) %>% 
	dplyr::select(!contains("_INTROG_")) %>% 
	dplyr::select(starts_with("hub_"), !starts_with("hub_"))
Rinker_et_al_join_anno <- ALL_1KGP_phase3_hub %>% right_join(Rinker_et_al_join) %>% 
	dplyr::select(!contains("_INTROG_")) %>% 
	dplyr::select(starts_with("hub_"), !starts_with("hub_"))

# join archaic SNP info
Rinker_et_al_join_anno <- Rinker_et_al_join_anno %>% left_join(archaics_all_hub) %>% 
	dplyr::select(starts_with("hub_"), !starts_with("hub_"))

# save hub file
write_tsv(Rinker_et_al_join_anno, gzfile("../../results/preprocess_1KGP_SNPs/Rinker_et_al_join_B_stat_mask_1KGP_archaic_hub.txt.gz"))

# get BED file
Rinker_et_al_join_B_stat_mask_1KGP_archaic_hub <- as_tibble(fread("../../results/preprocess_1KGP_SNPs/Rinker_et_al_join_B_stat_mask_1KGP_archaic_hub.txt.gz"))
Rinker_et_al_join_B_stat_mask_1KGP_archaic_bed <- "../../results/preprocess_1KGP_SNPs/Rinker_et_al_join_B_stat_mask_1KGP_archaic_hub.bed"
hub_to_bed(Rinker_et_al_join_B_stat_mask_1KGP_archaic_hub, Rinker_et_al_join_B_stat_mask_1KGP_archaic_bed)
