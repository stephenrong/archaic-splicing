#!/bin/R

library(tidyverse)
library(data.table)
# source("../get_helper.R")

# load files
print("load files")
ALL_1KGP_phase3_MAF0.01_hapR2_hub <- as_tibble(fread("../../results/preprocess_1KGP_SNPs/ALL_1KGP_phase3_MAF0.01_hapR2_hub.txt.gz"))
final_v2_variants_B_stat_mask_1KGP_archaic_gnomAD_hapR2_hub <- as_tibble(fread("../../results/preprocess_1KGP_SNPs/final_v2_variants_B_stat_mask_1KGP_archaic_gnomAD_hapR2_hub.txt.gz"))

# filter 1KGP variants
ALL_1KGP_phase3_MAF0.01_hapR2_hub <- ALL_1KGP_phase3_MAF0.01_hapR2_hub %>% 
	filter(hub_variant_REF %in% c("A", "C", "T", "G")) %>%  # SNP only
	filter(hub_variant_ALT %in% c("A", "C", "T", "G")) %>%  # SNP only
	filter(!(hub_variant_CHROM %in% c("X", "Y")))  # autosomes only
write_tsv(ALL_1KGP_phase3_MAF0.01_hapR2_hub, gzfile("../../results/enrichment_GTEx_QTLs_v2/ALL_1KGP_phase3_MAF0.01_hapR2_filtered_hub.txt.gz"))

# # save in VCF format
# hub_to_vcf(ALL_1KGP_phase3_MAF0.01_hapR2_hub, "ALL_1KGP_phase3_MAF0.01_hapR2_hub", "../../results/enrichment_GTEx_QTLs_v2/ALL_1KGP_phase3_MAF0.01_hapR2_filtered_hub.vcf", BSgenome.Hsapiens.UCSC.hg19)

# # save in BED format
# hub_to_bed(ALL_1KGP_phase3_MAF0.01_hapR2_hub, "../../results/enrichment_GTEx_QTLs_v2/ALL_1KGP_phase3_MAF0.01_hapR2_filtered_hub.bed")

# get not archaically introgressed
print("get nonAI")
ALL_1KGP_phase3_MAF0.01_hapR2_hub <- ALL_1KGP_phase3_MAF0.01_hapR2_hub %>% 
	filter(!(hub_variant_ID %in% filter(final_v2_variants_B_stat_mask_1KGP_archaic_gnomAD_hapR2_hub, hub_in_final_study_introgressed)$hub_variant_ID))
  
# add index
ALL_1KGP_phase3_MAF0.01_hapR2_hub <- ALL_1KGP_phase3_MAF0.01_hapR2_hub %>% 
	mutate(index = row_number())
write_tsv(ALL_1KGP_phase3_MAF0.01_hapR2_hub, gzfile("../../results/enrichment_GTEx_QTLs_v2/ALL_1KGP_phase3_MAF0.01_hapR2_notAI_hub.txt.gz"))

# count bins
print("count bins")
ALL_1KGP_phase3_MAF0.01_hapR2_hub_EUR_AC_bin <- ALL_1KGP_phase3_MAF0.01_hapR2_hub %>% 
	group_by(EUR_AC_bin) %>% tally()
ALL_1KGP_phase3_MAF0.01_hapR2_hub_EUR_hapR2tag_bin <- ALL_1KGP_phase3_MAF0.01_hapR2_hub %>% 
	group_by(EUR_hapR2tag_bin) %>% tally()

print("save bins")
write_tsv(ALL_1KGP_phase3_MAF0.01_hapR2_hub_EUR_AC_bin, "../../results/enrichment_GTEx_QTLs_v2/ALL_1KGP_phase3_MAF0.01_hapR2_notAI_hub_EUR_AC_bin.txt")
write_tsv(ALL_1KGP_phase3_MAF0.01_hapR2_hub_EUR_hapR2tag_bin, "../../results/enrichment_GTEx_QTLs_v2/ALL_1KGP_phase3_MAF0.01_hapR2_notAI_hub_EUR_hapR2tag_bin.txt")

# create an index
# 	this will be used 
print("create index")
ALL_1KGP_phase3_MAF0.01_hapR2_hub_EUR_AC_bin_index <- ALL_1KGP_phase3_MAF0.01_hapR2_hub %>% 
	group_by(EUR_AC_bin)  %>% summarise(index = list(index))
ALL_1KGP_phase3_MAF0.01_hapR2_hub_EUR_hapR2tag_bin_index <- ALL_1KGP_phase3_MAF0.01_hapR2_hub %>% 
	group_by(EUR_hapR2tag_bin)  %>% summarise(index = list(index))

print("save index")
saveRDS(ALL_1KGP_phase3_MAF0.01_hapR2_hub_EUR_AC_bin_index, file="../../results/enrichment_GTEx_QTLs_v2/ALL_1KGP_phase3_MAF0.01_hapR2_notAI_hub_EUR_AC_bin_index.rds")
saveRDS(ALL_1KGP_phase3_MAF0.01_hapR2_hub_EUR_hapR2tag_bin_index, file="../../results/enrichment_GTEx_QTLs_v2/ALL_1KGP_phase3_MAF0.01_hapR2_notAI_hub_EUR_hapR2tag_bin_index.rds")

# count cross bins
# 	this is to check that discretization isn't too fine
print("count cross bins")
ALL_1KGP_phase3_MAF0.01_hapR2_hub_EUR_AC_binxEUR_hapR2tag_bin <- ALL_1KGP_phase3_MAF0.01_hapR2_hub %>% 
	group_by(EUR_AC_bin, EUR_hapR2tag_bin) %>% tally()

print("save bins")
write_tsv(ALL_1KGP_phase3_MAF0.01_hapR2_hub_EUR_AC_binxEUR_hapR2tag_bin, "../../results/enrichment_GTEx_QTLs_v2/ALL_1KGP_phase3_MAF0.01_hapR2_notAI_hub_EUR_AC_binxEUR_hapR2tag_bin.txt")

# create an index
# 	this will be used 
print("create index")
ALL_1KGP_phase3_MAF0.01_hapR2_hub_EUR_AC_binxEUR_hapR2tag_index <- ALL_1KGP_phase3_MAF0.01_hapR2_hub %>% 
	group_by(EUR_AC_bin, EUR_hapR2tag_bin) %>% summarise(index = list(index))

print("save index")
saveRDS(ALL_1KGP_phase3_MAF0.01_hapR2_hub_EUR_AC_binxEUR_hapR2tag_index, file="../../results/enrichment_GTEx_QTLs_v2/ALL_1KGP_phase3_MAF0.01_hapR2_notAI_hub_EUR_AC_binxEUR_hapR2tag_index.rds")
