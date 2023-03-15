#!/bin/R

# Take collated QTL files, add 1KGP SNP info

library(tidyverse)
library(data.table)

# load 1KGP variants
ALL_1KGP_phase3_MAF0.01_hapR2_hub <- as_tibble(fread("../../results/preprocess_1KGP_SNPs/ALL_1KGP_phase3_MAF0.01_hapR2_hub.txt.gz"))

# load GTEx files
collate_GTEx_aQTLs_hg19 <- as_tibble(fread("../../results/preprocess_GTEx_QTLs/collate_GTEx_aQTLs_hg19.txt.gz")) %>% 
	mutate(hub_variant_CHROM = as.character(hub_variant_CHROM))
collate_GTEx_sQTLs_hg19 <- as_tibble(fread("../../results/preprocess_GTEx_QTLs/collate_GTEx_sQTLs_hg19.txt.gz")) %>% 
	mutate(hub_variant_CHROM = as.character(hub_variant_CHROM))
collate_GTEx_eQTLs_hg19 <- as_tibble(fread("../../results/preprocess_GTEx_QTLs/collate_GTEx_eQTLs_hg19.txt.gz")) %>% 
	mutate(hub_variant_CHROM = as.character(hub_variant_CHROM))

# join 1KGP variants
join_in1KGP_GTEx_aQTLs_hg19 <- collate_GTEx_aQTLs_hg19 %>% left_join(ALL_1KGP_phase3_MAF0.01_hapR2_hub)
join_in1KGP_GTEx_sQTLs_hg19 <- collate_GTEx_sQTLs_hg19 %>% left_join(ALL_1KGP_phase3_MAF0.01_hapR2_hub)
join_in1KGP_GTEx_eQTLs_hg19 <- collate_GTEx_eQTLs_hg19 %>% left_join(ALL_1KGP_phase3_MAF0.01_hapR2_hub)

# filter 1KGP variants
filter_in1KGP_GTEx_aQTLs_hg19 <- join_in1KGP_GTEx_aQTLs_hg19 %>% 
	filter(hub_variant_REF %in% c("A", "C", "T", "G")) %>%  # SNP only
	filter(hub_variant_ALT %in% c("A", "C", "T", "G")) %>%  # SNP only
	filter(!(hub_variant_CHROM %in% c("X", "Y")))  # autosomes only
filter_in1KGP_GTEx_sQTLs_hg19 <- join_in1KGP_GTEx_sQTLs_hg19 %>% 
	filter(hub_variant_REF %in% c("A", "C", "T", "G")) %>%  # SNP only
	filter(hub_variant_ALT %in% c("A", "C", "T", "G")) %>%  # SNP only
	filter(!(hub_variant_CHROM %in% c("X", "Y")))  # autosomes only
filter_in1KGP_GTEx_eQTLs_hg19 <- join_in1KGP_GTEx_eQTLs_hg19 %>% 
	filter(hub_variant_REF %in% c("A", "C", "T", "G")) %>%  # SNP only
	filter(hub_variant_ALT %in% c("A", "C", "T", "G")) %>%  # SNP only
	filter(!(hub_variant_CHROM %in% c("X", "Y")))  # autosomes only

# # filter by European MAF
# filter_in1KGP_GTEx_aQTLs_hg19 <- filter_in1KGP_GTEx_aQTLs_hg19 %>% 
# 	filter(EUR_AF >= 0.01, EUR_AF <= 0.99)
# filter_in1KGP_GTEx_sQTLs_hg19 <- filter_in1KGP_GTEx_sQTLs_hg19 %>% 
# 	filter(EUR_AF >= 0.01, EUR_AF <= 0.99)
# filter_in1KGP_GTEx_eQTLs_hg19 <- filter_in1KGP_GTEx_eQTLs_hg19 %>% 
# 	filter(EUR_AF >= 0.01, EUR_AF <= 0.99)

# note 02-07-2023: 
# not necessary, since switched 
# from ALL_1KGP_phase3_hapR2_hub to
# ALL_1KGP_phase3_MAF0.01_hapR2_hub

# save filter 1KGP variants
write_tsv(filter_in1KGP_GTEx_aQTLs_hg19, "../../results/enrichment_GTEx_QTLs_v2/filter_in1KGP_GTEx_aQTLs_hg19.txt.gz")
write_tsv(filter_in1KGP_GTEx_sQTLs_hg19, "../../results/enrichment_GTEx_QTLs_v2/filter_in1KGP_GTEx_sQTLs_hg19.txt.gz")
write_tsv(filter_in1KGP_GTEx_eQTLs_hg19, "../../results/enrichment_GTEx_QTLs_v2/filter_in1KGP_GTEx_eQTLs_hg19.txt.gz")
