#!/bin/R

# Filter 1KGP SNP VEP info, filter to MAF >= 1%

library(tidyverse)
library(data.table)

# load 
ALL_1KGP_phase3_MAF0.01_hapR2_hub_vep <- as_tibble(fread("../../results/vep_annotations_SNPs/ALL_1KGP_phase3_MAF0.01_hapR2_hub_vep.txt.gz"))

# filter 1KGP variants
ALL_1KGP_phase3_MAF0.01_hapR2_hub_vep <- ALL_1KGP_phase3_MAF0.01_hapR2_hub_vep %>% 
	filter(hub_variant_REF %in% c("A", "C", "T", "G")) %>%  # SNP only
	filter(hub_variant_ALT %in% c("A", "C", "T", "G")) %>%  # SNP only
	filter(!(hub_variant_CHROM %in% c("X", "Y")))  # autosomes only

# save filter 1KGP variants
write_tsv(ALL_1KGP_phase3_MAF0.01_hapR2_hub_vep, "../../results/enrichment_GTEx_QTLs_v2/ALL_1KGP_phase3_MAF0.01_hapR2_hub_vep.txt.gz")
