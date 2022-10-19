#!/bin/R

library(tidyverse)
library(data.table)

source("../get_helper.R")

# load table
ALL_1KGP_phase3_hub <- as_tibble(fread("../../results/preprocess_1KGP_SNPs/ALL_1KGP_phase3_hub.txt.gz"))

# output bed
# hub_to_bed(ALL_1KGP_phase3_hub, "../../results/annotate_splice_prediction/ALL_1KGP_phase3_hub.bed")

for (chr in as.character(c(1:22, "X"))) {
	print(chr)
	ALL_1KGP_phase3_hub_temp <- ALL_1KGP_phase3_hub %>% filter(hub_variant_CHROM == chr)
	hub_to_bed(ALL_1KGP_phase3_hub_temp, paste("../../results/annotate_splice_prediction/ALL_1KGP_phase3_hub_chr", chr, ".bed", sep=""))
}
