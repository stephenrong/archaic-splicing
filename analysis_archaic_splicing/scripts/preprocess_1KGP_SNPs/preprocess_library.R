#!/bin/R

# Preprocess old version of variants under consideration (get B, get mask, add 1 KGP SNP info, add archaic SNP info, convert SNP to BED [add gnomAD SNP info])

library(tidyverse)
library(data.table)
library(plyranges)
library(wrapr)
library(vcfR)
library(rtracklayer)
source("../get_helper.R")

# load
merge_variants_B_stat_mask_hub <- as_tibble(fread("../../data/premapsy_variants/merge_variants_library_hub.txt.gz"))

# convert to GRange
print("convert to GRange")
merge_variants_B_stat_mask_hub <- merge_variants_B_stat_mask_hub %>% 
	mutate(hub_variant_POS = as.numeric(hub_variant_POS))

merge_variants_B_stat_mask_gr <- merge_variants_B_stat_mask_hub %>% 
	dplyr::rename(
		seqnames = hub_variant_CHROM
	) %>% 
	mutate(
		seqnames = paste("chr", seqnames, sep=""), 
		start = hub_variant_POS, 
		end = hub_variant_POS) %>% 
	GRanges()

# get B statistics
print("get B statistics")
B_stats_hg19_gr <- GRanges(as_tibble(fread("../../results/preprocess_1KGP_SNPs/B_stats_hg19_gr.txt.gz")))
# 	possible for liftOver to result in multiple overlapping regions
# 	in this case, just pick the first such overlapping region
archaic_overlap <- findOverlaps( merge_variants_B_stat_mask_gr, B_stats_hg19_gr) %>% 
	as_tibble() %>% filter(!duplicated(queryHits))
archaic_overlap_index <- archaic_overlap$queryHits
archaic_overlap_B_statistic <- B_stats_hg19_gr$B_statistic[archaic_overlap$subjectHits]
archaic_overlap_B_statistic_bin <- B_stats_hg19_gr$B_statistic_bin[archaic_overlap$subjectHits]
# 	now append B statistic information to table based on index
merge_variants_B_stat_mask_hub$B_statistic <- NA
merge_variants_B_stat_mask_hub$B_statistic_bin <- NA
merge_variants_B_stat_mask_hub$B_statistic[archaic_overlap_index] <- archaic_overlap_B_statistic
merge_variants_B_stat_mask_hub$B_statistic_bin[archaic_overlap_index] <- archaic_overlap_B_statistic_bin

# get mask overlap
print("get inter mask overlap")
archaic_mask_inter <- GRanges(as_tibble(fread("../../results/preprocess_1KGP_SNPs/archaic_mask_inter.bed.txt.gz")))
seqlevels(archaic_mask_inter) <- paste("chr", seqlevels(archaic_mask_inter), sep="")
mask_overlap <- findOverlaps( merge_variants_B_stat_mask_gr, archaic_mask_inter) %>% 
	as_tibble() %>% filter(!duplicated(queryHits))
mask_index <- mask_overlap$queryHits
merge_variants_B_stat_mask_hub$archaic_mask_inter <- FALSE 
merge_variants_B_stat_mask_hub$archaic_mask_inter[mask_index] <- TRUE

# get mask overlap
print("get union mask overlap")
archaic_mask_union <- GRanges(as_tibble(fread("../../results/preprocess_1KGP_SNPs/archaic_mask_union.bed.txt.gz")))
seqlevels(archaic_mask_union) <- paste("chr", seqlevels(archaic_mask_union), sep="")
mask_overlap <- findOverlaps( merge_variants_B_stat_mask_gr, archaic_mask_union) %>% 
	as_tibble() %>% filter(!duplicated(queryHits))
mask_index <- mask_overlap$queryHits
merge_variants_B_stat_mask_hub$archaic_mask_union <- FALSE 
merge_variants_B_stat_mask_hub$archaic_mask_union[mask_index] <- TRUE

# save file
print("save file")
merge_variants_B_stat_mask_hub <- merge_variants_B_stat_mask_hub %>% 
	dplyr::select(starts_with("hub_"), !starts_with("hub_"))
write_tsv(merge_variants_B_stat_mask_hub, gzfile("../../results/preprocess_1KGP_SNPs/merge_variants_B_stat_mask_hub.txt.gz"))

# # temp, remove
# merge_variants_B_stat_mask_hub <- as_tibble(fread("../../results/preprocess_1KGP_SNPs/merge_variants_B_stat_mask_hub.txt.gz"))

# add 1KG and archaic info
print("add 1KG info")
ALL_1KGP_phase3_hub <- as_tibble(fread("../../results/preprocess_1KGP_SNPs/ALL_1KGP_phase3_hub.txt.gz"))
merge_variants_B_stat_mask_1KGP_hub <- merge_variants_B_stat_mask_hub %>% left_join(ALL_1KGP_phase3_hub)
write_tsv(merge_variants_B_stat_mask_1KGP_hub, gzfile("../../results/preprocess_1KGP_SNPs/merge_variants_B_stat_mask_1KGP_hub.txt.gz"))

# # temp, remove
# merge_variants_B_stat_mask_1KGP_hub <- as_tibble(fread("../../results/preprocess_1KGP_SNPs/merge_variants_B_stat_mask_1KGP_hub.txt.gz"))

# add archaic info
print("add archaic info")
archaics_all_hub <- as_tibble(fread("../../results/preprocess_1KGP_SNPs/archaics_all_hub.txt.gz"))
merge_variants_B_stat_mask_1KGP_archaic_hub <- merge_variants_B_stat_mask_1KGP_hub %>% left_join(archaics_all_hub)
write_tsv(merge_variants_B_stat_mask_1KGP_archaic_hub, gzfile("../../results/preprocess_1KGP_SNPs/merge_variants_B_stat_mask_1KGP_archaic_hub.txt.gz"))

# gnomad sites
# 	temp, remove
merge_variants_B_stat_mask_1KGP_archaic_hub <- as_tibble(fread("../../results/preprocess_1KGP_SNPs/merge_variants_B_stat_mask_1KGP_archaic_hub.txt.gz"))

# 	table to bed
print("gnomAD table to bed")
merge_variants_B_stat_mask_1KGP_archaic_bed <- "../../results/preprocess_1KGP_SNPs/merge_variants_B_stat_mask_1KGP_archaic_hub.bed"
hub_to_bed(merge_variants_B_stat_mask_1KGP_archaic_hub, merge_variants_B_stat_mask_1KGP_archaic_bed)
