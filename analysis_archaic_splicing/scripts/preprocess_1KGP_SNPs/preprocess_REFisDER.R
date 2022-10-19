#!/bin/sh

library(tidyverse)
library(data.table)
library(plyranges)
library(wrapr)
library(vcfR)
library(rtracklayer)
source("../get_helper.R")

# load
hg19_REFisDER_vcf <- read.vcfR("../../../final-splicing-archaic-variants/data/annotate_genome/fasta/hg19_REFisDER/hg19_REFisDER_hub.vcf.gz")

# drop INFO col
hg19_REFisDER_hub <- hg19_REFisDER_vcf@fix %>% as_tibble() %>% 
	dplyr::select(-INFO)

# # filter to SNPs, unnecessary
# print("filter to SNPs")
# table(hg19_REFisDER_hub$ALT)
# hg19_REFisDER_hub <- hg19_REFisDER_hub %>% 
# 	filter(
# 		REF %in% c("A", "C", "T", "G"),
# 		ALT %in% c("A", "C", "T", "G")
# 	)

# convert to GRange
print("convert to GRange")
hg19_REFisDER_hub <- hg19_REFisDER_hub %>% 
	mutate(POS = as.numeric(POS))

hg19_REFisDER_gr <- hg19_REFisDER_hub %>% 
	dplyr::rename(
		seqnames = CHROM
	) %>% 
	mutate(
		seqnames = paste("chr", seqnames, sep=""), 
		start = POS, 
		end = POS) %>% 
	GRanges()

# get B statistics
print("get B statistics")
B_stats_hg19_gr <- GRanges(as_tibble(fread("../../results/preprocess_1KGP_SNPs/B_stats_hg19_gr.txt.gz")))
# 	possible for liftOver to result in multiple overlapping regions
# 	in this case, just pick the first such overlapping region
archaic_overlap <- findOverlaps(hg19_REFisDER_gr, B_stats_hg19_gr) %>% 
	as_tibble() %>% filter(!duplicated(queryHits))
archaic_overlap_index <- archaic_overlap$queryHits
archaic_overlap_B_statistic <- B_stats_hg19_gr$B_statistic[archaic_overlap$subjectHits]
archaic_overlap_B_statistic_bin <- B_stats_hg19_gr$B_statistic_bin[archaic_overlap$subjectHits]
# 	now append B statistic information to table based on index
hg19_REFisDER_hub$B_statistic <- NA
hg19_REFisDER_hub$B_statistic_bin <- NA
hg19_REFisDER_hub$B_statistic[archaic_overlap_index] <- archaic_overlap_B_statistic
hg19_REFisDER_hub$B_statistic_bin[archaic_overlap_index] <- archaic_overlap_B_statistic_bin

# get mask overlap
print("get inter mask overlap")
archaic_mask_inter <- GRanges(as_tibble(fread("../../results/preprocess_1KGP_SNPs/archaic_mask_inter.bed.txt.gz")))
seqlevels(archaic_mask_inter) <- paste("chr", seqlevels(archaic_mask_inter), sep="")
mask_overlap <- findOverlaps(hg19_REFisDER_gr, archaic_mask_inter) %>% 
	as_tibble() %>% filter(!duplicated(queryHits))
mask_index <- mask_overlap$queryHits
hg19_REFisDER_hub$archaic_mask_inter <- FALSE 
hg19_REFisDER_hub$archaic_mask_inter[mask_index] <- TRUE

# get mask overlap
print("get union mask overlap")
archaic_mask_union <- GRanges(as_tibble(fread("../../results/preprocess_1KGP_SNPs/archaic_mask_union.bed.txt.gz")))
seqlevels(archaic_mask_union) <- paste("chr", seqlevels(archaic_mask_union), sep="")
mask_overlap <- findOverlaps(hg19_REFisDER_gr, archaic_mask_union) %>% 
	as_tibble() %>% filter(!duplicated(queryHits))
mask_index <- mask_overlap$queryHits
hg19_REFisDER_hub$archaic_mask_union <- FALSE 
hg19_REFisDER_hub$archaic_mask_union[mask_index] <- TRUE

# convert to hub
print("convert to hub")
hg19_REFisDER_hub <- hg19_REFisDER_hub %>% 
	mutate(
		hub_reference_genome = "hg19", 
		hub_variant_CHROM = CHROM, 
		hub_variant_POS = POS, 
		hub_variant_ID = NA, 
		hub_variant_REF = REF, 
		hub_variant_ALT = ALT
	)

# add ancestral
print("add ancestral")
hg19_REFisDER_hub <- hg19_REFisDER_hub %>% 
	mutate_hub_variant_ANC() %>% 
	mutate_hub_variant_DER() %>% 
	mutate_hub_variant_ID() %>% 
	split_ancestral() %>% 
	standard_sort()

# # temp, remove
# hg19_REFisDER_hub <- as_tibble(fread("../../results/preprocess_1KGP_SNPs/hg19_REFisDER_hub.txt.gz"))

# drop excess cols
print("drop excess cols")
hg19_REFisDER_hub <- hg19_REFisDER_hub %>% 
	dplyr::select(-c(CHROM, POS, ID, REF, ALT, QUAL, FILTER))

# save file
print("save file")
hg19_REFisDER_hub <- hg19_REFisDER_hub %>% 
	dplyr::select(starts_with("hub_"), !starts_with("hub_"))
write_tsv(hg19_REFisDER_hub, gzfile("../../results/preprocess_1KGP_SNPs/hg19_REFisDER_hub.txt.gz"))

# # temp, remove
# hg19_REFisDER_hub <- as_tibble(fread("../../results/preprocess_1KGP_SNPs/hg19_REFisDER_hub.txt.gz"))

# add 1KG and archaic info
print("add 1KG info")
ALL_1KGP_phase3_hub <- as_tibble(fread("../../results/preprocess_1KGP_SNPs/ALL_1KGP_phase3_hub.txt.gz"))
hg19_REFisDER_1KGP_hub <- hg19_REFisDER_hub %>% left_join(ALL_1KGP_phase3_hub)
write_tsv(hg19_REFisDER_1KGP_hub, gzfile("../../results/preprocess_1KGP_SNPs/hg19_REFisDER_1KGP_hub.txt.gz"))

# temp, remove
# hg19_REFisDER_1KGP_hub <- as_tibble(fread("../../results/preprocess_1KGP_SNPs/hg19_REFisDER_1KGP_hub.txt.gz"))

# add archaic info
print("add archaic info")
archaics_all_hub <- as_tibble(fread("../../results/preprocess_1KGP_SNPs/archaics_all_hub.txt.gz"))
hg19_REFisDER_1KGP_archaic_hub <- hg19_REFisDER_1KGP_hub %>% left_join(archaics_all_hub)
write_tsv(hg19_REFisDER_1KGP_archaic_hub, gzfile("../../results/preprocess_1KGP_SNPs/hg19_REFisDER_1KGP_archaic_hub.txt.gz"))

# gnomad sites
# 	temp, remove
hg19_REFisDER_1KGP_archaic_hub <- as_tibble(fread("../../results/preprocess_1KGP_SNPs/hg19_REFisDER_1KGP_archaic_hub.txt.gz"))

# 	table to bed
print("gnomAD table to bed")
hg19_REFisDER_1KGP_archaic_bed <- "../../results/preprocess_1KGP_SNPs/hg19_REFisDER_1KGP_archaic_hub.bed"
hub_to_bed(hg19_REFisDER_1KGP_archaic_hub, hg19_REFisDER_1KGP_archaic_bed)
