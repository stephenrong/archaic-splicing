#!/bin/R

# Preprocess 1KGP SNP data (convert to hub format, get AC, get B, save all, save AF >= 0.9)

library(tidyverse)
library(data.table)
library(plyranges)
library(wrapr)
library(vcfR)
library(rtracklayer)
source("../get_helper.R")

# load
ALL_1KGP_phase3_vcf <- read.vcfR("../../../final-splicing-archaic-variants/data/annotate_popgen_1000G/1000Genomes_phase3/1000Genomes_phase3_filterAF/ALL.phase3_shapeit2_mvncall_integrated_v1b.20130502.genotypes.dropgeno.concat.norm.vcf.gz")

# clean up AFs
print("clean up AFs")
col_list <- c("AC", "AF", "AN", "NS", "DP", "EAS_AF", "AMR_AF", "AFR_AF", "EUR_AF", "SAS_AF")
ALL_1KGP_phase3_hub <- ALL_1KGP_phase3_vcf@fix %>% as_tibble() %>% 
	filter(
		REF %in% c("A", "C", "T", "G"),
		ALT %in% c("A", "C", "T", "G")
	) %>% 
	mutate(INFO = gsub(";AA=.*|;VT=.*", "", INFO)) %>% 
	separate(col=INFO, into=col_list, sep=";")
for (col in col_list) {
	let(c(VAR=col), 
		ALL_1KGP_phase3_hub <- ALL_1KGP_phase3_hub %>% 
			mutate(VAR = as.numeric(gsub(".*=", "", VAR))))
}

# recover ACs
print("recover ACs")
EAS_AN <- 1008
AMR_AN <- 694
AFR_AN <- 1322
EUR_AN <- 1006
SAS_AN <- 978
ALL_1KGP_phase3_hub <- ALL_1KGP_phase3_hub %>% 
	mutate(
		# EAS_AN = EAS_AN, 
		EAS_AC = round(EAS_AN*EAS_AF),
		# AMR_AN = AMR_AN, 
		AMR_AC = round(AMR_AN*AMR_AF), 
		# AFR_AN = AFR_AN, 
		AFR_AC = round(AFR_AN*AFR_AF), 
		# EUR_AN = EUR_AN, 
		EUR_AC = round(EUR_AN*EUR_AF), 
		# SAS_AN = SAS_AN, 
		SAS_AC = round(SAS_AN*SAS_AF)
	)

# convert to GRange
print("convert to GRange")
ALL_1KGP_phase3_hub <- ALL_1KGP_phase3_hub %>% 
	mutate(POS = as.numeric(POS))

ALL_1KGP_phase3_gr <- ALL_1KGP_phase3_hub %>% 
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
ALL_1KGP_phase3_overlap <- findOverlaps(ALL_1KGP_phase3_gr, B_stats_hg19_gr) %>% 
	as_tibble() %>% filter(!duplicated(queryHits))
ALL_1KGP_phase3_overlap_index <- ALL_1KGP_phase3_overlap$queryHits
ALL_1KGP_phase3_overlap_B_statistic <- B_stats_hg19_gr$B_statistic[ALL_1KGP_phase3_overlap$subjectHits]
ALL_1KGP_phase3_overlap_B_statistic_bin <- B_stats_hg19_gr$B_statistic_bin[ALL_1KGP_phase3_overlap$subjectHits]
# 	now append B statistic information to table based on index
ALL_1KGP_phase3_hub$B_statistic <- NA
ALL_1KGP_phase3_hub$B_statistic_bin <- NA
ALL_1KGP_phase3_hub$B_statistic[ALL_1KGP_phase3_overlap_index] <- ALL_1KGP_phase3_overlap_B_statistic
ALL_1KGP_phase3_hub$B_statistic_bin[ALL_1KGP_phase3_overlap_index] <- ALL_1KGP_phase3_overlap_B_statistic_bin

# get mask overlap
print("get inter mask overlap")
archaic_mask_inter <- GRanges(as_tibble(fread("../../results/preprocess_1KGP_SNPs/archaic_mask_inter.bed.txt.gz")))
seqlevels(archaic_mask_inter) <- paste("chr", seqlevels(archaic_mask_inter), sep="")
mask_overlap <- findOverlaps(ALL_1KGP_phase3_gr, archaic_mask_inter) %>% 
	as_tibble() %>% filter(!duplicated(queryHits))
mask_index <- mask_overlap$queryHits
ALL_1KGP_phase3_hub$archaic_mask_inter <- FALSE 
ALL_1KGP_phase3_hub$archaic_mask_inter[mask_index] <- TRUE

# get mask overlap
print("get union mask overlap")
archaic_mask_union <- GRanges(as_tibble(fread("../../results/preprocess_1KGP_SNPs/archaic_mask_union.bed.txt.gz")))
seqlevels(archaic_mask_union) <- paste("chr", seqlevels(archaic_mask_union), sep="")
mask_overlap <- findOverlaps(ALL_1KGP_phase3_gr, archaic_mask_union) %>% 
	as_tibble() %>% filter(!duplicated(queryHits))
mask_index <- mask_overlap$queryHits
ALL_1KGP_phase3_hub$archaic_mask_union <- FALSE 
ALL_1KGP_phase3_hub$archaic_mask_union[mask_index] <- TRUE

# convert to hub format
print("convert to hub format")
ALL_1KGP_phase3_hub <- ALL_1KGP_phase3_hub %>% 
	mutate(
		hub_reference_genome = "hg19", 
		hub_variant_CHROM = CHROM, 
		hub_variant_POS = POS, 
		hub_variant_ID = NA, 
		hub_variant_REF = REF, 
		hub_variant_ALT = ALT
	)

print("get ancestral allele")
ALL_1KGP_phase3_hub <- ALL_1KGP_phase3_hub %>% 
	mutate_hub_variant_ANC() %>% 
	mutate_hub_variant_DER() %>% 
	mutate_hub_variant_ID() %>% 
	split_ancestral() %>% 
	standard_sort()  # would usually use clean_up_hub 
					 #   but is too memory intensive

# remove AC == 0
ALL_1KGP_phase3_hub <- ALL_1KGP_phase3_hub %>% 
	filter(!(AC == 0))  # variants present in 1KGP

# remove duplicates
ALL_1KGP_phase3_hub <- ALL_1KGP_phase3_hub %>% 
	mutate(unique = !duplicated(hub_variant_ID)) %>% 
	filter(unique) %>% dplyr::select(-unique)  # second is usually multi-allelic

# save file
print("save file")
ALL_1KGP_phase3_hub <- ALL_1KGP_phase3_hub %>% 
	dplyr::select(starts_with("hub_"), !starts_with("hub_")) # %>% mutate(index = row_number())
write_tsv(ALL_1KGP_phase3_hub, gzfile("../../results/preprocess_1KGP_SNPs/ALL_1KGP_phase3_hub.txt.gz"))

# create filters
print("create filters")
ALL_1KGP_phase3_hub <- as_tibble(fread("../../results/preprocess_1KGP_SNPs/ALL_1KGP_phase3_hub.txt.gz"))

ALL_1KGP_phase3_greater0.9_hub <- ALL_1KGP_phase3_hub %>% filter(AF >= 0.9)
write_tsv(ALL_1KGP_phase3_greater0.9_hub, gzfile("../../results/preprocess_1KGP_SNPs/ALL_1KGP_phase3_greater0.9_hub.txt.gz"))

ALL_1KGP_phase3_greater0.5_hub <- ALL_1KGP_phase3_hub %>% filter(AF >= 0.5)
write_tsv(ALL_1KGP_phase3_greater0.5_hub, gzfile("../../results/preprocess_1KGP_SNPs/ALL_1KGP_phase3_greater0.5_hub.txt.gz"))

ALL_1KGP_phase3_greater0.1_hub <- ALL_1KGP_phase3_hub %>% filter(AF >= 0.1)
write_tsv(ALL_1KGP_phase3_greater0.1_hub, gzfile("../../results/preprocess_1KGP_SNPs/ALL_1KGP_phase3_greater0.1_hub.txt.gz"))

ALL_1KGP_phase3_MAF0.01_hub <- ALL_1KGP_phase3_hub %>% filter(EUR_AF >= 0.01, EUR_AF <= 0.99)
write_tsv(ALL_1KGP_phase3_MAF0.01_hub, gzfile("../../results/preprocess_1KGP_SNPs/ALL_1KGP_phase3_MAF0.01_hub.txt.gz"))

# table to BED
print("table to BED")
ALL_1KGP_phase3_greater0.9_hub <- as_tibble(fread("../../results/preprocess_1KGP_SNPs/ALL_1KGP_phase3_greater0.9_hub.txt.gz"))
ALL_1KGP_phase3_greater0.9_bed <- "../../results/preprocess_1KGP_SNPs/ALL_1KGP_phase3_greater0.9_hub.bed"
hub_to_bed(ALL_1KGP_phase3_greater0.9_hub, ALL_1KGP_phase3_greater0.9_bed)

ALL_1KGP_phase3_greater0.5_hub <- as_tibble(fread("../../results/preprocess_1KGP_SNPs/ALL_1KGP_phase3_greater0.5_hub.txt.gz"))
ALL_1KGP_phase3_greater0.5_bed <- "../../results/preprocess_1KGP_SNPs/ALL_1KGP_phase3_greater0.5_hub.bed"
hub_to_bed(ALL_1KGP_phase3_greater0.5_hub, ALL_1KGP_phase3_greater0.5_bed)

ALL_1KGP_phase3_greater0.1_hub <- as_tibble(fread("../../results/preprocess_1KGP_SNPs/ALL_1KGP_phase3_greater0.1_hub.txt.gz"))
ALL_1KGP_phase3_greater0.1_bed <- "../../results/preprocess_1KGP_SNPs/ALL_1KGP_phase3_greater0.1_hub.bed"
hub_to_bed(ALL_1KGP_phase3_greater0.1_hub, ALL_1KGP_phase3_greater0.1_bed)

ALL_1KGP_phase3_MAF0.01_hub <- as_tibble(fread("../../results/preprocess_1KGP_SNPs/ALL_1KGP_phase3_MAF0.01_hub.txt.gz"))
ALL_1KGP_phase3_MAF0.01_bed <- "../../results/preprocess_1KGP_SNPs/ALL_1KGP_phase3_MAF0.01_hub.bed"
hub_to_bed(ALL_1KGP_phase3_MAF0.01_hub, ALL_1KGP_phase3_MAF0.01_bed)
