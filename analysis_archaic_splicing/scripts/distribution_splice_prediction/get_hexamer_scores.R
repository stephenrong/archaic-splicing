#!/bin/R

library(tidyverse)
library(data.table)
library(BSgenome.Hsapiens.UCSC.hg19)
library(BSgenome.Hsapiens.UCSC.hg38)
hg19 <- BSgenome.Hsapiens.UCSC.hg19
hg38 <- BSgenome.Hsapiens.UCSC.hg38

# load table of hexamer splicing scores
splicing_scores <- as_tibble(fread("../../data/annotate_hexamer_scores/chasin_rosenberg_scores.txt"))
splicing_scores_EI <- splicing_scores$Hexamer_EI_Ke_et_al_2011
splicing_scores_exonicA3SS <- splicing_scores$A3SS_Exonic_Rosenberg_et_al_2015
splicing_scores_exonicA5SS <- splicing_scores$A5SS_Exonic_Rosenberg_et_al_2015
# splicing_scores_intronicA3SS <- splicing_scores$A3SS_Intronic_Rosenberg_et_al_2015
# splicing_scores_intronicA5SS <- splicing_scores$A5SS_Intronic_Rosenberg_et_al_2015

# create lookup table of hexamer scores
names(splicing_scores_EI) <- splicing_scores$Hexamer
names(splicing_scores_exonicA3SS) <- splicing_scores$Hexamer
names(splicing_scores_exonicA5SS) <- splicing_scores$Hexamer
# names(splicing_scores_intronicA3SS) <- splicing_scores$Hexamer
# names(splicing_scores_intronicA5SS) <- splicing_scores$Hexamer

# create function that will
# get elevenmer sequence and compute delta
get_hexamer_scores <- function(input, ANCDER=FALSE) {
	# which reference genome
	reference <- c("hg19" = hg19, "hg38" = hg38)[[input$ss_reference_genome[[1]]]]

	# get wt sequence centered around variant
	temp <- input %>% mutate(ss_wt_elevenmer = 
		toupper(as.character(getSeq(reference, 
		names=paste("chr", ss_variant_CHROM, sep=""), 
		start=(ss_variant_POS-5), end=(ss_variant_POS+5), strand=ss_strand))))
	# get mt sequence centered around variant
	temp <- temp %>% 
		mutate(ss_mt_elevenmer = 
		ifelse(ss_strand == "+", 
			paste(substr(ss_wt_elevenmer, 1, 5), 
				ss_variant_ALT, substr(ss_wt_elevenmer, 7, 11), sep=""),
			paste(substr(ss_wt_elevenmer, 1, 5), 
				chartr("ACGTacgt", "TGCAtgca", ss_variant_ALT), 
				substr(ss_wt_elevenmer, 7, 11), sep=""))) # %>% ungroup()

	# break up wt and mt seqs into hexamers, take the mean, and the mt-wt difference
	# 	do for each for the five splicing scores (Chasin EI, 
	# 	Rosenberg exonic A3SS, A5SS, and Rosenberg intronic A3SS, A5SS)
	temp <- temp %>% rowwise() %>% 
		mutate(ss_EI_mt = 
			(splicing_scores_EI[[substr(ss_mt_elevenmer, 1, 6)]] + 
			splicing_scores_EI[[substr(ss_mt_elevenmer, 2, 7)]] + 
			splicing_scores_EI[[substr(ss_mt_elevenmer, 3, 8)]] + 
			splicing_scores_EI[[substr(ss_mt_elevenmer, 4, 9)]] + 
			splicing_scores_EI[[substr(ss_mt_elevenmer, 5, 10)]] + 
			splicing_scores_EI[[substr(ss_mt_elevenmer, 6, 11)]])/6) %>% 
		mutate(ss_EI_wt = 
			(splicing_scores_EI[[substr(ss_wt_elevenmer, 1, 6)]] + 
			splicing_scores_EI[[substr(ss_wt_elevenmer, 2, 7)]] + 
			splicing_scores_EI[[substr(ss_wt_elevenmer, 3, 8)]] + 
			splicing_scores_EI[[substr(ss_wt_elevenmer, 4, 9)]] + 
			splicing_scores_EI[[substr(ss_wt_elevenmer, 5, 10)]] + 
			splicing_scores_EI[[substr(ss_wt_elevenmer, 6, 11)]])/6) %>% 
		mutate(ss_EI = ss_EI_mt-ss_EI_wt) %>% 
		mutate(ss_exonicA3SS_mt = 
			(splicing_scores_exonicA3SS[[substr(ss_mt_elevenmer, 1, 6)]] + 
			splicing_scores_exonicA3SS[[substr(ss_mt_elevenmer, 2, 7)]] + 
			splicing_scores_exonicA3SS[[substr(ss_mt_elevenmer, 3, 8)]] + 
			splicing_scores_exonicA3SS[[substr(ss_mt_elevenmer, 4, 9)]] + 
			splicing_scores_exonicA3SS[[substr(ss_mt_elevenmer, 5, 10)]] + 
			splicing_scores_exonicA3SS[[substr(ss_mt_elevenmer, 6, 11)]])/6) %>% 
		mutate(ss_exonicA3SS_wt = 
			(splicing_scores_exonicA3SS[[substr(ss_wt_elevenmer, 1, 6)]] + 
			splicing_scores_exonicA3SS[[substr(ss_wt_elevenmer, 2, 7)]] + 
			splicing_scores_exonicA3SS[[substr(ss_wt_elevenmer, 3, 8)]] + 
			splicing_scores_exonicA3SS[[substr(ss_wt_elevenmer, 4, 9)]] + 
			splicing_scores_exonicA3SS[[substr(ss_wt_elevenmer, 5, 10)]] + 
			splicing_scores_exonicA3SS[[substr(ss_wt_elevenmer, 6, 11)]])/6) %>% 
		mutate(ss_exonicA3SS = ss_exonicA3SS_mt-ss_exonicA3SS_wt) %>% 
		mutate(ss_exonicA5SS_mt = 
			(splicing_scores_exonicA5SS[[substr(ss_mt_elevenmer, 1, 6)]] + 
			splicing_scores_exonicA5SS[[substr(ss_mt_elevenmer, 2, 7)]] + 
			splicing_scores_exonicA5SS[[substr(ss_mt_elevenmer, 3, 8)]] + 
			splicing_scores_exonicA5SS[[substr(ss_mt_elevenmer, 4, 9)]] + 
			splicing_scores_exonicA5SS[[substr(ss_mt_elevenmer, 5, 10)]] + 
			splicing_scores_exonicA5SS[[substr(ss_mt_elevenmer, 6, 11)]])/6) %>% 
		mutate(ss_exonicA5SS_wt = 
			(splicing_scores_exonicA5SS[[substr(ss_wt_elevenmer, 1, 6)]] + 
			splicing_scores_exonicA5SS[[substr(ss_wt_elevenmer, 2, 7)]] + 
			splicing_scores_exonicA5SS[[substr(ss_wt_elevenmer, 3, 8)]] + 
			splicing_scores_exonicA5SS[[substr(ss_wt_elevenmer, 4, 9)]] + 
			splicing_scores_exonicA5SS[[substr(ss_wt_elevenmer, 5, 10)]] + 
			splicing_scores_exonicA5SS[[substr(ss_wt_elevenmer, 6, 11)]])/6) # %>% 
		# mutate(ss_exonicA5SS = ss_exonicA5SS_mt-ss_exonicA5SS_wt) %>% 
		# mutate(ss_intronicA3SS_mt = 
		# 	(splicing_scores_intronicA3SS[[substr(ss_mt_elevenmer, 1, 6)]] + 
		# 	splicing_scores_intronicA3SS[[substr(ss_mt_elevenmer, 2, 7)]] + 
		# 	splicing_scores_intronicA3SS[[substr(ss_mt_elevenmer, 3, 8)]] + 
		# 	splicing_scores_intronicA3SS[[substr(ss_mt_elevenmer, 4, 9)]] + 
		# 	splicing_scores_intronicA3SS[[substr(ss_mt_elevenmer, 5, 10)]] + 
		# 	splicing_scores_intronicA3SS[[substr(ss_mt_elevenmer, 6, 11)]])/6) %>% 
		# mutate(ss_intronicA3SS_wt = 
		# 	(splicing_scores_intronicA3SS[[substr(ss_wt_elevenmer, 1, 6)]] + 
		# 	splicing_scores_intronicA3SS[[substr(ss_wt_elevenmer, 2, 7)]] + 
		# 	splicing_scores_intronicA3SS[[substr(ss_wt_elevenmer, 3, 8)]] + 
		# 	splicing_scores_intronicA3SS[[substr(ss_wt_elevenmer, 4, 9)]] + 
		# 	splicing_scores_intronicA3SS[[substr(ss_wt_elevenmer, 5, 10)]] + 
		# 	splicing_scores_intronicA3SS[[substr(ss_wt_elevenmer, 6, 11)]])/6) %>% 
		# mutate(ss_intronicA3SS = ss_intronicA3SS_mt-ss_intronicA3SS_wt) %>% 
		# mutate(ss_intronicA5SS_mt = 
		# 	(splicing_scores_intronicA5SS[[substr(ss_mt_elevenmer, 1, 6)]] + 
		# 	splicing_scores_intronicA5SS[[substr(ss_mt_elevenmer, 2, 7)]] + 
		# 	splicing_scores_intronicA5SS[[substr(ss_mt_elevenmer, 3, 8)]] + 
		# 	splicing_scores_intronicA5SS[[substr(ss_mt_elevenmer, 4, 9)]] + 
		# 	splicing_scores_intronicA5SS[[substr(ss_mt_elevenmer, 5, 10)]] + 
		# 	splicing_scores_intronicA5SS[[substr(ss_mt_elevenmer, 6, 11)]])/6) %>% 
		# mutate(ss_intronicA5SS_wt = 
		# 	(splicing_scores_intronicA5SS[[substr(ss_wt_elevenmer, 1, 6)]] + 
		# 	splicing_scores_intronicA5SS[[substr(ss_wt_elevenmer, 2, 7)]] + 
		# 	splicing_scores_intronicA5SS[[substr(ss_wt_elevenmer, 3, 8)]] + 
		# 	splicing_scores_intronicA5SS[[substr(ss_wt_elevenmer, 4, 9)]] + 
		# 	splicing_scores_intronicA5SS[[substr(ss_wt_elevenmer, 5, 10)]] + 
		# 	splicing_scores_intronicA5SS[[substr(ss_wt_elevenmer, 6, 11)]])/6) %>% 
		# mutate(ss_intronicA5SS = ss_intronicA5SS_mt-ss_intronicA5SS_wt) %>% 
		ungroup() %>% 
	# also get absolute values of predicted splicing effect
		mutate(ss_abs_EI = abs(ss_EI)) %>% 
		mutate(ss_abs_exonicA3SS = abs(ss_exonicA3SS)) %>% 
		mutate(ss_abs_exonicA5SS = abs(ss_exonicA5SS)) # %>% 
		# mutate(ss_abs_intronicA3SS = abs(ss_intronicA3SS)) %>% 
		# mutate(ss_abs_intronicA5SS = abs(ss_intronicA5SS))

	# flip if ANCDER
	if (ANCDER) {
		temp <- temp %>% 
			mutate(ss_EI = ifelse( 
				hub_variant_ALT==hub_variant_DER, ss_EI, -ss_EI)) %>% 
			mutate(ss_exonicA3SS = ifelse(
				hub_variant_ALT==hub_variant_DER, ss_exonicA3SS, -ss_exonicA3SS)) %>% 
			mutate(ss_exonicA5SS = ifelse(
				hub_variant_ALT==hub_variant_DER, ss_exonicA5SS, -ss_exonicA5SS)) # %>% 
			# mutate(ss_intronicA3SS = ifelse(
			# 	hub_variant_ALT==hub_variant_DER, ss_intronicA3SS, -ss_intronicA3SS)) %>% 
			# mutate(ss_intronicA5SS = ifelse(
			# 	hub_variant_ALT==hub_variant_DER, ss_intronicA5SS, -ss_intronicA5SS))
	}
	return(temp)
}
