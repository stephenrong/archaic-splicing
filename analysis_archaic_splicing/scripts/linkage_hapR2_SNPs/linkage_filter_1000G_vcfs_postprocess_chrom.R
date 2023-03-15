#!/bin/R

# load packages
library(tidyverse)
library(data.table)

# loop chromosomes
for (i in 12:1) {
	print(i)

	# get file name parts
	file_list <- list.files(path = "../../results/linkage_hapR2_SNPs/", pattern = paste0("ALL.chr", i, ".phase3_shapeit2_mvncall_integrated_v5a.20130502.hapR2.hap.ld.gz.part*"))

	# loop chromosome parts
	hapR2_summ_list <- list()
	iter <- 1
	for (file_name in file_list) {
		hapR2_orig <- fread(paste0("../../results/linkage_hapR2_SNPs/", file_name))
		names(hapR2_orig) <- c("CHR", "POS1", "POS2", "N_CHR", "R^2", "D", "Dprime")

		# remove if same POS
		hapR2_orig <- hapR2_orig %>% 
			filter(POS1 != POS2)

		# remove extra cols
		hapR2_orig <- hapR2_orig[,c(1,2,3,5)]

		# flip and append SNPs
		hapR2_orig_flip <- hapR2_orig[,c(1,3,2,4)]
		names(hapR2_orig_flip) <- names(hapR2_orig)
		hapR2_orig <- bind_rows(hapR2_orig, hapR2_orig_flip)
		# hapR2_orig_flip <- NULL

		# hapR2_orig <- data.table(hapR2_orig)
		hapR2_summ_temp <- hapR2_orig[,.(n=sum(`R^2`>0.2), s=sum(`R^2`)), by=.(CHR, POS1)]
		# hapR2_orig <- NULL

		hapR2_summ_list[[iter]] <- hapR2_summ_temp
		# hapR2_summ_temp <- NULL

		iter <- iter + 1
	}

	# save summary
	# 	parts work because we are computing sums
	hapR2_summ_final_temp <- bind_rows(hapR2_summ_list)
	hapR2_summ_final <- hapR2_summ_final_temp[,.(n=sum(n), s=sum(s)), by=.(CHR, POS1)]
	write_tsv(hapR2_summ_final, gzfile(paste("../../results/linkage_hapR2_SNPs/ALL.chr", i, ".phase3_shapeit2_mvncall_integrated_v5a.20130502.hapR2.count.gz", sep="")))
	# hapR2_summ_final <- NULL
}
