#!/bin/R

# load packages
library(tidyverse)

# process data
for (i in 1:22) {
	# load hapR2 info
	file_name <- paste("../../results/linkage_hapR2_SNPs/ALL.chr", i, ".phase3_shapeit2_mvncall_integrated_v5a.20130502.hapR2.hap.ld.gz", sep="")
	hapR2_orig <- as_tibble(fread(file_name))

	# remove if same POS
	hapR2_orig <- hapR2_orig %>% 
		filter(POS1 != POS2)

	# flip and append SNPs
	hapR2_orig_flip <- hapR2_orig[,c(1,3,2,4:7)]
	names(hapR2_orig_flip) <- names(hapR2_orig)
	hapR2_orig <- bind_rows(hapR2_orig, hapR2_orig_flip)

	# summarise per first SNP
	hapR2_summ_temp <- hapR2_orig %>% 
		group_by(CHR, POS1) %>% 
		summarise(n=n()) %>% 
		ungroup() %>% 
		arrange(CHR, POS1)

	# save summary
	write_tsv(hapR2_summ_temp, gzfile(paste("../../results/linkage_hapR2_SNPs/ALL.chr", i, ".phase3_shapeit2_mvncall_integrated_v5a.20130502.hapR2.count.gz", sep="")))

	# clear memory
	hapR2_orig <- NULL
	hapR2_orig_flip <- NULL
}

# collate data
hapR2_summ_final <- NULL
for (i in 1:22) {
	# load summary
	hapR2_summ_temp <- as_tibble(fread(paste("../../results/linkage_hapR2_SNPs/ALL.chr", i, ".phase3_shapeit2_mvncall_integrated_v5a.20130502.hapR2.count.gz", sep="")))

	# append summary
	if (is.null(hapR2_summ_final)) {
		hapR2_summ_final <- hapR2_summ_temp
	} else {
		hapR2_summ_final <- bind_rows(hapR2_summ_final, hapR2_summ_temp)
	}
}

# save summary
write_tsv(hapR2_summ_final, gzfile(paste("../../results/linkage_hapR2_SNPs/ALL.chr.phase3_shapeit2_mvncall_integrated_v5a.20130502.hapR2.collated.gz", sep="")))
