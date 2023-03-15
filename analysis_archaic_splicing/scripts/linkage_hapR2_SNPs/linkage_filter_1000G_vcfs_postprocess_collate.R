#!/bin/R

# load packages
library(tidyverse)
library(data.table)

# collate data
hapR2_summ_final <- NULL
for (i in 1:22) {
	print(i)
	
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
