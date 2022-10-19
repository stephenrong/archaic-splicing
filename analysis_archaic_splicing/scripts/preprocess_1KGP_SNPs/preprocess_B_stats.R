#!/bin/R

library(tidyverse)
library(data.table)
library(plyranges)
library(rtracklayer)

# convert to usable format, such as granges
chr_list <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX")

B_stats_hg18_gr <- NULL
for (chr in chr_list) {
	print(chr)
	chr_bkgd <- as_tibble(fread(paste("../../data/annotate_B_statistics/bkgd/", chr, ".bkgd", sep="")))
	chr_bkgd$V3 <- cumsum(chr_bkgd$V2)  # get ends from cumulative sum
	chr_bkgd$V4 <- c(1, (chr_bkgd$V3[-length(chr_bkgd$V3)]+1))  # get starts from ends minus one
	chr_bkgd$V5 <- chr  # get seqnames from loop iteration
	chr_bkgd <- chr_bkgd %>%  # get columns needed for conversion to GRange
		dplyr::select(V5, V4, V3, V1)  # seqnames, start, end, B_statistic
	names(chr_bkgd) <- c("seqnames", "start", "end", "B_statistic")
	chr_bkgd$B_statistic_bin <-  # get B_statistic decile bins
		as.numeric(cut_number(chr_bkgd$B_statistic, 10))
	chr_bkgd <- chr_bkgd %>% GRanges()
	if (is.null(B_stats_hg18_gr)) {
		B_stats_hg18_gr <- chr_bkgd
	} else {
		B_stats_hg18_gr <- c(B_stats_hg18_gr, chr_bkgd)
	}
}

# lift over from hg18 to hg19
hg18ToHg19 <- import.chain("../../data/annotate_lift_over/hg18ToHg19.over.chain")
B_stats_hg19_gr <- B_stats_hg18_gr %>% liftOver(hg18ToHg19) %>% unlist() %>% 
	as_tibble() %>% filter(seqnames %in% chr_list) %>%  # corner case of mapping to chrY
	mutate(seqnames = factor(seqnames, levels=chr_list)) %>% 
	arrange(seqnames, start, end)
write_tsv(B_stats_hg19_gr, gzfile("../../results/preprocess_1KGP_SNPs/B_stats_hg19_gr.txt.gz"))
