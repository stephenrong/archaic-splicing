#!/bin/R

library(tidyverse)
library(data.table)
library(rstatix)
library(ggpubr)
library(Hmisc)

# load tables
mapsy_variant_table <- as_tibble(fread("../../results/mapsy_to_variant_table_updated/Neanderthal_updated_mapsy_to_variant_table.txt.gz"))
mapsy_variant_table <- mapsy_variant_table %>% filter(!(exon_seqnames %in% c("chrX")))

final_v2_spliceai <- as_tibble(fread("../../results/annotate_splice_prediction/final_v2_variants_B_stat_mask_1KGP_archaic_gnomAD_introgr_hub_spliceai_gs_raw_dedup.txt.gz"))
final_v2_spliceai <- final_v2_spliceai %>% filter(!(hub_variant_CHROM %in% c("X")))
final_v2_spliceai <- final_v2_spliceai %>% filter(!is.na(spliceai_max))  # remove sites where spliceai does not make prediction
# final_v2_spliceai <- final_v2_spliceai %>% mutate(spliceai_max = ifelse(is.na(spliceai_max), 0, spliceai_max))  # fill in sites where spliceai does not make prediction

# join tables
mapsy_variant_table_spliceai <- left_join(mapsy_variant_table, final_v2_spliceai[,c(4,124:134)])

# visualize in bins
mapsy_variant_table_spliceai <- mapsy_variant_table_spliceai %>% 
	mutate(mpralm.ANCDER.logFC.bin = cut_width(mpralm.ANCDER.logFC, 2, boundary=1))
write_tsv(mapsy_variant_table_spliceai, gzfile("../../results/distribution_splice_prediction/mapsy_variant_table_spliceai_gs_raw.txt.gz"))

ggplot(mapsy_variant_table_spliceai %>% filter(!is.na(mpralm.ANCDER.logFC.bin))) + 
	geom_boxplot(aes(x=mpralm.ANCDER.logFC.bin, y=spliceai_max), fill="#beaed4") + 
	theme_pubr() + xlab(expression(paste("MaPSy functional score (log"[2], " FC) bin"))) + ylab("SpliceAI max score") + theme(aspect.ratio=0.3) + 
	theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
ggsave("../../results/distribution_splice_prediction/distribution_spliceai_gs_raw.pdf", scale=0.55)

# visualize by sigclass
ggplot(mapsy_variant_table_spliceai %>% filter(!is.na(mpralm.sigclass))) + 
	geom_boxplot(aes(x=mpralm.sigclass, y=spliceai_max), fill="#beaed4") + 
	theme_pubr() + xlab(expression(paste("MaPSy functional score (log"[2], " FC) bin"))) + ylab("SpliceAI max score") + theme(aspect.ratio=0.3) + 
	theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
ggsave("../../results/distribution_splice_prediction/distribution_spliceai_gs_raw_sigclass.pdf", scale=0.55)

# visualize in scatterplot
ggplot(mapsy_variant_table_spliceai %>% filter(!is.na(mpralm.ANCDER.logFC))) + 
	geom_hline(yintercept=0, color="red", alpha=0.5) + 
	geom_point(aes(x=mpralm.ANCDER.logFC, y=spliceai_max), fill="#beaed4", alpha=0.05) + 
	geom_smooth(aes(x=mpralm.ANCDER.logFC, y=spliceai_max), method="lm", formula=y~poly(x, 5, raw=TRUE)) + 
	theme_pubr() + xlab(expression(paste("MaPSy functional score (log"[2], " FC)"))) + ylab("SpliceAI max score") + theme(aspect.ratio=0.3) + 
	theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
ggsave("../../results/distribution_splice_prediction/distribution_spliceai_gs_raw_scatter.pdf", scale=0.55)

# evaluate lm models
aic_fits <- lapply(1:10, function(x) {AIC(lm(spliceai_max ~ poly(mpralm.ANCDER.logFC, x, raw=TRUE), data=mapsy_variant_table_spliceai))})
bic_fits <- lapply(1:10, function(x) {BIC(lm(spliceai_max ~ poly(mpralm.ANCDER.logFC, x, raw=TRUE), data=mapsy_variant_table_spliceai))})
pdf("../../results/distribution_splice_prediction/distribution_spliceai_gs_raw_scatter-BIC.pdf", height=4, width=4)
par(mar=c(5,6,4,2) + 0.1)
plot(1:10, aic_fits, type="b", pch=19, col="red", xlab="Polynomial regression degree", ylab="Bayesian Information\nCriterion (BIC)", main="SpliceAI max score")
dev.off()
