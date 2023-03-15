#!/bin/R

library(tidyverse)
library(data.table)
library(rstatix)
library(ggpubr)
library(Hmisc)

# load tables
mapsy_variant_table <- as_tibble(fread("../../results/mapsy_to_variant_table_updated/Neanderthal_updated_mapsy_to_variant_table.txt.gz"))
mapsy_variant_table <- mapsy_variant_table %>% filter(!(exon_seqnames %in% c("chrX")))

# get hexamer scores
source("get_hexamer_scores.R")
mapsy_variant_table_hexamer <- mapsy_variant_table %>% 
	mutate(
		ss_reference_genome = hub_reference_genome,
		ss_variant_CHROM = hub_variant_CHROM, 
		ss_variant_POS = hub_variant_POS, 
		ss_variant_REF = hub_variant_REF,
		ss_variant_ALT = hub_variant_ALT,
		ss_variant_ANC = hub_variant_ANC,
		ss_variant_DER = hub_variant_DER,
		ss_strand = exon_strand
	) %>% get_hexamer_scores(ANCDER=TRUE)

# visualize in bins
mapsy_variant_table_hexamer <- mapsy_variant_table_hexamer %>% 
	mutate(mpralm.ANCDER.logFC.bin = cut_width(mpralm.ANCDER.logFC, 2, boundary=1))
write_tsv(mapsy_variant_table_hexamer, gzfile("../../results/distribution_splice_prediction/mapsy_variant_table_hexamer.txt.gz"))

ggplot(mapsy_variant_table_hexamer %>% filter(!is.na(mpralm.ANCDER.logFC.bin))) + 
	geom_boxplot(aes(x=mpralm.ANCDER.logFC.bin, y=ss_EI), fill="#beaed4") + 
	theme_pubr() + xlab(expression(paste("MaPSy functional score (log"[2], " FC) bin"))) + ylab("Delta EI score") + theme(aspect.ratio=0.3) + 
	theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
ggsave("../../results/distribution_splice_prediction/distribution_ss_EI.pdf", scale=0.55)

ggplot(mapsy_variant_table_hexamer %>% filter(!is.na(mpralm.ANCDER.logFC.bin))) + 
	geom_boxplot(aes(x=mpralm.ANCDER.logFC.bin, y=ss_exonicA3SS), fill="#beaed4") + 
	theme_pubr() + xlab(expression(paste("MaPSy functional score (log"[2], " FC) bin"))) + ylab("Delta exonic A3SS score") + theme(aspect.ratio=0.3) + 
	theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
ggsave("../../results/distribution_splice_prediction/distribution_ss_exonicA3SS.pdf", scale=0.55)

ggplot(mapsy_variant_table_hexamer %>% filter(!is.na(mpralm.ANCDER.logFC.bin))) + 
	geom_boxplot(aes(x=mpralm.ANCDER.logFC.bin, y=ss_exonicA5SS), fill="#beaed4") + 
	theme_pubr() + xlab(expression(paste("MaPSy functional score (log"[2], " FC) bin"))) + ylab("Delta exonic A5SS score") + theme(aspect.ratio=0.3) + 
	theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
ggsave("../../results/distribution_splice_prediction/distribution_ss_exonicA5SS.pdf", scale=0.55)

# ggplot(mapsy_variant_table_hexamer %>% filter(!is.na(mpralm.ANCDER.logFC.bin))) + 
# 	geom_boxplot(aes(x=mpralm.ANCDER.logFC.bin, y=ss_intronicA3SS), fill="#beaed4") + 
# 	theme_pubr() + xlab(expression(paste("MaPSy functional score (log"[2], " FC) bin"))) + ylab("Delta intronic A3SS score") + theme(aspect.ratio=0.3) + 
# 	theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
# ggsave("../../results/distribution_splice_prediction/distribution_ss_intronicA3SS.pdf", scale=0.55)

# ggplot(mapsy_variant_table_hexamer %>% filter(!is.na(mpralm.ANCDER.logFC.bin))) + 
# 	geom_boxplot(aes(x=mpralm.ANCDER.logFC.bin, y=ss_intronicA5SS), fill="#beaed4") + 
# 	theme_pubr() + xlab(expression(paste("MaPSy functional score (log"[2], " FC) bin"))) + ylab("Delta intronic A5SS score") + theme(aspect.ratio=0.3) + 
# 	theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
# ggsave("../../results/distribution_splice_prediction/distribution_ss_intronicA5SS.pdf", scale=0.55)

# visualize by sigclass
ggplot(mapsy_variant_table_hexamer %>% filter(!is.na(mpralm.sigclass))) + 
	geom_boxplot(aes(x=mpralm.sigclass, y=ss_EI), fill="#beaed4") + 
	theme_pubr() + xlab(expression(paste("MaPSy functional score (log"[2], " FC) bin"))) + ylab("Delta EI score") + theme(aspect.ratio=0.3) + 
	theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
ggsave("../../results/distribution_splice_prediction/distribution_ss_EI_sigclass.pdf", scale=0.55)

ggplot(mapsy_variant_table_hexamer %>% filter(!is.na(mpralm.sigclass))) + 
	geom_boxplot(aes(x=mpralm.sigclass, y=ss_exonicA3SS), fill="#beaed4") + 
	theme_pubr() + xlab(expression(paste("MaPSy functional score (log"[2], " FC) bin"))) + ylab("Delta exonic A3SS score") + theme(aspect.ratio=0.3) + 
	theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
ggsave("../../results/distribution_splice_prediction/distribution_ss_exonicA3SS_sigclass.pdf", scale=0.55)

ggplot(mapsy_variant_table_hexamer %>% filter(!is.na(mpralm.sigclass))) + 
	geom_boxplot(aes(x=mpralm.sigclass, y=ss_exonicA5SS), fill="#beaed4") + 
	theme_pubr() + xlab(expression(paste("MaPSy functional score (log"[2], " FC) bin"))) + ylab("Delta exonic A5SS score") + theme(aspect.ratio=0.3) + 
	theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
ggsave("../../results/distribution_splice_prediction/distribution_ss_exonicA5SS_sigclass.pdf", scale=0.55)

# ggplot(mapsy_variant_table_hexamer %>% filter(!is.na(mpralm.sigclass))) + 
# 	geom_boxplot(aes(x=mpralm.sigclass, y=ss_intronicA3SS), fill="#beaed4") + 
# 	theme_pubr() + xlab(expression(paste("MaPSy functional score (log"[2], " FC) bin"))) + ylab("Delta intronic A3SS score") + theme(aspect.ratio=0.3) + 
# 	theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
# ggsave("../../results/distribution_splice_prediction/distribution_ss_intronicA3SS_sigclass.pdf", scale=0.55)

# ggplot(mapsy_variant_table_hexamer %>% filter(!is.na(mpralm.sigclass))) + 
# 	geom_boxplot(aes(x=mpralm.sigclass, y=ss_intronicA5SS), fill="#beaed4") + 
# 	theme_pubr() + xlab(expression(paste("MaPSy functional score (log"[2], " FC) bin"))) + ylab("Delta intronic A5SS score") + theme(aspect.ratio=0.3) + 
# 	theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
# ggsave("../../results/distribution_splice_prediction/distribution_ss_intronicA5SS_sigclass.pdf", scale=0.55)

# visualize in scatterplot
ggplot(mapsy_variant_table_hexamer %>% filter(!is.na(mpralm.ANCDER.logFC))) + 
	geom_hline(yintercept=0, color="red", alpha=0.5) + 
	geom_point(aes(x=mpralm.ANCDER.logFC, y=ss_EI), fill="#beaed4", alpha=0.05) + 
	# geom_smooth(aes(x=mpralm.ANCDER.logFC, y=ss_EI), method="lm", formula=y~poly(x, 1, raw=TRUE), color="green") + 
	geom_smooth(aes(x=mpralm.ANCDER.logFC, y=ss_EI), method="lm", formula=y~poly(x, 4, raw=TRUE), color="blue") + 
	theme_pubr() + xlab(expression(paste("MaPSy functional score (log"[2], " FC)"))) + ylab("Delta EI score") + theme(aspect.ratio=0.3) + 
	theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
ggsave("../../results/distribution_splice_prediction/distribution_ss_EI_scatter.pdf", scale=0.55)

ggplot(mapsy_variant_table_hexamer %>% filter(!is.na(mpralm.ANCDER.logFC))) + 
	geom_hline(yintercept=0, color="red", alpha=0.5) + 
	geom_point(aes(x=mpralm.ANCDER.logFC, y=ss_exonicA3SS), fill="#beaed4", alpha=0.05) + 
	# geom_smooth(aes(x=mpralm.ANCDER.logFC, y=ss_exonicA3SS), method="lm", formula=y~poly(x, 1, raw=TRUE), color="green") + 
	geom_smooth(aes(x=mpralm.ANCDER.logFC, y=ss_exonicA3SS), method="lm", formula=y~poly(x, 4, raw=TRUE), color="blue") + 
	theme_pubr() + xlab(expression(paste("MaPSy functional score (log"[2], " FC)"))) + ylab("Delta exonic A3SS score") + theme(aspect.ratio=0.3) + 
	theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
ggsave("../../results/distribution_splice_prediction/distribution_ss_exonicA3SS_scatter.pdf", scale=0.55)

ggplot(mapsy_variant_table_hexamer %>% filter(!is.na(mpralm.ANCDER.logFC))) + 
	geom_hline(yintercept=0, color="red", alpha=0.5) + 
	geom_point(aes(x=mpralm.ANCDER.logFC, y=ss_exonicA5SS), fill="#beaed4", alpha=0.05) + 
	# geom_smooth(aes(x=mpralm.ANCDER.logFC, y=ss_exonicA5SS), method="lm", formula=y~poly(x, 1, raw=TRUE), color="green") + 
	geom_smooth(aes(x=mpralm.ANCDER.logFC, y=ss_exonicA5SS), method="lm", formula=y~poly(x, 4, raw=TRUE), color="blue") + 
	theme_pubr() + xlab(expression(paste("MaPSy functional score (log"[2], " FC)"))) + ylab("Delta exonic A5SS score") + theme(aspect.ratio=0.3) + 
	theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
ggsave("../../results/distribution_splice_prediction/distribution_ss_exonicA5SS_scatter.pdf", scale=0.55)

# ggplot(mapsy_variant_table_hexamer %>% filter(!is.na(mpralm.ANCDER.logFC))) + 
# 	geom_hline(yintercept=0, color="red", alpha=0.5) + 
# 	geom_point(aes(x=mpralm.ANCDER.logFC, y=ss_intronicA3SS), fill="#beaed4", alpha=0.05) + 
# 	# geom_smooth(aes(x=mpralm.ANCDER.logFC, y=ss_intronicA3SS), method="lm", formula=y~poly(x, 1, raw=TRUE), color="green") + 
# 	geom_smooth(aes(x=mpralm.ANCDER.logFC, y=ss_intronicA3SS), method="lm", formula=y~poly(x, 4, raw=TRUE), color="blue") + 
# 	theme_pubr() + xlab(expression(paste("MaPSy functional score (log"[2], " FC)"))) + ylab("Delta intronic A3SS score") + theme(aspect.ratio=0.3) + 
# 	theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
# ggsave("../../results/distribution_splice_prediction/distribution_ss_intronicA3SS_scatter.pdf", scale=0.55)

# ggplot(mapsy_variant_table_hexamer %>% filter(!is.na(mpralm.ANCDER.logFC))) + 
# 	geom_hline(yintercept=0, color="red", alpha=0.5) + 
# 	geom_point(aes(x=mpralm.ANCDER.logFC, y=ss_intronicA5SS), fill="#beaed4", alpha=0.05) + 
# 	# geom_smooth(aes(x=mpralm.ANCDER.logFC, y=ss_intronicA5SS), method="lm", formula=y~poly(x, 1, raw=TRUE), color="green") + 
# 	geom_smooth(aes(x=mpralm.ANCDER.logFC, y=ss_intronicA5SS), method="lm", formula=y~poly(x, 4, raw=TRUE), color="blue") + 
# 	theme_pubr() + xlab(expression(paste("MaPSy functional score (log"[2], " FC)"))) + ylab("Delta intronic A5SS score") + theme(aspect.ratio=0.3) + 
# 	theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
# ggsave("../../results/distribution_splice_prediction/distribution_ss_intronicA5SS_scatter.pdf", scale=0.55)

# evaluate lm models
aic_fits <- lapply(1:10, function(x) {AIC(lm(ss_EI ~ poly(mpralm.ANCDER.logFC, x, raw=TRUE), data=mapsy_variant_table_hexamer))})
bic_fits <- lapply(1:10, function(x) {BIC(lm(ss_EI ~ poly(mpralm.ANCDER.logFC, x, raw=TRUE), data=mapsy_variant_table_hexamer))})
pdf("../../results/distribution_splice_prediction/distribution_ss_EI_scatter-BIC.pdf", height=4, width=4)
par(mar=c(5,6,4,2) + 0.1)
plot(1:10, aic_fits, type="b", pch=19, col="red", xlab="Polynomial regression degree", ylab="Bayesian Information\nCriterion (BIC)", main="Delta EI score")
dev.off()

aic_fits <- lapply(1:10, function(x) {AIC(lm(ss_exonicA3SS ~ poly(mpralm.ANCDER.logFC, x, raw=TRUE), data=mapsy_variant_table_hexamer))})
bic_fits <- lapply(1:10, function(x) {BIC(lm(ss_exonicA3SS ~ poly(mpralm.ANCDER.logFC, x, raw=TRUE), data=mapsy_variant_table_hexamer))})
pdf("../../results/distribution_splice_prediction/distribution_ss_exonicA3SS_scatter-BIC.pdf", height=4, width=4)
par(mar=c(5,6,4,2) + 0.1)
plot(1:10, aic_fits, type="b", pch=19, col="red", xlab="Polynomial regression degree", ylab="Bayesian Information\nCriterion (BIC)", main="Delta exonic A3SS score")
dev.off()

aic_fits <- lapply(1:10, function(x) {AIC(lm(ss_exonicA5SS ~ poly(mpralm.ANCDER.logFC, x, raw=TRUE), data=mapsy_variant_table_hexamer))})
bic_fits <- lapply(1:10, function(x) {BIC(lm(ss_exonicA5SS ~ poly(mpralm.ANCDER.logFC, x, raw=TRUE), data=mapsy_variant_table_hexamer))})
pdf("../../results/distribution_splice_prediction/distribution_ss_exonicA5SS_scatter-BIC.pdf", height=4, width=4)
par(mar=c(5,6,4,2) + 0.1)
plot(1:10, aic_fits, type="b", pch=19, col="red", xlab="Polynomial regression degree", ylab="Bayesian Information\nCriterion (BIC)", main="Delta exonic A5SS score")
dev.off()

# visualize in by position
mapsy_variant_table_hexamer <- mapsy_variant_table_hexamer %>% 
	mutate(variant_exon_pos_in_from_3ss = ifelse(exon_strand == "+", variant_start-exon_start+1, (exon_end-variant_end+1))) %>% 
	mutate(variant_exon_pos_in_from_5ss = exon_width-variant_exon_pos_in_from_3ss+1) %>% 
	mutate(`Variant within 5 bp of 3' ss or 5' ss` = ifelse((variant_exon_pos_in_from_3ss<5)|(variant_exon_pos_in_from_5ss<5), "TRUE", "FALSE"))

ggplot(mapsy_variant_table_hexamer %>% filter(!is.na(mpralm.ANCDER.logFC))) + 
	geom_hline(yintercept=0, color="red", alpha=0.2) + 
	geom_point(aes(x=mpralm.ANCDER.logFC, y=ss_EI, color=`Variant within 5 bp of 3' ss or 5' ss`), alpha=0.2) + 
	# geom_smooth(aes(x=mpralm.ANCDER.logFC, y=ss_EI), method="lm", formula=y~poly(x, 1, raw=TRUE)), color="green") + 
	geom_smooth(aes(x=mpralm.ANCDER.logFC, y=ss_EI, color=`Variant within 5 bp of 3' ss or 5' ss`), method="lm", formula=y~poly(x, 1, raw=TRUE), color="blue") + 
	theme_pubr() + xlab(expression(paste("MaPSy functional score (log"[2], " FC)"))) + ylab("Delta EI score") + theme(aspect.ratio=0.3) + 
	theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + 
	facet_wrap(~`Variant within 5 bp of 3' ss or 5' ss`)
ggsave("../../results/distribution_splice_prediction/distribution_ss_EI_scatter.pdf", scale=0.7)

ggplot(mapsy_variant_table_hexamer %>% filter(!is.na(mpralm.ANCDER.logFC))) + 
	geom_hline(yintercept=0, color="red", alpha=0.2) + 
	geom_point(aes(x=mpralm.ANCDER.logFC, y=ss_exonicA3SS, color=`Variant within 5 bp of 3' ss or 5' ss`), alpha=0.2) + 
	# geom_smooth(aes(x=mpralm.ANCDER.logFC, y=ss_exonicA3SS), method="lm", formula=y~poly(x, 1, raw=TRUE)), color="green") + 
	geom_smooth(aes(x=mpralm.ANCDER.logFC, y=ss_exonicA3SS, color=`Variant within 5 bp of 3' ss or 5' ss`), method="lm", formula=y~poly(x, 1, raw=TRUE), color="blue") + 
	theme_pubr() + xlab(expression(paste("MaPSy functional score (log"[2], " FC)"))) + ylab("Delta exonic A3SS score") + theme(aspect.ratio=0.3) + 
	theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + 
	facet_wrap(~`Variant within 5 bp of 3' ss or 5' ss`)
ggsave("../../results/distribution_splice_prediction/distribution_ss_exonicA3SS_scatter.pdf", scale=0.7)

ggplot(mapsy_variant_table_hexamer %>% filter(!is.na(mpralm.ANCDER.logFC))) + 
	geom_hline(yintercept=0, color="red", alpha=0.2) + 
	geom_point(aes(x=mpralm.ANCDER.logFC, y=ss_exonicA5SS, color=`Variant within 5 bp of 3' ss or 5' ss`), alpha=0.2) + 
	# geom_smooth(aes(x=mpralm.ANCDER.logFC, y=ss_exonicA5SS), method="lm", formula=y~poly(x, 1, raw=TRUE)), color="green") + 
	geom_smooth(aes(x=mpralm.ANCDER.logFC, y=ss_exonicA5SS, color=`Variant within 5 bp of 3' ss or 5' ss`), method="lm", formula=y~poly(x, 1, raw=TRUE), color="blue") + 
	theme_pubr() + xlab(expression(paste("MaPSy functional score (log"[2], " FC)"))) + ylab("Delta exonic A5SS score") + theme(aspect.ratio=0.3) + 
	theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + 
	facet_wrap(~`Variant within 5 bp of 3' ss or 5' ss`)
ggsave("../../results/distribution_splice_prediction/distribution_ss_exonicA5SS_scatter.pdf", scale=0.7)

# ggplot(mapsy_variant_table_hexamer %>% filter(!is.na(mpralm.ANCDER.logFC))) + 
# 	geom_hline(yintercept=0, color="red", alpha=0.2) + 
# 	geom_point(aes(x=mpralm.ANCDER.logFC, y=ss_intronicA3SS, color=`Variant within 5 bp of 3' ss or 5' ss`), alpha=0.2) + 
# 	# geom_smooth(aes(x=mpralm.ANCDER.logFC, y=ss_intronicA3SS), method="lm", formula=y~poly(x, 1, raw=TRUE)), color="green") + 
# 	geom_smooth(aes(x=mpralm.ANCDER.logFC, y=ss_intronicA3SS, color=`Variant within 5 bp of 3' ss or 5' ss`), method="lm", formula=y~poly(x, 1, raw=TRUE), color="blue") + 
# 	theme_pubr() + xlab(expression(paste("MaPSy functional score (log"[2], " FC)"))) + ylab("Delta intronic A3SS score") + theme(aspect.ratio=0.3) + 
# 	theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + 
# 	facet_wrap(~`Variant within 5 bp of 3' ss or 5' ss`)
# ggsave("../../results/distribution_splice_prediction/distribution_ss_intronicA3SS_scatter.pdf", scale=0.7)

# ggplot(mapsy_variant_table_hexamer %>% filter(!is.na(mpralm.ANCDER.logFC))) + 
# 	geom_hline(yintercept=0, color="red", alpha=0.2) + 
# 	geom_point(aes(x=mpralm.ANCDER.logFC, y=ss_intronicA5SS, color=`Variant within 5 bp of 3' ss or 5' ss`), alpha=0.2) + 
# 	# geom_smooth(aes(x=mpralm.ANCDER.logFC, y=ss_intronicA5SS), method="lm", formula=y~poly(x, 1, raw=TRUE)), color="green") + 
# 	geom_smooth(aes(x=mpralm.ANCDER.logFC, y=ss_intronicA5SS, color=`Variant within 5 bp of 3' ss or 5' ss`), method="lm", formula=y~poly(x, 1, raw=TRUE), color="blue") + 
# 	theme_pubr() + xlab(expression(paste("MaPSy functional score (log"[2], " FC)"))) + ylab("Delta intronic A5SS score") + theme(aspect.ratio=0.3) + 
# 	theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + 
# 	facet_wrap(~`Variant within 5 bp of 3' ss or 5' ss`)
# ggsave("../../results/distribution_splice_prediction/distribution_ss_intronicA5SS_scatter.pdf", scale=0.7)
