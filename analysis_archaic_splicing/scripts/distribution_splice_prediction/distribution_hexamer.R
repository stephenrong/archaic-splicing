#!/bin/R

library(tidyverse)
library(data.table)
library(rstatix)
library(ggpubr)
library(Hmisc)

# load tables
mapsy_variant_table <- as_tibble(fread("../../results/mapsy_to_variant_table_updated/Neanderthal_updated_mapsy_to_variant_table.txt.gz"))

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

ggplot(mapsy_variant_table_hexamer %>% filter(!is.na(mpralm.ANCDER.logFC.bin))) + 
	geom_boxplot(aes(x=mpralm.ANCDER.logFC.bin, y=ss_intronicA3SS), fill="#beaed4") + 
	theme_pubr() + xlab(expression(paste("MaPSy functional score (log"[2], " FC) bin"))) + ylab("Delta intronic A3SS score") + theme(aspect.ratio=0.3) + 
	theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
ggsave("../../results/distribution_splice_prediction/distribution_ss_intronicA3SS.pdf", scale=0.55)

ggplot(mapsy_variant_table_hexamer %>% filter(!is.na(mpralm.ANCDER.logFC.bin))) + 
	geom_boxplot(aes(x=mpralm.ANCDER.logFC.bin, y=ss_intronicA5SS), fill="#beaed4") + 
	theme_pubr() + xlab(expression(paste("MaPSy functional score (log"[2], " FC) bin"))) + ylab("Delta intronic A5SS score") + theme(aspect.ratio=0.3) + 
	theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
ggsave("../../results/distribution_splice_prediction/distribution_ss_intronicA5SS.pdf", scale=0.55)
