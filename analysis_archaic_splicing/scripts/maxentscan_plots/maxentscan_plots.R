#!/bin/R

library(tidyverse)
library(data.table)
library(VarCon)
library(BSgenome.Hsapiens.UCSC.hg19)
hg19 <- BSgenome.Hsapiens.UCSC.hg19


# Genome-wide density

# load gene annotations
canonical_exons <- as_tibble(fread("../../../final-mapsy-geisinger/data/known_canonical/gencode_v32_lift37_basic_canonical_exons.txt.gz"))
canonical_exons_internal <- canonical_exons %>% 
	group_by(canonical_transcript_id) %>% 
	filter(!(canonical_exon_number %in% c(1, max(canonical_exon_number)))) %>% 
	ungroup()

# hard code common sequences
three_A <- "ATGgtgagt"
three_B <- "AAGgtatgg"
five_A <- "ctccttgcctcttcttgtagTCC"
five_B <- "taatttcatatttcccccagGCA"

maxent_three_A <- as.numeric(calculateMaxEntScanScore(three_A, 5))
maxent_three_B <- as.numeric(calculateMaxEntScanScore(three_B, 5))
maxent_five_A <- as.numeric(calculateMaxEntScanScore(five_A, 3))
maxent_five_B <- as.numeric(calculateMaxEntScanScore(five_B, 3))

color_three_A <- "#00AEEf"
color_three_B <- "#C43974"
color_five_A <- "#C79A34"
color_five_B <- "#893B90"

# helper function
expandRange = function(x, upstream, downstream) {
	strand_is_minus = strand(x) == "-"
	on_plus = which(!strand_is_minus)
	on_minus = which(strand_is_minus)
	start(x)[on_plus] = start(x)[on_plus] - upstream
	start(x)[on_minus] = start(x)[on_minus] - downstream
	end(x)[on_plus] = end(x)[on_plus] + downstream
	end(x)[on_minus] = end(x)[on_minus] + upstream
	x
}

# get three prime splice sites
canonical_exons_internal_3ss <- canonical_exons_internal %>% 
	dplyr::select(seqnames, start, end, width, strand) %>% 
	mutate(start = ifelse(strand == "+", start, end+1)) %>% 
	mutate(end = ifelse(strand == "+", start, end+1)-1) %>% 
	mutate(width = end - start)

canonical_exons_internal_3ss_ext <- canonical_exons_internal_3ss %>% 
	GRanges() %>% expandRange(upstream = 20, downstream = 3)

canonical_exons_internal_3ss_ext_seqs <- getSeq(hg19, canonical_exons_internal_3ss_ext)
canonical_exons_internal_3ss_ext_maxent <- calculateMaxEntScanScore(canonical_exons_internal_3ss_ext_seqs, 3)
canonical_exons_internal_3ss_ext_tb <- as_tibble(canonical_exons_internal_3ss_ext)
canonical_exons_internal_3ss_ext_tb$seqs <- as.character(canonical_exons_internal_3ss_ext_seqs)
canonical_exons_internal_3ss_ext_tb$maxent <- as.numeric(canonical_exons_internal_3ss_ext_maxent)
write_tsv(canonical_exons_internal_3ss_ext_tb, gzfile("../../results/maxentscan_plots/canonical_exons_internal_3ss_ext_tb.txt"))

# get five prime splice sites
canonical_exons_internal_5ss <- canonical_exons_internal %>% 
	dplyr::select(seqnames, start, end, width, strand) %>% 
	mutate(start = ifelse(strand == "+", end+1, start)) %>% 
	mutate(end = ifelse(strand == "+", end+1, start)-1) %>% 
	mutate(width = end - start)

canonical_exons_internal_5ss_ext <- canonical_exons_internal_5ss %>% 
	GRanges() %>% expandRange(upstream = 3, downstream = 6)

canonical_exons_internal_5ss_ext_seqs <- getSeq(hg19, canonical_exons_internal_5ss_ext)
canonical_exons_internal_5ss_ext_maxent <- calculateMaxEntScanScore(canonical_exons_internal_5ss_ext_seqs, 5)
canonical_exons_internal_5ss_ext_tb <- as_tibble(canonical_exons_internal_5ss_ext)
canonical_exons_internal_5ss_ext_tb$seqs <- as.character(canonical_exons_internal_5ss_ext_seqs)
canonical_exons_internal_5ss_ext_tb$maxent <- as.numeric(canonical_exons_internal_5ss_ext_maxent)
write_tsv(canonical_exons_internal_5ss_ext_tb, gzfile("../../results/maxentscan_plots/canonical_exons_internal_5ss_ext_tb.txt"))

# visualize
ggplot(canonical_exons_internal_3ss_ext_tb) + theme_classic() + 
	geom_density(aes(maxent)) + 
	geom_segment(aes(x=maxent_five_A, xend=maxent_five_A, y=0.19, yend=0.17), arrow = arrow(length = unit(0.2, "cm"), type="closed"), color=color_five_A) + 
	geom_segment(aes(x=maxent_five_B, xend=maxent_five_B, y=0.18, yend=0.16), arrow = arrow(length = unit(0.2, "cm"), type="closed"), color=color_five_B) + 
	xlab("MaxEntScan 3ss score") + ylab("Genome-wide density") + theme(aspect.ratio=0.9)
ggsave("../../results/maxentscan_plots/canonical_exons_internal_3ss_ext_tb.pdf", scale=0.3)
# ggsave("../../results/maxentscan_plots/canonical_exons_internal_3ss_ext_tb.png", scale=0.3)

ggplot(canonical_exons_internal_5ss_ext_tb) + theme_classic() + 
	geom_density(aes(maxent)) + 
	geom_segment(aes(x=maxent_three_A, xend=maxent_three_A, y=0.26, yend=0.24), arrow = arrow(length = unit(0.2, "cm"), type="closed"), color=color_three_A) + 
	geom_segment(aes(x=maxent_three_B, xend=maxent_three_B, y=0.28, yend=0.26), arrow = arrow(length = unit(0.2, "cm"), type="closed"), color=color_three_B) + 
	xlab("MaxEntScan 5ss score") + ylab("Genome-wide density") + theme(aspect.ratio=0.9)
ggsave("../../results/maxentscan_plots/canonical_exons_internal_5ss_ext_tb.pdf", scale=0.3)
# ggsave("../../results/maxentscan_plots/canonical_exons_internal_5ss_ext_tb.png", scale=0.3)


# MaPSy-wide density

# load gene annotations
mapsywide_exons <- as_tibble(fread("../../results/mapsy_to_variant_table_updated/Neanderthal_updated_mapsy_to_variant_table.txt.gz"))
mapsywide_exons_internal <- mapsywide_exons %>% 
	dplyr::select(exon_seqnames, exon_start, exon_end, exon_width, exon_strand)
names(mapsywide_exons_internal) <- c("seqnames", "start", "end", "width", "strand")

# hard code common sequences
three_A <- "ATGgtgagt"
three_B <- "AAGgtatgg"
five_A <- "ctccttgcctcttcttgtagTCC"
five_B <- "taatttcatatttcccccagGCA"

maxent_three_A <- as.numeric(calculateMaxEntScanScore(three_A, 5))
maxent_three_B <- as.numeric(calculateMaxEntScanScore(three_B, 5))
maxent_five_A <- as.numeric(calculateMaxEntScanScore(five_A, 3))
maxent_five_B <- as.numeric(calculateMaxEntScanScore(five_B, 3))

color_three_A <- "#00AEEf"
color_three_B <- "#C43974"
color_five_A <- "#C79A34"
color_five_B <- "#893B90"

# helper function
expandRange = function(x, upstream, downstream) {
	strand_is_minus = strand(x) == "-"
	on_plus = which(!strand_is_minus)
	on_minus = which(strand_is_minus)
	start(x)[on_plus] = start(x)[on_plus] - upstream
	start(x)[on_minus] = start(x)[on_minus] - downstream
	end(x)[on_plus] = end(x)[on_plus] + downstream
	end(x)[on_minus] = end(x)[on_minus] + upstream
	x
}

# get three prime splice sites
mapsywide_exons_internal_3ss <- mapsywide_exons_internal %>% 
	dplyr::select(seqnames, start, end, width, strand) %>% 
	mutate(start = ifelse(strand == "+", start, end+1)) %>% 
	mutate(end = ifelse(strand == "+", start, end+1)-1) %>% 
	mutate(width = end - start)

mapsywide_exons_internal_3ss_ext <- mapsywide_exons_internal_3ss %>% 
	GRanges() %>% expandRange(upstream = 20, downstream = 3)

mapsywide_exons_internal_3ss_ext_seqs <- getSeq(hg19, mapsywide_exons_internal_3ss_ext)
mapsywide_exons_internal_3ss_ext_maxent <- calculateMaxEntScanScore(mapsywide_exons_internal_3ss_ext_seqs, 3)
mapsywide_exons_internal_3ss_ext_tb <- as_tibble(mapsywide_exons_internal_3ss_ext)
mapsywide_exons_internal_3ss_ext_tb$seqs <- as.character(mapsywide_exons_internal_3ss_ext_seqs)
mapsywide_exons_internal_3ss_ext_tb$maxent <- as.numeric(mapsywide_exons_internal_3ss_ext_maxent)
write_tsv(mapsywide_exons_internal_3ss_ext_tb, gzfile("../../results/maxentscan_plots/mapsywide_exons_internal_3ss_ext_tb.txt"))

# get five prime splice sites
mapsywide_exons_internal_5ss <- mapsywide_exons_internal %>% 
	dplyr::select(seqnames, start, end, width, strand) %>% 
	mutate(start = ifelse(strand == "+", end+1, start)) %>% 
	mutate(end = ifelse(strand == "+", end+1, start)-1) %>% 
	mutate(width = end - start)

mapsywide_exons_internal_5ss_ext <- mapsywide_exons_internal_5ss %>% 
	GRanges() %>% expandRange(upstream = 3, downstream = 6)

mapsywide_exons_internal_5ss_ext_seqs <- getSeq(hg19, mapsywide_exons_internal_5ss_ext)
mapsywide_exons_internal_5ss_ext_maxent <- calculateMaxEntScanScore(mapsywide_exons_internal_5ss_ext_seqs, 5)
mapsywide_exons_internal_5ss_ext_tb <- as_tibble(mapsywide_exons_internal_5ss_ext)
mapsywide_exons_internal_5ss_ext_tb$seqs <- as.character(mapsywide_exons_internal_5ss_ext_seqs)
mapsywide_exons_internal_5ss_ext_tb$maxent <- as.numeric(mapsywide_exons_internal_5ss_ext_maxent)
write_tsv(mapsywide_exons_internal_5ss_ext_tb, gzfile("../../results/maxentscan_plots/mapsywide_exons_internal_5ss_ext_tb.txt"))

# visualize
ggplot(canonical_exons_internal_3ss_ext_tb) + theme_classic() + 
	geom_density(aes(maxent), color="black") + 
	geom_density(aes(maxent), color="red", data=mapsywide_exons_internal_3ss_ext_tb) + 
	geom_segment(aes(x=maxent_five_A, xend=maxent_five_A, y=0.19, yend=0.17), arrow = arrow(length = unit(0.2, "cm"), type="closed"), color=color_five_A) + 
	geom_segment(aes(x=maxent_five_B, xend=maxent_five_B, y=0.18, yend=0.16), arrow = arrow(length = unit(0.2, "cm"), type="closed"), color=color_five_B) + 
	xlab("MaxEntScan 3ss score") + ylab("Genome-wide density") + theme(aspect.ratio=0.9)
ggsave("../../results/maxentscan_plots/mapsywide_exons_internal_3ss_ext_tb.pdf", scale=0.3)
# ggsave("../../results/maxentscan_plots/mapsywide_exons_internal_3ss_ext_tb.png", scale=0.3)

ggplot(canonical_exons_internal_5ss_ext_tb) + theme_classic() + 
	geom_density(aes(maxent), color="black") + 
	geom_density(aes(maxent), color="red", data=mapsywide_exons_internal_5ss_ext_tb) + 
	geom_segment(aes(x=maxent_three_A, xend=maxent_three_A, y=0.26, yend=0.24), arrow = arrow(length = unit(0.2, "cm"), type="closed"), color=color_three_A) + 
	geom_segment(aes(x=maxent_three_B, xend=maxent_three_B, y=0.28, yend=0.26), arrow = arrow(length = unit(0.2, "cm"), type="closed"), color=color_three_B) + 
	xlab("MaxEntScan 5ss score") + ylab("Genome-wide density") + theme(aspect.ratio=0.9)
ggsave("../../results/maxentscan_plots/mapsywide_exons_internal_5ss_ext_tb.pdf", scale=0.3)
# ggsave("../../results/maxentscan_plots/mapsywide_exons_internal_5ss_ext_tb.png", scale=0.3)
