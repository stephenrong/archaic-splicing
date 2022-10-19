#!/bin/R

library(tidyverse)
library(data.table)
library(ggpubr)

# load mapsy data
mapsy_table <- as_tibble(fread("../../../final-postmapsy-archaic-EndToEnd2/results/postprocess_Neanderthal_order/Neanderthal_order_variant_mpralm.txt.gz"))
mapsy_table_short <- mapsy_table %>% 
	# filter(filter_unique_seqs) %>% 
	# mutate(variant_exon_pair = paste(variant_id, exon_gene_name, exon_transcript_name, exon_exon_rank, sep="-")) %>% 
	mutate(construct_common_pair = paste(construct_type, Common_id_three, Common_id_five, sep="-")) %>% 
	dplyr::select(source, variant_id, exon_gene_name, exon_transcript_name, exon_exon_rank, construct_common_pair, mpralm.logFC)
mapsy_table_short_wider <- mapsy_table_short %>% 
	pivot_wider(names_from=construct_common_pair, values_from=mpralm.logFC)
write_tsv(mapsy_table_short_wider, gzfile("../../results/validate_half_exons/validate_half_exons_order/Neanderthal_order_mapsy_table_short.txt.gz"))

# visualize short half exons 3A and 3B
ggscatter(mapsy_table_short_wider, 
	x = "short_full_exon_variant_three-NA-NA", 
	y = "short_half_exon_variant_three-NA-Common_five_05", 
	alpha=0.5, color="#377eb8", 
	add = "reg.line",  # Add regression line
	add.params = list(color="black", fill="lightgray"), 
	# conf.int = TRUE, # Add confidence interval
	cor.coef = TRUE, # Add correlation coefficient
) + theme(aspect.ratio=1) + 
	xlab("Short full exon log2 FC") + 
	ylab("Short half exon 3A log2 FC") + 
	geom_vline(xintercept=0, color="lightgrey", size=0.33) + 
	geom_hline(yintercept=0, color="lightgrey", size=0.33) + 
	geom_vline(xintercept=0, color="lightgrey", size=0.33) + 
	geom_hline(yintercept=0, color="lightgrey")
ggsave("../../results/validate_half_exons/Neanderthal_order_short_full_short_half_common_five_05.pdf", scale=0.5)

ggscatter(mapsy_table_short_wider, 
	x = "short_full_exon_variant_three-NA-NA", 
	y = "short_half_exon_variant_three-NA-Common_five_06", 
	alpha=0.5, color="#377eb8", 
	add = "reg.line",  # Add regression line
	add.params = list(color="black", fill="lightgray"), 
	# conf.int = TRUE, # Add confidence interval
	cor.coef = TRUE, # Add correlation coefficient
) + theme(aspect.ratio=1) + 
	xlab("Short full exon log2 FC") + 
	ylab("Short half exon 3B log2 FC") + 
	geom_vline(xintercept=0, color="lightgrey", size=0.33) + 
	geom_hline(yintercept=0, color="lightgrey", size=0.33) + 
	geom_vline(xintercept=0, color="lightgrey", size=0.33) + 
	geom_hline(yintercept=0, color="lightgrey")
ggsave("../../results/validate_half_exons/Neanderthal_order_short_full_short_half_common_five_06.pdf", scale=0.5)

ggscatter(mapsy_table_short_wider, 
	x = "short_half_exon_variant_three-NA-Common_five_05", 
	y = "short_half_exon_variant_three-NA-Common_five_06", 
	alpha=0.5, color="#377eb8", 
	add = "reg.line",  # Add regression line
	add.params = list(color="black", fill="lightgray"), 
	# conf.int = TRUE, # Add confidence interval
	cor.coef = TRUE, # Add correlation coefficient
) + theme(aspect.ratio=1) + 
	xlab("Short half exon 3A log2 FC") + 
	ylab("Short half exon 3B log2 FC") + 
	geom_vline(xintercept=0, color="lightgrey", size=0.33) + 
	geom_hline(yintercept=0, color="lightgrey", size=0.33) + 
	geom_vline(xintercept=0, color="lightgrey", size=0.33) + 
	geom_hline(yintercept=0, color="lightgrey")
ggsave("../../results/validate_half_exons/Neanderthal_order_short_half_common_five_05_short_half_common_five_06.pdf", scale=0.5)

# visualize short half exons 5A and 5B
ggscatter(mapsy_table_short_wider, 
	x = "short_full_exon_variant_three-NA-NA", 
	y = "short_half_exon_variant_five_legacy-Common_three_05-NA", 
	alpha=0.5, color="#377eb8", 
	add = "reg.line",  # Add regression line
	add.params = list(color="black", fill="lightgray"), 
	# conf.int = TRUE, # Add confidence interval
	cor.coef = TRUE, # Add correlation coefficient
) + theme(aspect.ratio=1) + 
	xlab("Short full exon log2 FC") + 
	ylab("Short half exon 5A log2 FC") + 
	geom_vline(xintercept=0, color="lightgrey", size=0.33) + 
	geom_hline(yintercept=0, color="lightgrey", size=0.33) + 
	geom_vline(xintercept=0, color="lightgrey", size=0.33) + 
	geom_hline(yintercept=0, color="lightgrey")
ggsave("../../results/validate_half_exons/Neanderthal_order_short_full_short_half_common_three_05.pdf", scale=0.5)

ggscatter(mapsy_table_short_wider, 
	x = "short_full_exon_variant_three-NA-NA", 
	y = "short_half_exon_variant_five_legacy-Common_three_06-NA", 
	alpha=0.5, color="#377eb8", 
	add = "reg.line",  # Add regression line
	add.params = list(color="black", fill="lightgray"), 
	# conf.int = TRUE, # Add confidence interval
	cor.coef = TRUE, # Add correlation coefficient
) + theme(aspect.ratio=1) + 
	xlab("Short full exon log2 FC") + 
	ylab("Short half exon 5B log2 FC") + 
	geom_vline(xintercept=0, color="lightgrey", size=0.33) + 
	geom_hline(yintercept=0, color="lightgrey", size=0.33) + 
	geom_vline(xintercept=0, color="lightgrey", size=0.33) + 
	geom_hline(yintercept=0, color="lightgrey")
ggsave("../../results/validate_half_exons/Neanderthal_order_short_full_short_half_common_three_06.pdf", scale=0.5)

ggscatter(mapsy_table_short_wider, 
	x = "short_half_exon_variant_five_legacy-Common_three_05-NA", 
	y = "short_half_exon_variant_five_legacy-Common_three_06-NA", 
	alpha=0.5, color="#377eb8", 
	add = "reg.line",  # Add regression line
	add.params = list(color="black", fill="lightgray"), 
	# conf.int = TRUE, # Add confidence interval
	cor.coef = TRUE, # Add correlation coefficient
) + theme(aspect.ratio=1) + 
	xlab("Short half exon 5A log2 FC") + 
	ylab("Short half exon 5B log2 FC") + 
	geom_vline(xintercept=0, color="lightgrey", size=0.33) + 
	geom_hline(yintercept=0, color="lightgrey", size=0.33) + 
	geom_vline(xintercept=0, color="lightgrey", size=0.33) + 
	geom_hline(yintercept=0, color="lightgrey")
ggsave("../../results/validate_half_exons/Neanderthal_order_short_half_common_three_05_short_half_common_three_06.pdf", scale=0.5)

# visualize long half exons 3A and 3B
ggscatter(mapsy_table_short_wider, 
	x = "long_half_exon_variant_three-NA-Common_five_05", 
	y = "long_half_exon_variant_three-NA-Common_five_06", 
	alpha=0.5, color="#377eb8", 
	add = "reg.line",  # Add regression line
	add.params = list(color="black", fill="lightgray"), 
	# conf.int = TRUE, # Add confidence interval
	cor.coef = TRUE, # Add correlation coefficient
) + theme(aspect.ratio=1) + 
	xlab("Long half exon 3A log2 FC") + 
	ylab("Long half exon 3B log2 FC") + 
	geom_vline(xintercept=0, color="lightgrey", size=0.33) + 
	geom_hline(yintercept=0, color="lightgrey", size=0.33) + 
	geom_vline(xintercept=0, color="lightgrey", size=0.33) + 
	geom_hline(yintercept=0, color="lightgrey")
ggsave("../../results/validate_half_exons/Neanderthal_order_long_half_common_five_05_long_half_common_five_06.pdf", scale=0.5)

# visualize long half exons 5A and 5B
ggscatter(mapsy_table_short_wider, 
	x = "long_half_exon_variant_five-Common_three_05-NA", 
	y = "long_half_exon_variant_five-Common_three_06-NA", 
	alpha=0.5, color="#377eb8", 
	add = "reg.line",  # Add regression line
	add.params = list(color="black", fill="lightgray"), 
	# conf.int = TRUE, # Add confidence interval
	cor.coef = TRUE, # Add correlation coefficient
) + theme(aspect.ratio=1) + 
	xlab("Long half exon 5A log2 FC") + 
	ylab("Long half exon 5B log2 FC") + 
	geom_vline(xintercept=0, color="lightgrey", size=0.33) + 
	geom_hline(yintercept=0, color="lightgrey", size=0.33) + 
	geom_vline(xintercept=0, color="lightgrey", size=0.33) + 
	geom_hline(yintercept=0, color="lightgrey")
ggsave("../../results/validate_half_exons/Neanderthal_order_long_half_common_three_05_long_half_common_three_06.pdf", scale=0.5)
