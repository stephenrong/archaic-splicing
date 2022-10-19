#!/bin/R

library(tidyverse)
library(data.table)
library(ggthemes)
library(MetBrewer)

# load files
region_enrichments_list <- list()
region_enrichments_list[["intergenic"]] <- as_tibble(fread("../../results/enrichment_GTEx_QTLs_v2/archaic_introgressed/archaic_introgressed_matched_EUR_AC_binxEUR_hapR2tag_bin_control_standard_VEP_consequence_intergenic_variant-enrichment_score_table.txt"))
region_enrichments_list[["downstream"]] <- as_tibble(fread("../../results/enrichment_GTEx_QTLs_v2/archaic_introgressed/archaic_introgressed_matched_EUR_AC_binxEUR_hapR2tag_bin_control_standard_VEP_consequence_downstream_gene_variant-enrichment_score_table.txt"))
region_enrichments_list[["upstream"]] <- as_tibble(fread("../../results/enrichment_GTEx_QTLs_v2/archaic_introgressed/archaic_introgressed_matched_EUR_AC_binxEUR_hapR2tag_bin_control_standard_VEP_consequence_upstream_gene_variant-enrichment_score_table.txt"))
region_enrichments_list[["3 prime UTR"]] <- as_tibble(fread("../../results/enrichment_GTEx_QTLs_v2/archaic_introgressed/archaic_introgressed_matched_EUR_AC_binxEUR_hapR2tag_bin_control_standard_VEP_consequence_3_prime_UTR_variant-enrichment_score_table.txt"))
region_enrichments_list[["5 prime UTR"]] <- as_tibble(fread("../../results/enrichment_GTEx_QTLs_v2/archaic_introgressed/archaic_introgressed_matched_EUR_AC_binxEUR_hapR2tag_bin_control_standard_VEP_consequence_5_prime_UTR_variant-enrichment_score_table.txt"))
region_enrichments_list[["intronic"]] <- as_tibble(fread("../../results/enrichment_GTEx_QTLs_v2/archaic_introgressed/archaic_introgressed_matched_EUR_AC_binxEUR_hapR2tag_bin_control_standard_VEP_consequence_intron_variant-enrichment_score_table.txt"))
region_enrichments_list[["synonymous"]] <- as_tibble(fread("../../results/enrichment_GTEx_QTLs_v2/archaic_introgressed/archaic_introgressed_matched_EUR_AC_binxEUR_hapR2tag_bin_control_standard_VEP_consequence_synonymous_variant-enrichment_score_table.txt"))
region_enrichments_list[["missense"]] <- as_tibble(fread("../../results/enrichment_GTEx_QTLs_v2/archaic_introgressed/archaic_introgressed_matched_EUR_AC_binxEUR_hapR2tag_bin_control_standard_VEP_consequence_missense_variant-enrichment_score_table.txt"))

# join files
region_enrichments_final <- bind_rows(region_enrichments_list, .id="region")
region_enrichments_final$region <- factor(region_enrichments_final$region, levels=c("intergenic", "downstream", "upstream", "3 prime UTR", "5 prime UTR", "intronic", "synonymous", "missense"))  # , "eQTL", "sQTL", "3'aQTL"))
write_tsv(region_enrichments_final, gzfile("../../results/enrichment_GTEx_QTLs_v2/figures/archaic_introgressed-region-enrichment_scores_table_pivot_hapR2_pointerror.txt.gz"))

# visualize
ggplot(region_enrichments_final) + 
	geom_hline(yintercept=0, linetype="dashed", color="grey70") + 
	geom_boxplot(aes(region, value, color=region)) + xlab("Variant effect") + ylab("log2(odds-ratio)") + 
	theme_base() + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + theme(legend.position = "none") + theme(aspect.ratio=0.3) + theme(plot.background=element_blank()) + scale_color_met_d("Monet")
ggsave("../../results/enrichment_GTEx_QTLs_v2/figures/archaic_introgressed-region-enrichment_scores_table_pivot_hapR2_pointerror.pdf", device=cairo_pdf, scale=0.7)
