# Neanderthal MaPSy Manuscript

Updated: 2023-03-15. Contact: Stephen Rong.


### Introduction

This repository contains the scripts for reproducing the analyses, tables, and figures in the manuscript.


### Contents

**scripts/**: Contains the scripts for reproducing analyses, tables, and figures in manuscript.
  * additional_analyses_1KGP_SpliceAI_diff: SpliceAI differences between 1KGP populations.
  * additional_analyses_constrained_genes: MaPSy and SpliceAI vs gnomAD constrained genes.
  * additional_analyses_enrichment_vep_prediction: Enrichments of Ensembl VEP proportions.
  * additional_analyses_fine_map_overlap: UKB and BBJ fine-mapped variant overlaps.
  * annotate_splice_prediction: SpliceAI annotations.
  * distribution_splice_prediction: MaPSy vs Ensembl VEP, SpliceAI, and hexamer ESE/ESS scores.
  * enrichment_GTEx_QTLs_v2: Enrichments of GTEx, VEP, and SpliceAI relative to matched controls.
  * enrichment_splice_prediction_gw_raw: Enrichment of SpliceAI proportions.
  * linkage_hapR2_SNPs: LD statistics used for enrichment analyses.
  * mapsy_to_variant_table_updated: Main comparison of proportion figures.
  * maxentscan_plots: MaxEntScan plots for half exon MaPSy construct validations.
  * preprocess_1KGP_SNPs: Generate final hominin evolution variant sets.
  * preprocess_GTEx_QTLs: Collate GTEx eVariants and sVariants.
  * supplementary_tables: Generate supplementary tables.
  * upset_plot_variants: UpSet plots of variant sets.
  * validate_half_exons: Half exon MaPSy construct validations.
  * vep_annotations_SNPs: Ensembl VEP annotations.
  * visualize_genomic_range: Visualizations of individual splicing variants.
  * get_helper.R: Functions for standardizing variant tables.

**results/**: Contains the result corresponding to scripts subdirectory of same name.
  * additional_analyses_1KGP_SpliceAI_diff
  * additional_analyses_constrained_genes
  * additional_analyses_enrichment_vep_prediction
  * additional_analyses_fine_map_overlap
  * annotate_splice_prediction: large files available at (https://doi.org/10.5281/zenodo.7158564).
  * distribution_splice_prediction
  * enrichment_GTEx_QTLs_v2: large files available at (https://doi.org/10.5281/zenodo.7158564).
  * enrichment_splice_prediction_gw_raw
  * linkage_hapR2_SNPs
  * mapsy_to_variant_table_updated
  * maxentscan_plots
  * preprocess_1KGP_SNPs: large files available at (https://doi.org/10.5281/zenodo.7158564).
  * preprocess_GTEx_QTLs
  * supplementary_tables
  * upset_plot_variants
  * validate_half_exons
  * vep_annotations_SNPs: large files available at (https://doi.org/10.5281/zenodo.7158564).
  * visualize_genomic_range

**data/**: Contains the input data needed for running the scripts.
  * annotate_B_statistics: B statistics for background selection from McVicker et al. (2009).
  * annotate_GTEx_eQTLs: GTEx v8 eQTL significant variant gene pairs.
  * annotate_GTEx_sQTLs: GTEx v8 sQTL significant variant phenotype pairs
  * annotate_hexamer_scores: Hexamer ESE/ESS scores from Ke et al. (2011) and Rosenberg et al. (2015).
  * annotate_lift_over: UCSC liftOver chains for hg18, hg19, and hg38.
  * annotate_spliceai: SpliceAI raw genome scores v1.3.
  * fasta: Homo sapiens GRCh37 reference genome.
  * finemap_overlap: UKB and BBJ fine-mapped variants from Kanai et al. (2021)
  * gnomAD_v2_constraint: gnomAD LOEUF metrics from Karczewski et al. (2020.)
  * premapsy_variants: MaPSy variants from an earlier analysis.

### Dependencies

R (version 4.0.2), with packages: tidyverse (1.3.1), data.table (1.14.2), plyranges (1.8.0), wrapr (2.0.9), vcfR (1.12.0), rtracklayer (1.48.0), EnsDb.Hsapiens.v86 (2.99.0), ggpubr (0.4.0), Hmisc (4.6-0), rstatix (0.7.0), MetBrewer (0.2.0), BSgenome.Hsapiens.UCSC.hg19 (1.4.3), BSgenome.Hsapiens.UCSC.hg38 (1.4.3), UpSetR (1.4.0), optparse (1.7.1), gtools (3.9.2), ggrepel (0.9.1), ggbio (1.36.0), viridis (0.6.2), Gviz (1.32.0), erma (1.4.0). Other software: ensembl-vep (107), htslib (1.10.2), CrossMap (0.5.4).