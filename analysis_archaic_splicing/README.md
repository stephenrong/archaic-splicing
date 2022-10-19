# Neanderthal MaPSy Manuscript

Updated: 2022-10-18. Contact: Stephen Rong.


### Introduction

This repository contains the scripts for reproducing the analyses, tables, and figures in the manuscript.


### Contents

**scripts/**: Contains the scripts for reproducing analyses, tables, and figures in manuscript.
  * additional_analyses_1KGP_SpliceAI_diff
  * additional_analyses_constrained_genes
  * additional_analyses_enrichment_vep_prediction
  * additional_analyses_fine_map_overlap
  * annotate_splice_prediction
  * distribution_splice_prediction
  * enrichment_GTEx_QTLs_v2
  * enrichment_splice_prediction_gw_raw
  * linkage_hapR2_SNPs
  * mapsy_to_variant_table_updated
  * preprocess_1KGP_SNPs
  * preprocess_GTEx_QTLs
  * supplementary_tables
  * upset_plot_variants
  * validate_half_exons
  * vep_annotations_SNPs
  * visualize_genomic_range

**results/**: Contains the result corresponding to scripts subdirectory of same name.
  * additional_analyses_1KGP_SpliceAI_diff
  * additional_analyses_constrained_genes
  * additional_analyses_enrichment_vep_prediction
  * additional_analyses_fine_map_overlap
  * annotate_splice_prediction: 
      large files available at (https://doi.org/10.5281/zenodo.7158564).
  * distribution_splice_prediction
  * enrichment_GTEx_QTLs_v2: 
      large files available at (https://doi.org/10.5281/zenodo.7158564).
  * enrichment_splice_prediction_gw_raw
  * linkage_hapR2_SNPs
  * mapsy_to_variant_table_updated
  * preprocess_1KGP_SNPs: 
      large files available at (https://doi.org/10.5281/zenodo.7158564).
  * preprocess_GTEx_QTLs
  * supplementary_tables
  * upset_plot_variants
  * validate_half_exons
  * vep_annotations_SNPs: 
      large files available at (https://doi.org/10.5281/zenodo.7158564).
  * visualize_genomic_range

**data/**: Contains the input data needed for running the scripts.
  * annotate_B_statistics: 
  * annotate_GTEx_eQTLs
  * annotate_GTEx_sQTLs
  * annotate_hexamer_scores
  * annotate_lift_over
  * annotate_spliceai
  * fasta
  * finemap_overlap
  * gnomAD_v2_constraint
  * premapsy_variants

### Dependencies

R (version 4.0.2), with packages: tidyverse (1.3.1), data.table (1.14.2), plyranges (1.8.0), wrapr (2.0.9), vcfR (1.12.0), rtracklayer (1.48.0), EnsDb.Hsapiens.v86 (2.99.0), ggpubr (0.4.0), Hmisc (4.6-0), rstatix (0.7.0), MetBrewer (0.2.0), BSgenome.Hsapiens.UCSC.hg19 (1.4.3), BSgenome.Hsapiens.UCSC.hg38 (1.4.3), UpSetR (1.4.0), optparse (1.7.1), gtools (3.9.2), ggrepel (0.9.1), ggbio (1.36.0), viridis (0.6.2), Gviz (1.32.0), erma (1.4.0). Other software: ensembl-vep (107), htslib (1.10.2), CrossMap (0.5.4).