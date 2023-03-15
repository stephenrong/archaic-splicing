#!/bin/sh

# Add EUR_AC_bin, hapR2_EUR_ntag_bin
Rscript linkage_filter_1KGP_vcfs_append.R

# Take collated QTL files, add 1KGP SNP info
Rscript filter_in1KGP_GTEx_QTLs.R

# Filter final library SNP info, filter AI, filter to MAF >= 1% in EUR
Rscript filter_final_variants_SNPs.R

# Filter 1KGP SNP info, filter non-AI, filter to MAF >= 1% in EUR, create index based on AC and LD
Rscript filter_and_index_notAI.R

# Filter 1KGP SNP info, filter non-AI, filter to MAF >= 1% in EUR, filter to not low recomb, filter to accessible, create index based on AC and LD
Rscript filter_and_index_notAI_recmap_access.R
