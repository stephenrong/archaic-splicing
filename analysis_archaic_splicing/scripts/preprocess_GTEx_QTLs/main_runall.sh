#!/bin/sh

# Convert sQTLs to hub format, liftover to hg19, clean up
Rscript preprocess_GTEx_sQTLs.R

# Convert eQTLs to hub format, liftover to hg19, clean up
Rscript preprocess_GTEx_eQTLs.R

# Combine all tissue and QTL specific files
Rscript collate_GTEx_QTLs.R
