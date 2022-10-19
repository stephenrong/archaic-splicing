#!/bin/sh

# overlap UK Biobank
Rscript get_overlap_variants_ukb.R

# ovelap Biobank Japan
Rscript preprocess_pip_bbj.R
Rscript get_overlap_variants_bbj.R
