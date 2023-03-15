#!/bin/sh

# SpliceAI annotate final variants
Rscript annotate_spliceai_final_part1.R
sh annotate_spliceai_final_part2.sh
Rscript annotate_spliceai_final_part3.R

# SpliceAI annotate 1KGP variants
Rscript annotate_spliceai_all1kgp_part1.R
sh annotate_spliceai_all1kgp_part2.sh
Rscript annotate_spliceai_all1kgp_part3_v2.R
