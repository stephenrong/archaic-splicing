#!/bin/R

library(tidyverse)
library(data.table)
library(plyranges)
library(rtracklayer)

# load
altai_denisovan <- import("../../../final-splicing-archaic-variants/data/additional_archaic_vcfs/altai_denisovan/chr_mask.bed.gz", sep="\t")
altai_neanderthal <- import("../../../final-splicing-archaic-variants/data/additional_archaic_vcfs/altai_neanderthal/chr_mask.bed.gz", sep="\t")
chagyrskaya_neanderthal <- import("../../../final-splicing-archaic-variants/data/additional_archaic_vcfs/chagyrskaya_neanderthal/chr_mask.bed.gz", sep="\t")
vindija_neanderthal <- import("../../../final-splicing-archaic-variants/data/additional_archaic_vcfs/vindija_neanderthal/chr_mask.bed.gz", sep="\t")

# intersect
archaic_mask_inter <- altai_denisovan %>% intersect(altai_neanderthal) %>% intersect(chagyrskaya_neanderthal) %>% intersect(vindija_neanderthal)
write_tsv(as_tibble(archaic_mask_inter), gzfile("../../results/preprocess_1KGP_SNPs/archaic_mask_inter.bed.txt.gz"))

# union
archaic_mask_union <- altai_denisovan %>% union(altai_neanderthal) %>% union(chagyrskaya_neanderthal) %>% union(vindija_neanderthal)
write_tsv(as_tibble(archaic_mask_union), gzfile("../../results/preprocess_1KGP_SNPs/archaic_mask_union.bed.txt.gz"))
