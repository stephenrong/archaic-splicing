#!/bin/R

# Param file for generating variable part of oligos.

oligo_length <- 190  # oligo length minus primers
exon_length_max <- 120  # maximum length of full exons
exon_length_min_full <- 10  # minimum length of full exons
exon_length_min_half <- 40  # minimum length of half exons
intron_extend_three <- 55  # intronic portion near 3'ss
intron_extend_five <- 15  # intronic portion near 5'ss
exon_variant_half <- 90  # maximum length into half exons
exon_variant_buffer <- 10  # length of buffer zone
exon_variant_common <- 20  # length of common zone
intron_variant_buffer <- 5  # buffer for intron near 3'ss
# we use a smaller buffer at intron ends than within exons

# check for consistency
if (!(oligo_length == intron_extend_three + exon_length_max + intron_extend_five)) {
  stop ('oligo_length != intron_extend_three + exon_length_max + intron_extend_five')
}
if (!(exon_length_max == exon_variant_half + exon_variant_buffer + exon_variant_common)) {
  stop ('exon_length_max 1= exon_variant_half + exon_variant_buffer + exon_variant_common')
}
