#!/bin/R
source("visualize_genomic_range.R")

# visualize
win_half_size_gtex=NA
win_half_size_splice=NA
win_half_size_seq=10
out_folder="../../results/visualize_genomic_range/"

vis_variant_id="12_52946336_A/G"
print(vis_variant_id)
visualize_genomic_range_mapsyless(vis_variant_id, win_half_size_gtex, win_half_size_splice, win_half_size_seq, out_folder)

vis_variant_id="14_92791244_T/A"
print(vis_variant_id)
visualize_genomic_range_mapsyless(vis_variant_id, win_half_size_gtex, win_half_size_splice, win_half_size_seq, out_folder)

vis_variant_id="18_55238820_A/G"
print(vis_variant_id)
visualize_genomic_range_mapsyless(vis_variant_id, win_half_size_gtex, win_half_size_splice, win_half_size_seq, out_folder)

vis_variant_id="3_195623905_A/C"
print(vis_variant_id)
visualize_genomic_range_mapsyless(vis_variant_id, win_half_size_gtex, win_half_size_splice, win_half_size_seq, out_folder)

vis_variant_id="4_169348039_C/T"
print(vis_variant_id)
visualize_genomic_range_mapsyless(vis_variant_id, win_half_size_gtex, win_half_size_splice, win_half_size_seq, out_folder)

vis_variant_id="4_169348476_G/A"
print(vis_variant_id)
visualize_genomic_range_mapsyless(vis_variant_id, win_half_size_gtex, win_half_size_splice, win_half_size_seq, out_folder)

vis_variant_id="8_22020259_G/A"
print(vis_variant_id)
visualize_genomic_range_mapsyless(vis_variant_id, win_half_size_gtex, win_half_size_splice, win_half_size_seq, out_folder)

vis_variant_id="9_94119997_C/T"
print(vis_variant_id)
visualize_genomic_range_mapsyless(vis_variant_id, win_half_size_gtex, win_half_size_splice, win_half_size_seq, out_folder)
