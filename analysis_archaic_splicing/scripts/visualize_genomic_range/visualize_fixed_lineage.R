#!/bin/R
source("visualize_genomic_range.R")

# visualize
win_half_size_gtex=NA
win_half_size_splice=NA
win_half_size_seq=10
out_folder="../../results/visualize_genomic_range/"

vis_variant_id="8_24342791_G/A"
print(vis_variant_id)
visualize_genomic_range(vis_variant_id, win_half_size_gtex, win_half_size_splice, win_half_size_seq, out_folder)

vis_variant_id="X_65247366_T/C"
print(vis_variant_id)
visualize_genomic_range_no_LD(vis_variant_id, win_half_size_gtex, win_half_size_splice, win_half_size_seq, out_folder)

vis_variant_id="2_26702482_T/C"
print(vis_variant_id)
visualize_genomic_range(vis_variant_id, win_half_size_gtex, win_half_size_splice, win_half_size_seq, out_folder)

vis_variant_id="16_20435240_G/A"
print(vis_variant_id)
visualize_genomic_range(vis_variant_id, win_half_size_gtex, win_half_size_splice, win_half_size_seq, out_folder)

vis_variant_id="19_55598746_C/T"
print(vis_variant_id)
visualize_genomic_range(vis_variant_id, win_half_size_gtex, win_half_size_splice, win_half_size_seq, out_folder)

vis_variant_id="3_180381728_C/T"
print(vis_variant_id)
visualize_genomic_range(vis_variant_id, win_half_size_gtex, win_half_size_splice, win_half_size_seq, out_folder)

vis_variant_id="9_133576259_G/A"
print(vis_variant_id)
visualize_genomic_range(vis_variant_id, win_half_size_gtex, win_half_size_splice, win_half_size_seq, out_folder)

vis_variant_id="16_30390786_C/T"
print(vis_variant_id)
visualize_genomic_range(vis_variant_id, win_half_size_gtex, win_half_size_splice, win_half_size_seq, out_folder)

vis_variant_id="7_87470991_T/C"
print(vis_variant_id)
visualize_genomic_range(vis_variant_id, win_half_size_gtex, win_half_size_splice, win_half_size_seq, out_folder)
