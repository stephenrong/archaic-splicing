#!/bin/R
source("visualize_genomic_range.R")

# visualize
win_half_size_gtex=NA
win_half_size_splice=NA
win_half_size_seq=10
out_folder="../../results/visualize_genomic_range/"

vis_variant_id="4_38805942_G/C"
print(vis_variant_id)
visualize_genomic_range(vis_variant_id, win_half_size_gtex, win_half_size_splice, win_half_size_seq, out_folder)

vis_variant_id="1_22174518_G/T"
print(vis_variant_id)
visualize_genomic_range(vis_variant_id, win_half_size_gtex, win_half_size_splice, win_half_size_seq, out_folder)

vis_variant_id="1_46877284_T/G"
print(vis_variant_id)
visualize_genomic_range(vis_variant_id, win_half_size_gtex, win_half_size_splice, win_half_size_seq, out_folder)

vis_variant_id="12_21329761_G/A"
print(vis_variant_id)
visualize_genomic_range(vis_variant_id, win_half_size_gtex, win_half_size_splice, win_half_size_seq, out_folder)

vis_variant_id="1_169429035_C/T"
print(vis_variant_id)
visualize_genomic_range_manual_spliceai(vis_variant_id, win_half_size_gtex, win_half_size_splice, win_half_size_seq, out_folder, 0.59)
