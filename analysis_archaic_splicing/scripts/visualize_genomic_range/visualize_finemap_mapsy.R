#!/bin/R
source("visualize_genomic_range.R")

# visualize
win_half_size_gtex=NA
win_half_size_splice=NA
win_half_size_seq=10
out_folder="../../results/visualize_genomic_range/"

vis_variant_id="12_21329761_G/A"
print(vis_variant_id)
visualize_genomic_range(vis_variant_id, win_half_size_gtex, win_half_size_splice, win_half_size_seq, out_folder)

vis_variant_id="11_406483_T/C"
print(vis_variant_id)
visualize_genomic_range(vis_variant_id, win_half_size_gtex, win_half_size_splice, win_half_size_seq, out_folder)

vis_variant_id="12_133393323_C/T"
print(vis_variant_id)
visualize_genomic_range(vis_variant_id, win_half_size_gtex, win_half_size_splice, win_half_size_seq, out_folder)

vis_variant_id="2_183703336_G/A"
print(vis_variant_id)
visualize_genomic_range(vis_variant_id, win_half_size_gtex, win_half_size_splice, win_half_size_seq, out_folder)

vis_variant_id="5_172196752_A/G"
print(vis_variant_id)
visualize_genomic_range(vis_variant_id, win_half_size_gtex, win_half_size_splice, win_half_size_seq, out_folder)

vis_variant_id="16_336898_G/A"
print(vis_variant_id)
visualize_genomic_range(vis_variant_id, win_half_size_gtex, win_half_size_splice, win_half_size_seq, out_folder)

vis_variant_id="17_76121318_G/A"
print(vis_variant_id)
visualize_genomic_range(vis_variant_id, win_half_size_gtex, win_half_size_splice, win_half_size_seq, out_folder)
