#!/bin/sh
sh linkage_filter_1000G_vcfs.sh
sh linkage_filter_1000G_vcfs_filtered.sh
sh linkage_filter_1000G_vcfs_hapR2.sh
sh linkage_filter_1000G_vcfs_bgzip.sh
Rscript linkage_filter_1000G_vcfs_postprocess.R