#!/bin/sh
# get ld r^2 statistics
for i in {1..22}
do
	sbatch linkage_filter_1000G_vcfs_filtered_chr${i}.sh
done
