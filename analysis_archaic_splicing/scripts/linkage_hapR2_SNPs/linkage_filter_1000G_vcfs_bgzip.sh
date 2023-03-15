#!/bin/sh
# bgzip filtered 1000G VCFs
for i in {1..22}
do
	sbatch linkage_filter_1000G_vcfs_bgzip_chr${i}.sh
done
