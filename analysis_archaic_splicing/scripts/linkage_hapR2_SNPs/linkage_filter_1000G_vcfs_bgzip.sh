#!/bin/sh
# bgzip filtered 1000G VCFs
for i in {1..22}
do
	bgzip ../../results/linkage_hapR2_SNPs/ALL.chr${i}.phase3_shapeit2_mvncall_integrated_v5a.20130502.hapR2.hap.ld >| ../../results/linkage_hapR2_SNPs/ALL.chr${i}.phase3_shapeit2_mvncall_integrated_v5a.20130502.hapR2.hap.ld.gz
done
