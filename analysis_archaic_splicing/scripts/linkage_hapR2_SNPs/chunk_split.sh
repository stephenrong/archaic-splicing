#!/bin/sh

zcat ALL.chr1.phase3_shapeit2_mvncall_integrated_v5a.20130502.hapR2.hap.ld.gz | split -l 500000000 - ALL.chr1.phase3_shapeit2_mvncall_integrated_v5a.20130502.hapR2.hap.ld.gz.part-
parallel bgzip {} ::: bgzip ALL.chr1.phase3_shapeit2_mvncall_integrated_v5a.20130502.hapR2.hap.ld.gz.part-*
zcat ALL.chr2.phase3_shapeit2_mvncall_integrated_v5a.20130502.hapR2.hap.ld.gz | split -l 500000000 - ALL.chr2.phase3_shapeit2_mvncall_integrated_v5a.20130502.hapR2.hap.ld.gz.part-
parallel bgzip {} ::: bgzip ALL.chr2.phase3_shapeit2_mvncall_integrated_v5a.20130502.hapR2.hap.ld.gz.part-*
zcat ALL.chr3.phase3_shapeit2_mvncall_integrated_v5a.20130502.hapR2.hap.ld.gz | split -l 500000000 - ALL.chr3.phase3_shapeit2_mvncall_integrated_v5a.20130502.hapR2.hap.ld.gz.part-
parallel bgzip {} ::: bgzip ALL.chr3.phase3_shapeit2_mvncall_integrated_v5a.20130502.hapR2.hap.ld.gz.part-*
zcat ALL.chr4.phase3_shapeit2_mvncall_integrated_v5a.20130502.hapR2.hap.ld.gz | split -l 500000000 - ALL.chr4.phase3_shapeit2_mvncall_integrated_v5a.20130502.hapR2.hap.ld.gz.part-
parallel bgzip {} ::: bgzip ALL.chr4.phase3_shapeit2_mvncall_integrated_v5a.20130502.hapR2.hap.ld.gz.part-*
zcat ALL.chr5.phase3_shapeit2_mvncall_integrated_v5a.20130502.hapR2.hap.ld.gz | split -l 500000000 - ALL.chr5.phase3_shapeit2_mvncall_integrated_v5a.20130502.hapR2.hap.ld.gz.part-
parallel bgzip {} ::: bgzip ALL.chr5.phase3_shapeit2_mvncall_integrated_v5a.20130502.hapR2.hap.ld.gz.part-*
zcat ALL.chr6.phase3_shapeit2_mvncall_integrated_v5a.20130502.hapR2.hap.ld.gz | split -l 500000000 - ALL.chr6.phase3_shapeit2_mvncall_integrated_v5a.20130502.hapR2.hap.ld.gz.part-
parallel bgzip {} ::: bgzip ALL.chr6.phase3_shapeit2_mvncall_integrated_v5a.20130502.hapR2.hap.ld.gz.part-*
zcat ALL.chr7.phase3_shapeit2_mvncall_integrated_v5a.20130502.hapR2.hap.ld.gz | split -l 500000000 - ALL.chr7.phase3_shapeit2_mvncall_integrated_v5a.20130502.hapR2.hap.ld.gz.part-
parallel bgzip {} ::: bgzip ALL.chr7.phase3_shapeit2_mvncall_integrated_v5a.20130502.hapR2.hap.ld.gz.part-*
zcat ALL.chr8.phase3_shapeit2_mvncall_integrated_v5a.20130502.hapR2.hap.ld.gz | split -l 500000000 - ALL.chr8.phase3_shapeit2_mvncall_integrated_v5a.20130502.hapR2.hap.ld.gz.part-
parallel bgzip {} ::: bgzip ALL.chr8.phase3_shapeit2_mvncall_integrated_v5a.20130502.hapR2.hap.ld.gz.part-*
zcat ALL.chr9.phase3_shapeit2_mvncall_integrated_v5a.20130502.hapR2.hap.ld.gz | split -l 500000000 - ALL.chr9.phase3_shapeit2_mvncall_integrated_v5a.20130502.hapR2.hap.ld.gz.part-
parallel bgzip {} ::: bgzip ALL.chr9.phase3_shapeit2_mvncall_integrated_v5a.20130502.hapR2.hap.ld.gz.part-*
zcat ALL.chr10.phase3_shapeit2_mvncall_integrated_v5a.20130502.hapR2.hap.ld.gz | split -l 500000000 - ALL.chr10.phase3_shapeit2_mvncall_integrated_v5a.20130502.hapR2.hap.ld.gz.part-
parallel bgzip {} ::: bgzip ALL.chr10.phase3_shapeit2_mvncall_integrated_v5a.20130502.hapR2.hap.ld.gz.part-*
zcat ALL.chr11.phase3_shapeit2_mvncall_integrated_v5a.20130502.hapR2.hap.ld.gz | split -l 500000000 - ALL.chr11.phase3_shapeit2_mvncall_integrated_v5a.20130502.hapR2.hap.ld.gz.part-
parallel bgzip {} ::: bgzip ALL.chr11.phase3_shapeit2_mvncall_integrated_v5a.20130502.hapR2.hap.ld.gz.part-*
zcat ALL.chr12.phase3_shapeit2_mvncall_integrated_v5a.20130502.hapR2.hap.ld.gz | split -l 500000000 - ALL.chr12.phase3_shapeit2_mvncall_integrated_v5a.20130502.hapR2.hap.ld.gz.part-
parallel bgzip {} ::: bgzip ALL.chr12.phase3_shapeit2_mvncall_integrated_v5a.20130502.hapR2.hap.ld.gz.part-*
zcat ALL.chr13.phase3_shapeit2_mvncall_integrated_v5a.20130502.hapR2.hap.ld.gz | split -l 500000000 - ALL.chr13.phase3_shapeit2_mvncall_integrated_v5a.20130502.hapR2.hap.ld.gz.part-
parallel bgzip {} ::: bgzip ALL.chr13.phase3_shapeit2_mvncall_integrated_v5a.20130502.hapR2.hap.ld.gz.part-*
zcat ALL.chr14.phase3_shapeit2_mvncall_integrated_v5a.20130502.hapR2.hap.ld.gz | split -l 500000000 - ALL.chr14.phase3_shapeit2_mvncall_integrated_v5a.20130502.hapR2.hap.ld.gz.part-
parallel bgzip {} ::: bgzip ALL.chr14.phase3_shapeit2_mvncall_integrated_v5a.20130502.hapR2.hap.ld.gz.part-*
zcat ALL.chr15.phase3_shapeit2_mvncall_integrated_v5a.20130502.hapR2.hap.ld.gz | split -l 500000000 - ALL.chr15.phase3_shapeit2_mvncall_integrated_v5a.20130502.hapR2.hap.ld.gz.part-
parallel bgzip {} ::: bgzip ALL.chr15.phase3_shapeit2_mvncall_integrated_v5a.20130502.hapR2.hap.ld.gz.part-*
zcat ALL.chr16.phase3_shapeit2_mvncall_integrated_v5a.20130502.hapR2.hap.ld.gz | split -l 500000000 - ALL.chr16.phase3_shapeit2_mvncall_integrated_v5a.20130502.hapR2.hap.ld.gz.part-
parallel bgzip {} ::: bgzip ALL.chr16.phase3_shapeit2_mvncall_integrated_v5a.20130502.hapR2.hap.ld.gz.part-*
zcat ALL.chr17.phase3_shapeit2_mvncall_integrated_v5a.20130502.hapR2.hap.ld.gz | split -l 500000000 - ALL.chr17.phase3_shapeit2_mvncall_integrated_v5a.20130502.hapR2.hap.ld.gz.part-
parallel bgzip {} ::: bgzip ALL.chr17.phase3_shapeit2_mvncall_integrated_v5a.20130502.hapR2.hap.ld.gz.part-*
zcat ALL.chr18.phase3_shapeit2_mvncall_integrated_v5a.20130502.hapR2.hap.ld.gz | split -l 500000000 - ALL.chr18.phase3_shapeit2_mvncall_integrated_v5a.20130502.hapR2.hap.ld.gz.part-
parallel bgzip {} ::: bgzip ALL.chr18.phase3_shapeit2_mvncall_integrated_v5a.20130502.hapR2.hap.ld.gz.part-*
zcat ALL.chr19.phase3_shapeit2_mvncall_integrated_v5a.20130502.hapR2.hap.ld.gz | split -l 500000000 - ALL.chr19.phase3_shapeit2_mvncall_integrated_v5a.20130502.hapR2.hap.ld.gz.part-
parallel bgzip {} ::: bgzip ALL.chr19.phase3_shapeit2_mvncall_integrated_v5a.20130502.hapR2.hap.ld.gz.part-*
zcat ALL.chr20.phase3_shapeit2_mvncall_integrated_v5a.20130502.hapR2.hap.ld.gz | split -l 500000000 - ALL.chr20.phase3_shapeit2_mvncall_integrated_v5a.20130502.hapR2.hap.ld.gz.part-
parallel bgzip {} ::: bgzip ALL.chr20.phase3_shapeit2_mvncall_integrated_v5a.20130502.hapR2.hap.ld.gz.part-*
zcat ALL.chr21.phase3_shapeit2_mvncall_integrated_v5a.20130502.hapR2.hap.ld.gz | split -l 500000000 - ALL.chr21.phase3_shapeit2_mvncall_integrated_v5a.20130502.hapR2.hap.ld.gz.part-
parallel bgzip {} ::: bgzip ALL.chr21.phase3_shapeit2_mvncall_integrated_v5a.20130502.hapR2.hap.ld.gz.part-*
zcat ALL.chr22.phase3_shapeit2_mvncall_integrated_v5a.20130502.hapR2.hap.ld.gz | split -l 500000000 - ALL.chr22.phase3_shapeit2_mvncall_integrated_v5a.20130502.hapR2.hap.ld.gz.part-
parallel bgzip {} ::: bgzip ALL.chr22.phase3_shapeit2_mvncall_integrated_v5a.20130502.hapR2.hap.ld.gz.part-*
