#!/bin/sh
# tabix -h ../../data/annotate_spliceai/spliceai_scores.raw.snv.hg19.vcf.gz -R ../../results/annotate_splice_prediction/ALL_1KGP_phase3_hub.bed | bgzip >| ../../results/annotate_splice_prediction/ALL_1KGP_phase3_spliceai_gs_raw.vcf.gz

for CHR in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X
do 
	sbatch --wrap "tabix -h ../../data/annotate_spliceai/spliceai_scores.raw.snv.hg19.vcf.gz -R ../../results/annotate_splice_prediction/ALL_1KGP_phase3_hub_chr${CHR}.bed | bgzip >| ../../results/annotate_splice_prediction/ALL_1KGP_phase3_spliceai_gs_raw_chr${CHR}.vcf.gz"
done

for CHR in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X
do 
	tabix -p vcf ../../results/annotate_splice_prediction/ALL_1KGP_phase3_spliceai_gs_raw_chr${CHR}.vcf.gz
done

bcftools merge \
	../../results/annotate_splice_prediction/ALL_1KGP_phase3_spliceai_gs_raw_chr1.vcf.gz \
	../../results/annotate_splice_prediction/ALL_1KGP_phase3_spliceai_gs_raw_chr2.vcf.gz \
	../../results/annotate_splice_prediction/ALL_1KGP_phase3_spliceai_gs_raw_chr3.vcf.gz \
	../../results/annotate_splice_prediction/ALL_1KGP_phase3_spliceai_gs_raw_chr4.vcf.gz \
	../../results/annotate_splice_prediction/ALL_1KGP_phase3_spliceai_gs_raw_chr5.vcf.gz \
	../../results/annotate_splice_prediction/ALL_1KGP_phase3_spliceai_gs_raw_chr6.vcf.gz \
	../../results/annotate_splice_prediction/ALL_1KGP_phase3_spliceai_gs_raw_chr7.vcf.gz \
	../../results/annotate_splice_prediction/ALL_1KGP_phase3_spliceai_gs_raw_chr8.vcf.gz \
	../../results/annotate_splice_prediction/ALL_1KGP_phase3_spliceai_gs_raw_chr9.vcf.gz \
	../../results/annotate_splice_prediction/ALL_1KGP_phase3_spliceai_gs_raw_chr10.vcf.gz \
	../../results/annotate_splice_prediction/ALL_1KGP_phase3_spliceai_gs_raw_chr11.vcf.gz \
	../../results/annotate_splice_prediction/ALL_1KGP_phase3_spliceai_gs_raw_chr12.vcf.gz \
	../../results/annotate_splice_prediction/ALL_1KGP_phase3_spliceai_gs_raw_chr13.vcf.gz \
	../../results/annotate_splice_prediction/ALL_1KGP_phase3_spliceai_gs_raw_chr14.vcf.gz \
	../../results/annotate_splice_prediction/ALL_1KGP_phase3_spliceai_gs_raw_chr15.vcf.gz \
	../../results/annotate_splice_prediction/ALL_1KGP_phase3_spliceai_gs_raw_chr16.vcf.gz \
	../../results/annotate_splice_prediction/ALL_1KGP_phase3_spliceai_gs_raw_chr17.vcf.gz \
	../../results/annotate_splice_prediction/ALL_1KGP_phase3_spliceai_gs_raw_chr18.vcf.gz \
	../../results/annotate_splice_prediction/ALL_1KGP_phase3_spliceai_gs_raw_chr19.vcf.gz \
	../../results/annotate_splice_prediction/ALL_1KGP_phase3_spliceai_gs_raw_chr20.vcf.gz \
	../../results/annotate_splice_prediction/ALL_1KGP_phase3_spliceai_gs_raw_chr21.vcf.gz \
	../../results/annotate_splice_prediction/ALL_1KGP_phase3_spliceai_gs_raw_chr22.vcf.gz \
	../../results/annotate_splice_prediction/ALL_1KGP_phase3_spliceai_gs_raw_chrX.vcf.gz \
	-o ../../results/annotate_splice_prediction/ALL_1KGP_phase3_spliceai_gs_raw.vcf.gz -O z
