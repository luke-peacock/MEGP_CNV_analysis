#!/bin/bash

chromosomes="
chr1
chr2
chr3
chr4
chr5
chr6
chr7
chr8
chr9
chr10
chr11
chr12
chr13
chr14
chr15
chr16
chr17
chr18
chr19
chr20
chr21
chr22
"

for chrom in $chromosomes; do

	/homes/athosnew/Genetics_Centre_Bioinformatics/resourses/vcftools/bin/vcftools --gzvcf /apittman-aramis/sgul/shares/aramis/Alan/SGUL_Exomes_R3_April_2021/Sample_characterisation/GRC_v3_16-3-2021.joint-called.vqsr_split_leftAligned.annotated.PASS.vcf.gz --remove-filtered-all --LROH --out /apittman-aramis/sgul/shares/aramis/Genetics_Centre_Bioinformatics/CNV_SGUL_Exomes/ROH_calls/vcftools/vcftools_roh_${chrom} --chr ${chrom}
	
done
