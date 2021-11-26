#!/bin/bash

# copy files from originals

cp /apittman-aramis/sgul/shares/aramis/Genetics_Centre_Bioinformatics/CNV_SGUL_Exomes/bedtools/CNV_ALL_OMIM_annotated_point01.csv /apittman-aramis/sgul/shares/aramis/Genetics_Centre_Bioinformatics/CNV_SGUL_Exomes/bedtools/CNV_ALL_OMIM_annotated_point01_header.csv
cp /apittman-aramis/sgul/shares/aramis/Genetics_Centre_Bioinformatics/CNV_SGUL_Exomes/bedtools/CNV_ROH_OMIM_annotated_point01.csv /apittman-aramis/sgul/shares/aramis/Genetics_Centre_Bioinformatics/CNV_SGUL_Exomes/bedtools/CNV_ROH_OMIM_annotated_point01_header.csv

# put headers on files

sed -i '1s/^/chrom_OMIM,OMIM_start,OMIM_end,OMIM_gene_name,OMIM_phenotype,CNV_chromosome,CNV_start,CNV_end,batch_ID,sample_ID,start.p,end.p,CNV_type,nexons,CNV_ID,BF,reads.expected,reads.observed,reads.ratio,Conrad.hg19,exons.hg19,OMIM_bp_overlap\n/' /apittman-aramis/sgul/shares/aramis/Genetics_Centre_Bioinformatics/CNV_SGUL_Exomes/bedtools/CNV_ALL_OMIM_annotated_point01_header.csv

sed -i '1s/^/chrom_OMIM,OMIM_start,OMIM_end,OMIM_gene_name,OMIM_phenotype,CNV_chromosome,CNV_start,CNV_end,ROH_chromosome,ROH_start,ROH_end,ROH_sample_ID,ROH_length,ROH_number_of_markers,ROH_quality_score,CNV_ID,batch_ID,sample_ID,start.p,end.p,CNV_type,nexons,BF,reads.expected,reads.observed,reads.ratio,Conrad.hg19,exons.hg19,roh_bp_overlap,OMIM_bp_overlap\n/' /apittman-aramis/sgul/shares/aramis/Genetics_Centre_Bioinformatics/CNV_SGUL_Exomes/bedtools/CNV_ROH_OMIM_annotated_point01_header.csv

############### ALL ####################
# for loop all CNV IDs to remove from CNV_ALL

cat /apittman-aramis/sgul/shares/aramis/Genetics_Centre_Bioinformatics/CNV_SGUL_Exomes/bedtools/CNV_ALL_OMIM_annotated_point01_header.csv | cut -d ',' -f15 | tail -n +2 > /apittman-aramis/sgul/shares/aramis/Genetics_Centre_Bioinformatics/CNV_SGUL_Exomes/bedtools/IDs_ALL_OMIM.txt

cp /apittman-aramis/sgul/shares/aramis/Genetics_Centre_Bioinformatics/CNV_SGUL_Exomes/SQL/combined_all_csv_full_ID.csv /apittman-aramis/sgul/shares/aramis/Genetics_Centre_Bioinformatics/CNV_SGUL_Exomes/SQL/combined_all_csv_full_ID_deldups.csv

# delete from all_CSVs file any rows containing IDs found in the OMIM_all_CNVs
grep -vFf /apittman-aramis/sgul/shares/aramis/Genetics_Centre_Bioinformatics/CNV_SGUL_Exomes/bedtools/IDs_ALL_OMIM.txt /apittman-aramis/sgul/shares/aramis/Genetics_Centre_Bioinformatics/CNV_SGUL_Exomes/SQL/combined_all_csv_full_ID_deldups.csv > /apittman-aramis/sgul/shares/aramis/Genetics_Centre_Bioinformatics/CNV_SGUL_Exomes/SQL/combined_all_csv_full_ID_rm_dups.csv

############### ROH #################

cat /apittman-aramis/sgul/shares/aramis/Genetics_Centre_Bioinformatics/CNV_SGUL_Exomes/bedtools/CNV_ROH_OMIM_annotated_point01_header.csv | cut -d ',' -f16 | tail -n +2 > /apittman-aramis/sgul/shares/aramis/Genetics_Centre_Bioinformatics/CNV_SGUL_Exomes/bedtools/IDs_ROH_OMIM.txt

cp /apittman-aramis/sgul/shares/aramis/Genetics_Centre_Bioinformatics/CNV_SGUL_Exomes/SQL/combined_roh_calls_2.csv /apittman-aramis/sgul/shares/aramis/Genetics_Centre_Bioinformatics/CNV_SGUL_Exomes/SQL/combined_roh_calls_2_deldups.csv

# delete from all_CSVs file any rows containing IDs found in the OMIM_all_CNVs
grep -vFf /apittman-aramis/sgul/shares/aramis/Genetics_Centre_Bioinformatics/CNV_SGUL_Exomes/bedtools/IDs_ROH_OMIM.txt /apittman-aramis/sgul/shares/aramis/Genetics_Centre_Bioinformatics/CNV_SGUL_Exomes/SQL/combined_roh_calls_2_deldups.csv > /apittman-aramis/sgul/shares/aramis/Genetics_Centre_Bioinformatics/CNV_SGUL_Exomes/SQL/combined_roh_calls_2_rm_dups.csv

# R script to rearrange columns, add new columns in and export (so the column order and names are identical) THEN combine the all CSVs and all OMIM and export
Rscript /apittman-aramis/sgul/shares/aramis/Genetics_Centre_Bioinformatics/CNV_SGUL_Exomes/scripts_misc/combine_OMIM_main_R.r &
pid=$!
wait ${pid}
# 


#while read ID; do
#
#	sed -i '/${ID}/d' /apittman-aramis/sgul/shares/aramis/Genetics_Centre_Bioinformatics/CNV_SGUL_Exomes/SQL/combined_all_csv_full_ID_deldups.csv
#
#done < /apittman-aramis/sgul/shares/aramis/Genetics_Centre_Bioinformatics/CNV_SGUL_Exomes/bedtools/IDs_ALL_OMIM.txt



