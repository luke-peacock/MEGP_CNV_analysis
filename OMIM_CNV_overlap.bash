#!/bin/bash
# cut -d "," -f 5- omim_phenotypes_modified.csv | sed 's/ /_/g' | sed 's/"//g' | sed 's/,/;/g' | less #to generate OMIM file final column

#sed -i '0,/string/{//d}' file 		#this deletes the first occurrence of a line with particular string

dir_path="/apittman-aramis/sgul/shares/aramis/Genetics_Centre_Bioinformatics/CNV_SGUL_Exomes"
bedtools_software="/homes/athosnew/Genetics_Centre_Bioinformatics/resourses/bedtools2/bin/bedtools"

#make bed of combined_roh_calls

#cat ${dir_path}/SQL/combined_roh_calls_2.csv | awk -F',' '{print $8, $9, $10, $1, $2, $3, $4, $5, $6, $7, $11, $12, $13, $14, $15, $16, $17, $18, $19, $20, $21, $22, $23, $24}' OFS='\t' | sed 's/ /NA/g' > ${dir_path}/SQL/combined_roh_calls_with_header.bed

#cat ${dir_path}/SQL/combined_roh_calls_with_header.bed | awk NR\>1 > ${dir_path}/SQL/combined_roh_calls_no_header.bed

#intersect all ROH CNVs with OMIM


#${bedtools_software} intersect -a ${dir_path}/OMIM/omim_phenotypes_modified_bed_withNAs.bed -b ${dir_path}/SQL/combined_roh_calls_no_header.bed -wo -F 0.0001 > ${dir_path}/bedtools/CNV_ROH_OMIM_annotated_point01.bed

#cat ${dir_path}/bedtools/CNV_ROH_OMIM_annotated_point01.bed | sed 's/\t/,/g' > ${dir_path}/bedtools/CNV_ROH_OMIM_annotated_point01.csv

#100%, 50%, 1%, 0.001%
#0.01 = 1%
#0.0001 = 0.01%



#make bed of combined_ALL_calls

cat ${dir_path}/SQL/combined_all_csv_full_ID.csv | awk -vFPAT='([^,]*)|("[^"]+")' -vOFS='\t' '{$9=gsub("\"","",$9)}; {$9="chr"$9}; {print $9, $7, $8, $1, $2, $3, $4, $5, $6, $10, $11, $12, $13, $14, $15, $16}' | sed 's/ /NA/g' > ${dir_path}/SQL/combined_all_csv_full_ID_with_header.bed

cat ${dir_path}/SQL/combined_all_csv_full_ID_with_header.bed | awk NR\>1 > ${dir_path}/SQL/combined_all_csv_full_ID_no_header.bed

#intersect all CNVs with OMIM


${bedtools_software} intersect -a ${dir_path}/OMIM/omim_phenotypes_modified_bed_withNAs.bed -b ${dir_path}/SQL/combined_all_csv_full_ID_no_header.bed -wo -F 0.0001 > ${dir_path}/bedtools/CNV_ALL_OMIM_annotated_point01.bed

cat ${dir_path}/bedtools/CNV_ALL_OMIM_annotated_point01.bed | sed 's/\t/,/g' > ${dir_path}/bedtools/CNV_ALL_OMIM_annotated_point01.csv


#/homes/athosnew/Genetics_Centre_Bioinformatics/resourses/bedtools2/bin/bedtools intersect -a /apittman-aramis/sgul/shares/aramis/Genetics_Centre_Bioinformatics/CNV_SGUL_Exomes/OMIM/omim_phenotypes_modified_bed_withNAs.bed -b /apittman-aramis/sgul/shares/aramis/Genetics_Centre_Bioinformatics/CNV_SGUL_Exomes/SQL/combined_all_csv_full_ID_no_header.bed -wo -F 0.0001 > /apittman-aramis/sgul/shares/aramis/Genetics_Centre_Bioinformatics/CNV_SGUL_Exomes/bedtools/CNV_ALL_OMIM_annotated_point01.bed


#cat /apittman-aramis/sgul/shares/aramis/Genetics_Centre_Bioinformatics/CNV_SGUL_Exomes/SQL/combined_all_csv_full_ID.csv | awk -vFPAT='([^,]*)|("[^"]+")' -vOFS='\t' '{print $8, $9, $10, $1, $2, $3, $4, $5, $6, $7, $11, $12, $13, $14, $15, $16}' | sed 's/ /NA/g' > /apittman-aramis/sgul/shares/aramis/Genetics_Centre_Bioinformatics/CNV_SGUL_Exomes/SQL/combined_all_csv_full_ID_with_header.bed