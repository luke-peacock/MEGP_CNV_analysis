#!/bin/bash

dir_path="/apittman-aramis/sgul/shares/aramis/Genetics_Centre_Bioinformatics/CNV_SGUL_Exomes"
bedtools_software="/homes/athosnew/Genetics_Centre_Bioinformatics/resourses/bedtools2/bin/bedtools"

#SQL query to fetch all sample IDs
sqlite3 ${dir_path}/SQL/CNV_data.db "SELECT DISTINCT ID from fullcombined;" > ${dir_path}/SQL/ID_all_samples.txt

#start while loop

while read sample; do

  #check if bed intersect output exists using if statement
  bed_file="${dir_path}/bedtools/cnvs_roh/${sample}_cnvs_in_roh.bed"
  if [ -f "${bed_file}" ]; then
	#ROH intersect bed file already exists
    continue
  else
	#if does not exist then awk select rows in ROH file for each sample ID AND put into temp sample ROH file
	grep "${sample}" ${dir_path}/ROH_calls/bcftools_roh_output | awk -F'\t' -v OFS="\t" '{print $3, $4, $5, $2, $6, $7, $8}' | sed '/chrX/d' | sed '/chrY/d' > ${dir_path}/ROH_calls/temp_"${sample}"_roh.bed #this is sample roh path

	#make temp bed file for the sample
	grep "${sample}" ${dir_path}/SQL/combined_all_csv.csv | sed 's/,/	/g' | sed 's/"//g' | awk -F'\t' -v OFS="\t" '{print $9, $7, $8, $10, $1, $2, $3, $4, $5, $6, $11, $12, $13, $14, $15, $16}' | sed 's/^/chr/' > ${dir_path}/ROH_calls/temp_"${sample}"_calls.bed #this is sample calls path

	#run intersect AND output to dir AND remove temp ROH file for sample
	${bedtools_software} intersect -a ${dir_path}/ROH_calls/temp_"${sample}"_roh.bed -b ${dir_path}/ROH_calls/temp_"${sample}"_calls.bed -wo -F 1 > ${dir_path}/bedtools/cnvs_roh/${sample}_cnvs_in_roh.bed
	
	rm ${dir_path}/ROH_calls/temp_"${sample}"_roh.bed
	rm ${dir_path}/ROH_calls/temp_"${sample}"_calls.bed
  
  fi

done < ${dir_path}/SQL/ID_all_samples.txt

#set columns of combined csv
echo "ROH_chromosome,ROH_start,ROH_end,ROH_sample_ID,ROH_length,ROH_number_of_markers,ROH_quality_score,CNV_chromosome,CNV_start,CNV_end,CNV_ID,batch_ID,sample_ID,start.p,end.p,CNV_type,nexons,BF,reads.expected,reads.observed,reads.ratio,Conrad.hg19,exons.hg19,bp_overlap" > ${dir_path}/SQL/combined_roh_calls.csv

#combine into one csv or bed file
for ID in "/apittman-aramis/sgul/shares/aramis/Genetics_Centre_Bioinformatics/CNV_SGUL_Exomes/bedtools/cnvs_roh/*.bed"; do
	
	cat ${ID} | sed 's/,/ /g' | sed 's/	/,/g' >> ${dir_path}/SQL/combined_roh_calls.csv

done

awk 'BEGIN {FS=","; OFS=","}; {gsub(/"/,"",$11)}; {new_var=""$4"_"$11"_"$16""}; {print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, new_var, $12, $13, $14, $15, $16, $17, $18, $19, $20, $21, $22, $23, $24}' /apittman-aramis/sgul/shares/aramis/Genetics_Centre_Bioinformatics/CNV_SGUL_Exomes/SQL/combined_roh_calls.csv > /apittman-aramis/sgul/shares/aramis/Genetics_Centre_Bioinformatics/CNV_SGUL_Exomes/SQL/combined_roh_calls_2.csv

#import to SQL with overwrite of old 

python3 ${dir_path}/scripts_misc/import_roh_combined_csv.py &
pid=$!
wait ${pid}

exit