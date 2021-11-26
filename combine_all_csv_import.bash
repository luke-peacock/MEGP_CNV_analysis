#!/bin/bash

#echo 'batch_ID,ID,start.p,end.p,type,nexons,start,end,chromosome,CNV_id,BF,reads.expected,reads.observed,reads.ratio,Conrad.hg19,exons.hg19' > /apittman-aramis/sgul/shares/aramis/Genetics_Centre_Bioinformatics/CNV_SGUL_Exomes/SQL/combined_all_csv.csv
:'
for csvs in /apittman-aramis/sgul/shares/aramis/Genetics_Centre_Bioinformatics/CNV_SGUL_Exomes/combined_CNV_csv/*.csv; do
	input=$csvs
	
	while IFS= read -r line
	do
		if $(grep -q 'start.p' <<< $line); then
			continue
		else
			#echo "Line = $line"
			batch_ID=$csvs
			batch_ID=$(echo $batch_ID | sed 's|/apittman-aramis/sgul/shares/aramis/Genetics_Centre_Bioinformatics/CNV_SGUL_Exomes/combined_CNV_csv/||' <<< $batch_ID)
			batch_ID=$(echo $batch_ID | sed 's/.csv//' <<< $batch_ID)
			
			echo "$batch_ID,$line" >> /apittman-aramis/sgul/shares/aramis/Genetics_Centre_Bioinformatics/CNV_SGUL_Exomes/SQL/combined_all_csv.csv
		fi
	done < "$input"

done
'
awk 'BEGIN {FS=","; OFS=","}; {gsub(/"/,"",$10)}; {gsub(/"/,"",$5)}; {new_var="\""$2"_"$10"_"$5"\""}; {print $1, $2, $3, $4, $5, $6, $7, $8, $9, new_var, $11, $12, $13, $14, $15, substr($0, index($0,$16))}' /apittman-aramis/sgul/shares/aramis/Genetics_Centre_Bioinformatics/CNV_SGUL_Exomes/SQL/combined_all_csv.csv > /apittman-aramis/sgul/shares/aramis/Genetics_Centre_Bioinformatics/CNV_SGUL_Exomes/SQL/combined_all_csv_full_ID.csv

#awk 'BEGIN {FS=","; OFS=","}; {gsub(/"/,"",$10)}; {gsub(/"/,"",$5)}; {new_var="\""$2"_"$10"_"$5"\""}; {print $1, $2, $3, $4, $5, $6, $7, $8, $9, new_var, $11, $12, $13, $14, $15, substr($0, index($0,$16))}' /apittman-aramis/sgul/shares/aramis/Genetics_Centre_Bioinformatics/CNV_SGUL_Exomes/SQL/combined_all_csv.csv > /apittman-aramis/sgul/shares/aramis/Genetics_Centre_Bioinformatics/CNV_SGUL_Exomes/SQL/combined_all_csv_full_ID.csv
#  | cut -d, -f18 --complement
#| cut -d, -f17 --complement /apittman-aramis/sgul/shares/aramis/Genetics_Centre_Bioinformatics/CNV_SGUL_Exomes/SQL/combined_all_csv.csv | less

echo "csv has been combined"

#import table to SQL
python3 /apittman-aramis/sgul/shares/aramis/Genetics_Centre_Bioinformatics/CNV_SGUL_Exomes/scripts_misc/import_combined_csv.py

echo "csv has been imported, check SQL"

