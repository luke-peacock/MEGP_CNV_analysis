#!/bin/bash

count_todo=$(ls /apittman-aramis/sgul/shares/aramis/Genetics_Centre_Bioinformatics/CNV_SGUL_Exomes/combined_CNV_csv/ | wc -l)
count_output=$(ls /apittman-aramis/sgul/shares/aramis/Genetics_Centre_Bioinformatics/CNV_SGUL_Exomes/frequency_CNVs/output | wc -l)
#count_output= expr $count_output - 2

n=20

while [ ${count_output} -lt ${count_todo} ]; do

  for batches in /apittman-aramis/sgul/shares/aramis/Genetics_Centre_Bioinformatics/CNV_SGUL_Exomes/combined_CNV_csv/*.csv; do
	
	i=`basename "$batches" .csv`;  
	batch=`echo -e $i | sed 's/.csv//'`; 
	
	echo ${batch}
	
	count_pre=$(ls /apittman-aramis/sgul/shares/aramis/Genetics_Centre_Bioinformatics/CNV_SGUL_Exomes/combined_CNV_csv/ | wc -l)
	count_post=$(ls /apittman-aramis/sgul/shares/aramis/Genetics_Centre_Bioinformatics/CNV_SGUL_Exomes/frequency_CNVs/output | wc -l)
	#count_post=`expr $count_post - 2`
	
	
	diff_count=`expr $count_pre - $count_post`
	
	count_todo=$(ls /apittman-aramis/sgul/shares/aramis/Genetics_Centre_Bioinformatics/CNV_SGUL_Exomes/combined_CNV_csv/ | wc -l)
	count_output=$(ls /apittman-aramis/sgul/shares/aramis/Genetics_Centre_Bioinformatics/CNV_SGUL_Exomes/frequency_CNVs/output | wc -l)
	
  if [ ${count_output} -ge ${count_todo} ]; then
				
		echo "breaking from while loop, all batches processed!"
		break 3
		
	else
	
	if [ -f /apittman-aramis/sgul/shares/aramis/Genetics_Centre_Bioinformatics/CNV_SGUL_Exomes/frequency_CNVs/output/${batch}_frequencies.csv ]; then 
		continue
	else
		echo "inside first if loop"
		valu=$( ps -u m1606864 | grep R | wc -l )
		#echo "val =" ${valu}
		if [ $valu -ge $n ]; then
			echo "failed 2nd if"
			sleep 30
			continue
		else
			subbed=$(ps aux | grep m1606864 | grep ${batch} | wc -l)
			if [ $subbed -le 1 ]; then
				echo "passed 3rd if"
				cat /apittman-aramis/sgul/shares/aramis/Genetics_Centre_Bioinformatics/CNV_SGUL_Exomes/SQL/combined_all_csv_full_ID.csv | grep "reads.expected" > /apittman-aramis/sgul/shares/aramis/Genetics_Centre_Bioinformatics/CNV_SGUL_Exomes/frequency_CNVs/temps/${batch}_temp.csv
	
				cat /apittman-aramis/sgul/shares/aramis/Genetics_Centre_Bioinformatics/CNV_SGUL_Exomes/SQL/combined_all_csv_full_ID.csv | grep ${batch} >> /apittman-aramis/sgul/shares/aramis/Genetics_Centre_Bioinformatics/CNV_SGUL_Exomes/frequency_CNVs/temps/${batch}_temp.csv
	
				Rscript /apittman-aramis/sgul/shares/aramis/Genetics_Centre_Bioinformatics/CNV_SGUL_Exomes/scripts_misc/frequency_extract.R --args ${batch} &
			
				echo "job submitted" $batch
			else
				echo "sleeping post 3rd if"
				sleep 2
			
				if [ ${count_output} -ge ${count_todo} ]; then
				
					echo "breaking from while loop, all batches processed!"
					break 7
				
				fi
			
			fi
			echo "line 52"
			sleep 5
		fi
	fi
  fi
  done
	sleep 10
done

echo "PASSED WHILE LOOP"

if [ ${count_output} -eq ${count_todo} ]; then
	echo "Finished frequency calculations; now merging files"
	
	echo "CNV_ID, lowest_frequency, number_CNVs_at_lowest_exon, highest_frequency, number_CNVs_at_highest_exon, mean_frequency, mean_number_of_exons" > /apittman-aramis/sgul/shares/aramis/Genetics_Centre_Bioinformatics/CNV_SGUL_Exomes/frequency_CNVs/all_freq_combined.csv
	
	for freqs in /apittman-aramis/sgul/shares/aramis/Genetics_Centre_Bioinformatics/CNV_SGUL_Exomes/frequency_CNVs/output/*_frequencies.csv; do
		
		cat $freqs | sed 1d >> /apittman-aramis/sgul/shares/aramis/Genetics_Centre_Bioinformatics/CNV_SGUL_Exomes/frequency_CNVs/all_freq_combined.csv
	
	done
fi

cat /apittman-aramis/sgul/shares/aramis/Genetics_Centre_Bioinformatics/CNV_SGUL_Exomes/frequency_CNVs/all_freq_combined.csv | sed 's/,0,/,NA,/g' > /apittman-aramis/sgul/shares/aramis/Genetics_Centre_Bioinformatics/CNV_SGUL_Exomes/frequency_CNVs/all_freq_combined_NAs.csv

echo "Script finished and frequencies combined and now importing."

python3 /apittman-aramis/sgul/shares/aramis/Genetics_Centre_Bioinformatics/CNV_SGUL_Exomes/scripts_misc/import_frequencies.py &
pid=$!
wait ${pid}

echo "Imported. Script finished. Now exiting..."


exit

