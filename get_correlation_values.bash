#!/bin/bash

for log in /apittman-aramis/sgul/shares/aramis/Genetics_Centre_Bioinformatics/CNV_SGUL_Exomes/logs/*_log.txt; do
	input=$log

	while IFS= read -r line
	do
		#	echo $line	
		if (grep "Now parsing /apittman-aramis/sgul/shares/aramis/Genetics_Centre_Bioinformatics/CNV_SGUL_Exomes/hg38_BAMs_to_convert/" <<< "$line"); then
			
			echo "$line" | sed 's,Now parsing /apittman-aramis/sgul/shares/aramis/Genetics_Centre_Bioinformatics/CNV_SGUL_Exomes/hg38_BAMs_to_convert/,,' | sed 's/_hg37aligned_sorted_rmduplicates.bam//' >> /apittman-aramis/sgul/shares/aramis/Genetics_Centre_Bioinformatics/CNV_SGUL_Exomes/QC/correlation_logs.txt
			
		fi
		
		if (grep "Correlation between reference and tests count is" <<< $line); then

			echo "$line" | sed 's/Correlation between reference and tests count is //' >> /apittman-aramis/sgul/shares/aramis/Genetics_Centre_Bioinformatics/CNV_SGUL_Exomes/QC/correlation_logs.txt

		fi
	
	done < "$input"

done

cat /apittman-aramis/sgul/shares/aramis/Genetics_Centre_Bioinformatics/CNV_SGUL_Exomes/QC/correlation_logs.txt | sed '/0.*/d' | sed 's/.*\///g' > /apittman-aramis/sgul/shares/aramis/Genetics_Centre_Bioinformatics/CNV_SGUL_Exomes/QC/correlation_logs_novals.txt

cat /apittman-aramis/sgul/shares/aramis/Genetics_Centre_Bioinformatics/CNV_SGUL_Exomes/QC/correlation_logs.txt | grep '\.' > /apittman-aramis/sgul/shares/aramis/Genetics_Centre_Bioinformatics/CNV_SGUL_Exomes/QC/correlation_logs_onlyvals.txt

paste -d"," /apittman-aramis/sgul/shares/aramis/Genetics_Centre_Bioinformatics/CNV_SGUL_Exomes/QC/correlation_logs_novals.txt /apittman-aramis/sgul/shares/aramis/Genetics_Centre_Bioinformatics/CNV_SGUL_Exomes/QC/correlation_logs_onlyvals.txt > /apittman-aramis/sgul/shares/aramis/Genetics_Centre_Bioinformatics/CNV_SGUL_Exomes/QC/correlation_logs_final.txt