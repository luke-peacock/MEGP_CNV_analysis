#!/bin/bash
######################						 UPDATE DATABASE SCRIPT 									######################
###################### Luke Peacock (lukepeacock98@gmail.com) and Dr Alan Pittman (apittman@sgul.ac.uk) ######################

echo "******************** PLEASE RUN THIS SCRIPT AS A BACKGROUND PROCESS USING COMMAND LINE AS BELOW ******************"
echo "******************** nohup bash /apittman-aramis/sgul/shares/aramis/Genetics_Centre_Bioinformatics/CNV_SGUL_Exomes/scripts_misc/update_db.bash ******************"

cat /apittman-aramis/sgul/shares/aramis/Genetics_Centre_Bioinformatics/CNV_SGUL_Exomes/scripts_misc/screen_update.txt

while true; do
    read -p "The SQL database must not be currently connected to before running this script. Ok to continue (y/n)?" yn
    case $yn in
        [Yy]* ) echo "Continuing with update..."; break;;
        [Nn]* ) exit;;
        * ) echo "Please answer yes or no.";;
	esac
done

echo "Prerequisites: (1) All samples to be in the database have been processed through the CNV_calling_pipeline.bash script; (2) Joint calling VCF path in this script is the latest version (Line 57)"
sleep 2
echo "This script will now perform a full update of CNV frequencies, OMIM annotation, ROH calling and import updated tables into the SQL database."
sleep 2
echo "Run this script as a background process as it may take several hours; The user who runs this script will receive an email on completion of the update"

###################### 						COMBINE CNV calls 											######################
rm /apittman-aramis/sgul/shares/aramis/Genetics_Centre_Bioinformatics/CNV_SGUL_Exomes/config*.csv
echo 'batch_ID,ID,start.p,end.p,type,nexons,start,end,chromosome,CNV_id,BF,reads.expected,reads.observed,reads.ratio,Conrad.hg19,exons.hg19' > /apittman-aramis/sgul/shares/aramis/Genetics_Centre_Bioinformatics/CNV_SGUL_Exomes/SQL/combined_all_csv.csv

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

#control to ensure no extra fields at end (has been problem in past)
cat /apittman-aramis/sgul/shares/aramis/Genetics_Centre_Bioinformatics/CNV_SGUL_Exomes/SQL/combined_all_csv.csv | sed 's/NA,NA,NA/NA,NA/g' > /apittman-aramis/sgul/shares/aramis/Genetics_Centre_Bioinformatics/CNV_SGUL_Exomes/SQL/combined_all_csv_temp.csv
sleep 2
mv /apittman-aramis/sgul/shares/aramis/Genetics_Centre_Bioinformatics/CNV_SGUL_Exomes/SQL/combined_all_csv_temp.csv /apittman-aramis/sgul/shares/aramis/Genetics_Centre_Bioinformatics/CNV_SGUL_Exomes/SQL/combined_all_csv.csv

#create CNV ID
awk 'BEGIN {FS=","; OFS=","}; {gsub(/"/,"",$10)}; {gsub(/"/,"",$5)}; {new_var="\""$2"_"$10"_"$5"\""}; {print $1, $2, $3, $4, $5, $6, $7, $8, $9, new_var, $11, $12, $13, $14, $15, substr($0, index($0,$16))}' /apittman-aramis/sgul/shares/aramis/Genetics_Centre_Bioinformatics/CNV_SGUL_Exomes/SQL/combined_all_csv.csv > /apittman-aramis/sgul/shares/aramis/Genetics_Centre_Bioinformatics/CNV_SGUL_Exomes/SQL/combined_all_csv_full_ID.csv

#control to ensure no extra fields at end (has been problem in past)
cat /apittman-aramis/sgul/shares/aramis/Genetics_Centre_Bioinformatics/CNV_SGUL_Exomes/SQL/combined_all_csv_full_ID.csv | sed 's/NA,NA,NA/NA,NA/g' > /apittman-aramis/sgul/shares/aramis/Genetics_Centre_Bioinformatics/CNV_SGUL_Exomes/SQL/combined_all_csv_full_ID_temp.csv
sleep 2
mv /apittman-aramis/sgul/shares/aramis/Genetics_Centre_Bioinformatics/CNV_SGUL_Exomes/SQL/combined_all_csv_full_ID_temp.csv /apittman-aramis/sgul/shares/aramis/Genetics_Centre_Bioinformatics/CNV_SGUL_Exomes/SQL/combined_all_csv_full_ID.csv

echo "csv has been combined"

#import table to SQL
python3 /apittman-aramis/sgul/shares/aramis/Genetics_Centre_Bioinformatics/CNV_SGUL_Exomes/scripts_misc/import_combined_csv.py &
pid=$!
wait ${pid}
echo "csv has been imported, check SQL"

###################### 						CALL ROHs from joint vcf calling 											######################
echo "Calling ROHs from joint-called vcf"
#OLD   ******* /homes/athosnew/Genetics_Centre_Bioinformatics/resourses/bcftools-1.9/bcftools roh -I -O r -o /apittman-aramis/sgul/shares/aramis/Genetics_Centre_Bioinformatics/CNV_SGUL_Exomes/ROH_calls/bcftools_roh_output /apittman-aramis/sgul/shares/aramis/Alan/SGUL_Exomes_R3_April_2021/Sample_characterisation/GRC_v3_16-3-2021.joint-called.vqsr_split_leftAligned.annotated.PASS.vcf.gz &
/homes/athosnew/Genetics_Centre_Bioinformatics/resourses/bcftools-1.9/bcftools roh -I -O r -o /apittman-aramis/sgul/shares/aramis/Genetics_Centre_Bioinformatics/CNV_SGUL_Exomes/ROH_calls/bcftools_roh_output /apittman-aramis/sgul/shares/aramis/Alan/SGUL_Exomes_R4_November_2021/GRC_Oct_2021.joint-called.vqsr_split_leftAligned.vcf.gz &
pid=$!
wait ${pid}

###################### 						Intersect ROHs vs all samples 											######################
echo "intersecting ROHs for all CNVs"
rm /apittman-aramis/sgul/shares/aramis/Genetics_Centre_Bioinformatics/CNV_SGUL_Exomes/bedtools/cnvs_roh/*.bed
bash /apittman-aramis/sgul/shares/aramis/Genetics_Centre_Bioinformatics/CNV_SGUL_Exomes/scripts_misc/allsamples_intersect_roh.bash &
pid=$!
wait ${pid}

###################### 						Intersect OMIM overlap 											######################
echo "intersecting OMIM for all CNVs"
bash /apittman-aramis/sgul/shares/aramis/Genetics_Centre_Bioinformatics/CNV_SGUL_Exomes/scripts_misc/OMIM_CNV_overlap.bash &
pid=$!
wait ${pid}

###################### 						Combine OMIM annotations 											######################
echo "Combining OMIM annotations"
bash /apittman-aramis/sgul/shares/aramis/Genetics_Centre_Bioinformatics/CNV_SGUL_Exomes/scripts_misc/combine_OMIM_main.bash &
pid=$!
wait ${pid}

###################### 						Collect frequencies for CNVs 											######################
rm /apittman-aramis/sgul/shares/aramis/Genetics_Centre_Bioinformatics/CNV_SGUL_Exomes/frequency_CNVs/output/*_frequencies.csv
echo "Calculating frequencies for each CNV"
bash /apittman-aramis/sgul/shares/aramis/Genetics_Centre_Bioinformatics/CNV_SGUL_Exomes/scripts_misc/batch_submit_frequencies.bash &
pid=$!
wait ${pid}

###################### 						UPDATE ALL TABLES SQLITE 											######################
echo "Updating tables in SQL database"
python3 /apittman-aramis/sgul/shares/aramis/Genetics_Centre_Bioinformatics/CNV_SGUL_Exomes/scripts_misc/update_db_python.py &
pid=$!
wait ${pid}
Rscript /apittman-aramis/sgul/shares/aramis/Genetics_Centre_Bioinformatics/CNV_SGUL_Exomes/scripts_misc/extract_PI_batch_inventory.R &
pid=$!
wait ${pid}

###################### 						SEND EMAIL 											######################
echo "Emailing user, update is complete."
echo "CNV Database Update finished" | mail -s "Finished updating CNV database :)" $USER@sgul.ac.uk




