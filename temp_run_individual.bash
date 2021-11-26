#!/bin/bash

# RUN COMMAND LINE AS (in square brackets need to be changed - delete brackets): nohup bash CNV_calling_pipeline.bash -b [batchID] >>& /apittman-aramis/sgul/shares/aramis/Genetics_Centre_Bioinformatics/CNV_SGUL_Exomes/logs/[batch_ID]_log.txt &
# find . -type f -size -100k | grep .bam | sed '/.bai/d' | less
#list batch samples below
list_of_samples="
J1913
"

#define variables:
hg37_ref="/homedirs-yjamshid/athosnew/Genetics_Centre_Bioinformatics/resourses/Genome_reference_files/human_g1k_v37.fasta"   #hg37 reference genome path

#define paths
working_dir="/apittman-aramis/sgul/shares/aramis/Genetics_Centre_Bioinformatics/CNV_SGUL_Exomes"

#path to directory containing original hg38 samples in directories (named by sample IDs)
#ARAMIS
#hg38_samples_dirpath="/apittman-aramis/sgul/shares/aramis/Genetics_Centre_Bioinformatics/exomes_hg38/Aligned/*"
#PORTHOS
#hg38_samples_dirpath="/homedirs-porthos/sgul/shares/incc/porthos/Genetics_Centre_Bioinformatics/HPC_exomes_hg38/Aligned/*"
hg38_samples_dirpath="/homedirs-porthos/sgul/shares/incc/porthos/Genetics_Centre_Bioinformatics/HPC_exomes_hg37_hg38_liftover/Aligned/*"
#ATHOS
#hg38_samples_dirpath="/homes/athosnew/Genetics_Centre_Bioinformatics/Exomes/Aligned/*"

resourses_path="/homedirs-yjamshid/athosnew/Genetics_Centre_Bioinformatics/resourses"

#program paths
BWAsoftware="/homes/athosnew/Genetics_Centre_Bioinformatics/resourses/bwa/bwa"
BWAindex="/homes/athosnew/Genetics_Centre_Bioinformatics/resourses/Genome_reference_files/human_g1k_v37.fasta"
samtoolsSoftware="/homes/athosnew/Genetics_Centre_Bioinformatics/resourses/samtools-1.8/samtools"
java="/homes/athosnew/Genetics_Centre_Bioinformatics/resourses/java/jre1.8.0_171/bin/java"
picard="/homes/athosnew/Genetics_Centre_Bioinformatics/resourses/picard-2.815/picard.jar"
gatk="/homes/athosnew/Genetics_Centre_Bioinformatics/resourses/gatk-4.0.4.0/gatk-package-4.0.4.0-local.jar"

# call batch_ID from command line
while getopts b: batch		# in command line: bash CNV_calling_pipeline.bash -b batchID
do
  case "${batch}"  in
    b) batch_no=${OPTARG};; 	#batch_no is variable for the batchID from command line

  esac
done

echo "Batch_no = ${batch_no}"

###### BLOCK OUT BELOW IF SAMPLES ARE IN HG37 #####
for sample_ID in $list_of_samples; do

	echo "Sample_ID ${sample_ID} currently going through pipeline..."

	#define sample_ID paths
		ID_to_copy_path=${hg38_samples_dirpath}/${sample_ID}/${sample_ID}_sorted_unique_recalibrated.bam #TO BE COPIED sample_ID path
		ID_copied_path=${working_dir}/hg38_BAMs_to_convert/${sample_ID}/${sample_ID}_sorted_unique_recalibrated.bam  #COPIED sample_ID path
		sample_dir=${working_dir}/hg38_BAMs_to_convert/${sample_ID}

	#echo of paths
		echo "Path of sample to be copied = ${ID_to_copy_path}";
		echo "Path of copied sample = ${ID_copied_path}";

	#copy bam file across
		echo "Making" ${sample_ID} "directory"
		mkdir ${working_dir}/hg38_BAMs_to_convert/${sample_ID}
		echo "Copying ${sample_ID} to working directory"
		cp ${ID_to_copy_path} ${ID_copied_path}

	#convert bam to fastq files (then delete bam file)
		echo "Converting ${sample_ID} BAM to FASTQ"
		#BAM_toFASTQ_command
		${java} -jar ${picard} SamToFastq INCLUDE_NON_PF_READS=true INCLUDE_NON_PRIMARY_ALIGNMENTS=false COMPRESSION_LEVEL=0 VALIDATION_STRINGENCY=SILENT QUIET=true \
		I=${sample_dir}/${sample_ID}_sorted_unique_recalibrated.bam FASTQ=${sample_dir}/${sample_ID}_1.fastq.gz SECOND_END_FASTQ=${sample_dir}/${sample_ID}_2.fastq.gz
		
		#rm ${sample_dir}/${sample_ID}_sorted_unique_recalibrated.bam

	#align to hg37 ref (output = SAM; delete fastq)
		echo "Aligning ${sample_ID} FASTQ to hg37 reference"
		#BWAcommand
		samHeader="@RG\\tID:${sample_ID}\\tSM:${sample_ID}\\tLB:${sample_ID}\\tPL:ILLUMINA"
		
		${BWAsoftware} mem ${BWAindex} ${sample_dir}/${sample_ID}_1.fastq.gz ${sample_dir}/${sample_ID}_2.fastq.gz -t 8 -R ${samHeader} -o ${sample_dir}/${sample_ID}_hg37aligned.sam
		
		#rm ${sample_dir}/${sample_ID}_1.fastq.gz
		#rm ${sample_dir}/${sample_ID}_2.fastq.gz

	#convert to bam - SAMtools? (output = BAM; delete SAM)
		echo "Converting ${sample_ID} SAM to BAM"
		#samtoolsViewCommand
		${samtoolsSoftware} view -Sb \
		${sample_dir}/${sample_ID}_hg37aligned.sam -o ${sample_dir}/${sample_ID}_hg37aligned.bam
		#rm ${sample_dir}/${sample_ID}_hg37aligned.sam
		
	#sort bam - samtools
		echo "Sorting ${sample_ID} BAM file"
		#samtoolsSortCommand
		${samtoolsSoftware} sort -m \
		5000000000 ${sample_dir}/${sample_ID}_hg37aligned.bam -o ${sample_dir}/${sample_ID}_hg37aligned_sorted.bam
		#rm ${sample_dir}/${sample_ID}_hg37aligned.bam
	
	#index bam
		echo "Indexing ${sample_ID} BAM file"
		#samtoolsINDEXCommand
		${samtoolsSoftware} index ${sample_dir}/${sample_ID}_hg37aligned_sorted.bam
		
	#remove PCR duplicates
	echo "Removing ${sample_ID} PCR duplicates"
	#picardMARKDUPLICATESCommand
	${java} -jar ${picard} MarkDuplicates \
		I=${sample_dir}/${sample_ID}_hg37aligned_sorted.bam \
		O=${sample_dir}/${sample_ID}_hg37aligned_sorted_rmduplicates.bam \
		M=${sample_dir}/${sample_ID}_marked_dup_metrics.txt
	
	#rm ${sample_dir}/${sample_ID}_hg37aligned_sorted.bam
	#rm ${sample_dir}/${sample_ID}_hg37aligned_sorted.bam.bai
	#rm ${sample_dir}/${sample_ID}_marked_dup_metrics.txt

	#index bam
	echo "Indexing ${sample_ID} BAM file"
	#samtoolsINDEXCommand
	${samtoolsSoftware} index ${sample_dir}/${sample_ID}_hg37aligned_sorted_rmduplicates.bam
	
	echo "Converted ${sample_ID} in ${batch_no} to hg37. Repeating for next sample"
done
###################################################

exit
