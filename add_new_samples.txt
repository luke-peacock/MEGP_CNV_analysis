######################################### HOW TO ADD NEW SAMPLES TO THE DATABASE #########################################
- For help: email lukepeacock98@gmail.com or apittman@sgul.ac.uk

# processing new batches of samples
1. Samples should be in bam format (either hg37 or hg38) - preferably in hg37
2. Use the SGUL_GeneticsCentre_Exome_Sample_Inventory.xlsx to separate the samples into batches of 11. Name these batches however you want (only use text, _ or - in the name)
   - Ideally the batch should not contain >11 samples since the power of ExomeDepth diminishes with >10 samples in the reference set.
   - ExomeDepth compares the sample being tested against all the other samples in the batch (reference samples) - this repeats until all samples are analysed
   - When the samples cannot be separated into batches of 11, try to keep them as close to 11 as possible. This doesn't need to be exact, but the power of ExomeDepth analysis will decrease as the batch has >11 samples.
3. Open (on aramis) S:/Genetics_Centre_Bioinformatics/CNV_SGUL_Exomes/CNV_calling_pipeline.bash using Notepad++
4. Delete the contents between the double quotation marks on Line 6
   (list_of_samples="
		    *DELETE TEXT HERE*
		    ")
5. Replace this text between the quotation marks on Line 6 with the list of sample IDs from SGUL_GeneticsCentre_Exome_Sample_Inventory.xlsx.
   - Each sample ID should be on a new line, with no characters before of after the sample ID (if you copy and paste directly from the SGUL_GeneticsCentre_Exome_Sample_Inventory.xlsx, this will be fine)
6. From Line 26-33 of the CNV_calling_pipeline.bash script, remove the hash before the location of your bam files. MAKE SURE the other paths are hashed out to avoid errors finding the bam files.
7. Code can be ignored/commented out by putting : ' (colon, space, apostrophe) before the text and ' (apostrophe) after the text to be ignored.
   - If your samples are in hg37 format, block out the text between (approximately) Line 56 and Line 152 - starting "for sample_ID in $list_of_samples; do" and ending at the next "done". Make sure lines 135-154 are not blocked out!
   - If your samples are in hg38 format, block out the text between (approximately) Line 135 and Line 154 - starting "for sample_ID in $list_of_samples; do" and ending at the next "done". Make sure lines 56-152 are not blocked out!
8. If you want an email each time the script finishes, type your email at the end of Line 177 where it says "example@sgul.ac.uk" (ensure the line is not hashed out)
8. That's it. You're ready to run the script - if all your samples are in the same format and in the same "/Aligned" directory, you just need to change the sample IDs for each batch.
9. The script should be run like this: (replace the [batchID] with the name you have given the batch - remove the square brackets!)
	nohup bash CNV_calling_pipeline.bash -b [batchID] >>& /apittman-aramis/sgul/shares/aramis/Genetics_Centre_Bioinformatics/CNV_SGUL_Exomes/logs/[batch_ID]_log.txt &

##### TROUBLESHOOTING ####
1. The most common problem currently with this script is a "command not found" error due to word splitting of the script code.
   - If the R script finished successfully, rerun the script with Line 169 hashed out
2. Most other errors are due to issues with the bam files. Make sure the bam file exists in "/apittman-aramis/sgul/shares/aramis/Genetics_Centre_Bioinformatics/CNV_SGUL_Exomes/hg38_BAMs_to_convert/{sample_ID}/{sample_ID}.bam"
   - If it doesn't exist, you'll need to remove this sample ID from the list of samples on Line 6 and re-run the script
   - If it exists but is very small in size (e.g. <1kb bam file), remove it from the list of samples on Line 6 and re-run the script
   - For other errors such as "no bam index found", you'll need to experiment by taking one of the samples out of the list on Line 6 - one at a time - and rerun the script.
	For help, email lukepeacock98@gmail.com or apittman@sgul.ac.uk
