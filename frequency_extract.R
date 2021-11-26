library(tidyverse)
library(DBI)
library(RSQLite)

batchID <- commandArgs(trailingOnly = TRUE)
batchID <- batchID[2]

path_csv <- paste("/apittman-aramis/sgul/shares/aramis/Genetics_Centre_Bioinformatics/CNV_SGUL_Exomes/frequency_CNVs/temps/", batchID, "_temp.csv", sep = "")
all_CNVs <- read.csv(path_csv)
del_freq <- as.data.frame(read.table("/apittman-aramis/sgul/shares/aramis/Genetics_Centre_Bioinformatics/CNV_SGUL_Exomes/bedtools/deletions_intersect_batches_overlap90.bed",header = FALSE, sep="\t",stringsAsFactors=FALSE, quote=""))
dup_freq <- as.data.frame(read.table("/apittman-aramis/sgul/shares/aramis/Genetics_Centre_Bioinformatics/CNV_SGUL_Exomes/bedtools/duplications_intersect_batches_overlap90.bed",header = FALSE, sep="\t",stringsAsFactors=FALSE, quote=""))

all_CNVs <- as.data.frame(all_CNVs)

CNV_db <- dbConnect(RSQLite::SQLite(), '/apittman-aramis/sgul/shares/aramis/Genetics_Centre_Bioinformatics/CNV_SGUL_Exomes/SQL/CNV_data.db')
all_IDs <- dbGetQuery(CNV_db, 'SELECT DISTINCT ID from fullcombined')
#all_IDs <- read.delim("/apittman-aramis/sgul/shares/aramis/Genetics_Centre_Bioinformatics/CNV_SGUL_Exomes/SQL/ID_all_samples.txt")

dbDisconnect(CNV_db)

df <- data.frame(matrix(ncol = 7, nrow = 0))
colnames(df) <- c("CNV_ID", "lowest_frequency", "number_CNVs_at_lowest_exon", "highest_frequency", "number_CNVs_at_highest_exon", "mean_frequency", "mean_number_of_exons")

file_x = paste("/apittman-aramis/sgul/shares/aramis/Genetics_Centre_Bioinformatics/CNV_SGUL_Exomes/frequency_CNVs/output/", batchID, "_frequencies.csv", sep = "")

for (CNV_row in 1:nrow(all_CNVs)) {
  
  row <- all_CNVs[CNV_row,]
  

  exons_string <- row$exons.hg19
  exons_string <- as.character(exons_string)
  exons_vector <- strsplit(exons_string, ",")[[1]]
  
  if (row[["type"]]=="duplication") {
      
      freq_vector <- vector()
      num_vector <- vector()
	  
      for (exon in exons_vector) {
        
		grep_exon <- noquote(exon)
		
        exon_row_freq <- dup_freq[grep(grep_exon, dup_freq$V4), ]
		
		number <- exon_row_freq$V5
		num_vector <- append(num_vector, number)
		
        freq <- (exon_row_freq$V5) / (nrow(all_IDs))
        
        freq_vector <- append(freq_vector, freq)
        
      }
      
      output_freq_min <- min(freq_vector)
	  output_num_min <- min(num_vector)
      output_freq_max <- max(freq_vector)
	  output_num_max <- max(num_vector)	  
      output_freq_mean <- mean(freq_vector)
	  output_num_mean <- mean(num_vector)	

		new_ID <- as.character(row$ID_CNV_id)
      #new_row <- c(new_ID, output_freq)
	  new_row <- data.frame(new_ID, output_freq_min, output_num_min, output_freq_max, output_num_max, output_freq_mean, output_num_mean)
	  names(new_row)<-c("CNV_ID", "lowest_frequency", "number_CNVs_at_lowest_exon", "highest_frequency", "number_CNVs_at_highest_exon", "mean_frequency", "mean_number_of_exons")
	  df <- rbind(df, new_row)
      #print(df)
	  
	  if (file.exists(file_x)) {
			stop()
		}
    }

    
    if (row[["type"]]=="deletion") {
    
      freq_vector <- vector()
	  num_vector <- vector()
      
      for (exon in exons_vector) {
		
		grep_exon <- as.character(noquote(exon))
		
        exon_row_freq <- del_freq[grep(grep_exon, del_freq$V4), ]
		
		number <- exon_row_freq$V5
		
        freq <- (exon_row_freq$V5) / (nrow(all_IDs))

        freq_vector <- append(freq_vector, freq)
		num_vector <- append(num_vector, number)
        
    
      }
		
      output_freq_min <- min(freq_vector)
	  output_num_min <- min(num_vector)
      output_freq_max <- max(freq_vector)
	  output_num_max <- max(num_vector)	  
      output_freq_mean <- mean(freq_vector)
	  output_num_mean <- mean(num_vector)

	  new_ID <- as.character(row$ID_CNV_id)
      #new_row <- c(new_ID, output_freq)
	  new_row <- data.frame(new_ID, output_freq_min, output_num_min, output_freq_max, output_num_max, output_freq_mean, output_num_mean)
	  names(new_row)<-c("CNV_ID", "lowest_frequency", "number_CNVs_at_lowest_exon", "highest_frequency", "number_CNVs_at_highest_exon", "mean_frequency", "mean_number_of_exons")
	  
      df <- rbind(df, new_row)
	  #print(df)
		
		if (file.exists(file_x)) {
			stop()
		}
    }

}
colnames(df) <- c("CNV_ID", "lowest_frequency", "number_CNVs_at_lowest_exon", "highest_frequency", "number_CNVs_at_highest_exon", "mean_frequency", "mean_number_of_exons")
write.csv(df, file = paste("/apittman-aramis/sgul/shares/aramis/Genetics_Centre_Bioinformatics/CNV_SGUL_Exomes/frequency_CNVs/output/", batchID, "_frequencies.csv", sep = ""), row.names = FALSE)