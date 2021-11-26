library(tidyverse)
library(data.table)

#import data frames
OMIM_all <- as.data.frame(read.csv("/apittman-aramis/sgul/shares/aramis/Genetics_Centre_Bioinformatics/CNV_SGUL_Exomes/bedtools/CNV_ALL_OMIM_annotated_point01_header.csv",header = TRUE, sep=",",stringsAsFactors=FALSE))
OMIM_roh <- as.data.frame(read.csv("/apittman-aramis/sgul/shares/aramis/Genetics_Centre_Bioinformatics/CNV_SGUL_Exomes/bedtools/CNV_ROH_OMIM_annotated_point01_header.csv",header = TRUE, sep=",",stringsAsFactors=FALSE))
full_all <- as.data.frame(read.csv("/apittman-aramis/sgul/shares/aramis/Genetics_Centre_Bioinformatics/CNV_SGUL_Exomes/SQL/combined_all_csv_full_ID_rm_dups.csv",header = TRUE, sep=",",stringsAsFactors=FALSE))
full_roh <- as.data.frame(read.csv("/apittman-aramis/sgul/shares/aramis/Genetics_Centre_Bioinformatics/CNV_SGUL_Exomes/SQL/combined_roh_calls_2_rm_dups.csv",header = TRUE, sep=",",stringsAsFactors=FALSE))

# fix and export ALL_CNVS ###########
full_all2 <- full_all[c("batch_ID", "ID", "start.p", "end.p", "type", "nexons", "start", "end", "chromosome", "ID_CNV_id_type", "BF", "reads.expected", "reads.observed", "reads.ratio", "Conrad.hg19", "exons.hg19")]
full_all2 <- full_all2 %>% add_column(chrom_OMIM = NA,
                         OMIM_start = NA,
                         OMIM_end = NA,
                         OMIM_gene_name = NA,
                         OMIM_phenotype = NA,
                         OMIM_bp_overlap = NA)

OMIM_all2 <- OMIM_all[c("batch_ID", "sample_ID", "start.p", "end.p", "CNV_type", "nexons", "CNV_start", "CNV_end", "CNV_chromosome", "CNV_ID", "BF", "reads.expected", "reads.observed", "reads.ratio", "Conrad.hg19", "exons.hg19", "chrom_OMIM", "OMIM_start", "OMIM_end", "OMIM_gene_name", "OMIM_phenotype", "OMIM_bp_overlap")]

# rename ALL_CNVS to have same column names as OMIM_all
setnames(full_all2, old = c("batch_ID", "ID", "start.p", "end.p", "type", "nexons", "start", "end", "chromosome", "ID_CNV_id_type", "BF", "reads.expected", "reads.observed", "reads.ratio", "Conrad.hg19", "exons.hg19", "chrom_OMIM", "OMIM_start", "OMIM_end", "OMIM_gene_name", "OMIM_phenotype", "OMIM_bp_overlap"), new = c("batch_ID", "sample_ID", "start.p", "end.p", "CNV_type", "nexons", "CNV_start", "CNV_end", "CNV_chromosome", "CNV_ID", "BF", "reads.expected", "reads.observed", "reads.ratio", "Conrad.hg19", "exons.hg19", "chrom_OMIM", "OMIM_start", "OMIM_end", "OMIM_gene_name", "OMIM_phenotype", "OMIM_bp_overlap"))

write.csv(full_all2, file = paste("/apittman-aramis/sgul/shares/aramis/Genetics_Centre_Bioinformatics/CNV_SGUL_Exomes/SQL/combined_all_csv_full_ID_reordered.csv", sep = ""), row.names = FALSE)
write.csv(OMIM_all2, file = paste("/apittman-aramis/sgul/shares/aramis/Genetics_Centre_Bioinformatics/CNV_SGUL_Exomes/bedtools/CNV_ALL_OMIM_annotated_point01_header_reordered.csv", sep = ""), row.names = FALSE)

# fix and export ROH_CNVs ############
full_roh2 <- full_roh[c("batch_ID", "sample_ID", "start.p", "end.p", "CNV_type", "nexons", "CNV_start", "CNV_end", "CNV_chromosome", "ROH_sample_ID_CNV_ID_CNV_type", "ROH_chromosome", "ROH_start", "ROH_end", "ROH_sample_ID", "ROH_length", "ROH_number_of_markers", "ROH_quality_score", "BF", "reads.expected", "reads.observed", "reads.ratio", "Conrad.hg19", "exons.hg19", "bp_overlap")]
full_roh2 <- full_roh2 %>% add_column(chrom_OMIM = NA,
                         OMIM_start = NA,
                         OMIM_end = NA,
                         OMIM_gene_name = NA,
                         OMIM_phenotype = NA,
                         OMIM_bp_overlap = NA)

OMIM_roh2 <- OMIM_roh[c("batch_ID", "sample_ID", "start.p", "end.p", "CNV_type", "nexons", "CNV_start", "CNV_end", "CNV_chromosome", "CNV_ID", "ROH_chromosome", "ROH_start", "ROH_end", "ROH_sample_ID", "ROH_length", "ROH_number_of_markers", "ROH_quality_score", "BF", "reads.expected", "reads.observed", "reads.ratio", "Conrad.hg19", "exons.hg19", "roh_bp_overlap", "chrom_OMIM", "OMIM_start", "OMIM_end", "OMIM_gene_name", "OMIM_phenotype", "OMIM_bp_overlap")]

# rename ALL_CNVS to have same column names as OMIM_all
setnames(full_roh2, old = c("batch_ID", "sample_ID", "start.p", "end.p", "CNV_type", "nexons", "CNV_start", "CNV_end", "CNV_chromosome", "ROH_sample_ID_CNV_ID_CNV_type", "ROH_chromosome", "ROH_start", "ROH_end", "ROH_sample_ID", "ROH_length", "ROH_number_of_markers", "ROH_quality_score", "BF", "reads.expected", "reads.observed", "reads.ratio", "Conrad.hg19", "exons.hg19", "bp_overlap", "chrom_OMIM", "OMIM_start", "OMIM_end", "OMIM_gene_name", "OMIM_phenotype", "OMIM_bp_overlap"), new = c("batch_ID", "sample_ID", "start.p", "end.p", "CNV_type", "nexons", "CNV_start", "CNV_end", "CNV_chromosome", "CNV_ID", "ROH_chromosome", "ROH_start", "ROH_end", "ROH_sample_ID", "ROH_length", "ROH_number_of_markers", "ROH_quality_score", "BF", "reads.expected", "reads.observed", "reads.ratio", "Conrad.hg19", "exons.hg19", "roh_bp_overlap", "chrom_OMIM", "OMIM_start", "OMIM_end", "OMIM_gene_name", "OMIM_phenotype", "OMIM_bp_overlap"))

write.csv(full_roh2, file = paste("/apittman-aramis/sgul/shares/aramis/Genetics_Centre_Bioinformatics/CNV_SGUL_Exomes/SQL/combined_roh_calls_2_reordered.csv", sep = ""), row.names = FALSE)
write.csv(OMIM_roh2, file = paste("/apittman-aramis/sgul/shares/aramis/Genetics_Centre_Bioinformatics/CNV_SGUL_Exomes/bedtools/CNV_ROH_OMIM_annotated_point01_header_reordered.csv", sep = ""), row.names = FALSE)


#for (row in 1:nrow(OMIM_all2)) {
#  row <- OMIM_all2[row,]
#
#  #delete df lines containing CNV ID
#  full_all2[!grepl(row$CNV_ID, full_all2$CNV_ID),]
#  
#}
#combine dataframes
combined_frame <- rbind(OMIM_all2, full_all2)
write.csv(combined_frame, file = paste("/apittman-aramis/sgul/shares/aramis/Genetics_Centre_Bioinformatics/CNV_SGUL_Exomes/SQL/combined_OMIM_main.csv", sep = ""), row.names = FALSE)

combined_frame_roh <- rbind(OMIM_roh2, full_roh2)
write.csv(combined_frame_roh, file = paste("/apittman-aramis/sgul/shares/aramis/Genetics_Centre_Bioinformatics/CNV_SGUL_Exomes/SQL/combined_OMIM_roh.csv", sep = ""), row.names = FALSE)

#write.csv(combined_frame, file = paste("/apittman-aramis/sgul/shares/aramis/Genetics_Centre_Bioinformatics/CNV_SGUL_Exomes/SQL/combined_frame.csv", sep = ""), row.names = TRUE)