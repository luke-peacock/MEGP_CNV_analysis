library(tidyverse)
library(data.table)

#import data frames
OMIM_all <- as.data.frame(read.csv("S://Genetics_Centre_Bioinformatics/CNV_SGUL_Exomes/bedtools/CNV_ALL_OMIM_annotated_point01_header.csv",header = TRUE, sep=",",stringsAsFactors=FALSE))
OMIM_roh <- as.data.frame(read.csv("S://Genetics_Centre_Bioinformatics/CNV_SGUL_Exomes/bedtools/CNV_ROH_OMIM_annotated_point01_header.csv",header = TRUE, sep=",",stringsAsFactors=FALSE))
full_all <- as.data.frame(read.csv("S://Genetics_Centre_Bioinformatics/CNV_SGUL_Exomes/SQL/combined_all_csv_full_ID.csv",header = TRUE, sep=",",stringsAsFactors=FALSE))
full_roh <- as.data.frame(read.csv("S://Genetics_Centre_Bioinformatics/CNV_SGUL_Exomes/SQL/combined_roh_calls_2.csv",header = TRUE, sep=",",stringsAsFactors=FALSE))

full_all2 <- full_all[c("batch_ID", "ID", "start.p", "end.p", "type", "nexons", "start", "end", "chromosome", "ID_CNV_id_type", "BF", "reads.expected", "reads.observed", "reads.ratio", "Conrad.hg19", "exons.hg19")]
full_all2 <- full_all2 %>% add_column(chrom_OMIM = NA,
                         OMIM_start = NA,
                         OMIM_end = NA,
                         OMIM_gene_name = NA,
                         OMIM_phenotype = NA,
                         OMIM_bp_overlap = NA)

OMIM_all2 <- OMIM_all[c("batch_ID", "sample_ID", "start.p", "end.p", "CNV_type", "nexons", "CNV_start", "CNV_end", "CNV_chromosome", "CNV_ID", "BF", "reads.expected", "reads.observed", "reads.ratio", "Conrad.hg19", "exons.hg19", "chrom_OMIM", "OMIM_start", "OMIM_end", "OMIM_gene_name", "OMIM_phenotype", "OMIM_bp_overlap")]

setnames(full_all2, old = c("batch_ID", "ID", "start.p", "end.p", "type", "nexons", "start", "end", "chromosome", "ID_CNV_id_type", "BF", "reads.expected", "reads.observed", "reads.ratio", "Conrad.hg19", "exons.hg19", "chrom_OMIM", "OMIM_start", "OMIM_end", "OMIM_gene_name", "OMIM_phenotype", "OMIM_bp_overlap"), new = c("batch_ID", "sample_ID", "start.p", "end.p", "CNV_type", "nexons", "CNV_start", "CNV_end", "CNV_chromosome", "CNV_ID", "BF", "reads.expected", "reads.observed", "reads.ratio", "Conrad.hg19", "exons.hg19", "chrom_OMIM", "OMIM_start", "OMIM_end", "OMIM_gene_name", "OMIM_phenotype", "OMIM_bp_overlap"))

for (row in 1:nrow(OMIM_all2)) {
  row <- OMIM_all2[row,]

  #delete df lines containing CNV ID
  full_all2[!grepl(row$CNV_ID, full_all2$CNV_ID),]
  
}
#combine dataframes
combined <- rbind(OMIM_all2, full_all2)




