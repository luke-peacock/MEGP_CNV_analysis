library(tidyverse)

samples_list <- list.files(path="S:/Genetics_Centre_Bioinformatics/CNV_SGUL_Exomes/CNV_calls")

calls_per_sample <- data.frame(matrix(0, ncol = 2, nrow = 0))

for (i in 1:length(samples_list)) {
  
  rows_inc_header <- nrow(read.csv(paste("S:/Genetics_Centre_Bioinformatics/CNV_SGUL_Exomes/CNV_calls/", samples_list[i], sep = "")))
  rows_no_header <- rows_inc_header - 1
  
  CNV_id <- samples_list[i]
  CNV_id <- str_extract(CNV_id, '.*(?=\\.csv)')
  
  row_to_append <- c(CNV_id, as.numeric(rows_no_header))
  
  calls_per_sample <- rbind(calls_per_sample, row_to_append)
  
}

colnames(calls_per_sample) <- c("sample_ID", "no_CNV_calls")

setwd("S:/Genetics_Centre_Bioinformatics/CNV_SGUL_Exomes/R/diagrams")

ggplot(calls_per_sample, aes(x=no_CNV_calls)) + geom_bar()

#geom_density(alpha=1) 

summary(calls_per_sample$no_CNV_calls)

p

#png("histo_reads_ratios_to_5.5.png", width = 1200, height = 800, res = 120)
#plot(p)
#dev.off()

library(dplyr)
summarise(calls_per_sample)
