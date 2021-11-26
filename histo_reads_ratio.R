
library(tidyverse)
setwd("S:/Genetics_Centre_Bioinformatics/CNV_SGUL_Exomes/R/diagrams")
df <- read.csv("S:/Genetics_Centre_Bioinformatics/CNV_SGUL_Exomes/SQL/combined_all_csv_full_ID.csv")

df_1.5 <- subset(df, reads.ratio>=0 & reads.ratio<5.5)

p <- ggplot(df_1.5, aes(x=reads.ratio, fill=type), xlim(0, 1.5)) +
  geom_histogram(binwidth = 0.01, alpha=0.8) +
  scale_fill_manual(values=c("#FF0000", "#7CFC00")) +
  labs(title="Histogram plotting reads ratio of all CNVs",x="Reads Ratio", y = "Count (no. of CNVs)") +
  theme(plot.title = element_text(hjust = 0.5))
  
  #geom_density(alpha=1) 

p

png("histo_reads_ratios_to_5.5.png", width = 1200, height = 800, res = 120)
plot(p)
dev.off()

#savePlot(filename="", type = "png")
