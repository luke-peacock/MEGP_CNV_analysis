
library(ExomeDepth)
library(GenomeInfoDb)
library(Rsamtools)


## Import Modified hg19 exons definitions
#  These have ~5% of the total exons removed due to poor mappability and repetative sequences


exons.hg19.mod <- as.data.frame(read.csv2(file = "/apittman-aramis/sgul/shares/aramis/Genetics_Centre_Bioinformatics/CNV_SGUL_Exomes/ExomeDepth_Exons_hg19_modified.csv", sep = ","))


# Read in the config.csv file containg our bam sample list

args <- commandArgs(trailingOnly = FALSE)
config_csv <- args[6]
#print(args)
#print(config_csv)

analysisConfig <- read.csv(paste("/apittman-aramis/sgul/shares/aramis/Genetics_Centre_Bioinformatics/CNV_SGUL_Exomes/", config_csv, sep = ""), 
							              header = TRUE, 
							              fill = TRUE)

list_of_bam_files <- as.vector(analysisConfig$list_of_bam_files)

print(list_of_bam_files)

sample_path <- paste("/apittman-aramis/sgul/shares/aramis/Genetics_Centre_Bioinformatics/CNV_SGUL_Exomes/hg38_BAMs_to_convert/", list_of_bam_files, "/", list_of_bam_files, "_hg37aligned_sorted_rmduplicates.bam", sep="")

print(sample_path)

# Perform read counting of exonic intervals - modified exon list


my.counts <- getBamCounts(bed.frame = exons.hg19.mod,
                          bam.files = sample_path,
                          include.chr = FALSE,)

print(head(my.counts))


# Create a dataframe and also a matrix of the exonic read counts


ExomeCount.dafr <- as(my.counts[, colnames(my.counts)], 'data.frame')

#check the fields and columns in dataframe:
print(head(ExomeCount.dafr))

# Create matrix of the bam counts only
ExomeCount.mat <- as.matrix(ExomeCount.dafr[, grep(names(ExomeCount.dafr), 
							pattern = '*.bam')])

##check the fields and columns in matrix:
print(head(ExomeCount.mat))



# Perform ExomeDepth CNV calling on each sample using a loop


nsamples <- ncol(ExomeCount.mat)

print(head(nsamples))


message('Now looping over all the samples innit')


for (i in 1:nsamples) {


my.test.data <- as.matrix(ExomeCount.mat[, i])

my.reference.set <- as.matrix(ExomeCount.mat[, -i])

my.choice <- select.reference.set(test.counts = my.test.data,
                                    reference.counts = my.reference.set,
                                    bin.length = (ExomeCount.dafr$end - ExomeCount.dafr$start)/1000,
                                    n.bins.reduced = 10000)


head(my.choice)

my.matrix <- as.matrix( ExomeCount.dafr[, my.choice$reference.choice, drop = FALSE])

my.reference.set <- apply(X = my.matrix,
                                MAR = 1,
                                FUN = sum)


## CNV calling

all.exons <- new('ExomeDepth',
                  test = ExomeCount.mat[,i],
                  reference = my.reference.set,
                  formula = 'cbind(test, reference) ~ 1')

#We can now call the CNV by running the underlying hidden Markov model:

all.exons <- CallCNVs(x = all.exons,
                      transition.probability = 10^-4,
                      chromosome = ExomeCount.dafr$chromosome,
                      start = ExomeCount.dafr$start,
                      end = ExomeCount.dafr$end,
                      name = ExomeCount.dafr$exon)

#check output
head(all.exons@CNV.calls)

#Annotating with Conrad Common CNVs

data(Conrad.hg19)

head(Conrad.hg19.common.CNVs)

all.exons <- AnnotateExtra(x = all.exons,
                          reference.annotation = Conrad.hg19.common.CNVs,
                          min.overlap = 0.5,
                          column.name = 'Conrad.hg19')



print(head(all.exons@CNV.calls))

#now annotating with exon/gene level information. 

data(exons.hg19)
 
exons.hg19.GRanges <- GenomicRanges::GRanges(seqnames = exons.hg19.mod$chromosome, 
											IRanges::IRanges(start=exons.hg19.mod$start,end=exons.hg19.mod$end), 
											names = exons.hg19.mod$name)

#here the minimum overlap should be very close to 0  
all.exons <- AnnotateExtra(x = all.exons, 
						reference.annotation = exons.hg19.GRanges, 
						min.overlap = 0.0001, 
						column.name = 'exons.hg19')

#Now save it in an easily readable format

output.file <- paste("/apittman-aramis/sgul/shares/aramis/Genetics_Centre_Bioinformatics/CNV_SGUL_Exomes/CNV_calls/", list_of_bam_files[i], 
					          '.csv', sep = '')

write.csv(file = output.file, 
		      x = all.exons@CNV.calls, 
		      row.names = FALSE)

}












