# MEGP_CNV_analysis

The Middle Eastern Genomes Project: Copy Number Variant Analysis to Solve Rare Genetic Diseases

Authors: Luke Peacock, Reza Maroofian, Ehsan Karimiani, Paria Najarzadeh, Chris Carroll, Yalda Jamshidi, Alan M Pittman

The Middle Eastern Genomes Project (MEGP) aims to investigate the genetic cause of disease in Middle-Eastern families with rare inherited conditions, including neurodevelopmental and musculoskeletal.
This cohort of nearly 2000 individuals is enriched for autozygosity, due to a higher incidence of consanguinity in the families under study, and provides a unique opportunity to investigate rare Mendelian disorders.

Aims:
-	To develop a bioinformatic framework for identifying Copy Number Variants (CNVs) in whole exome sequencing data.
-	To identify disease-causing CNVs in the Middle Eastern Genomes Project (MEGP) cohort.
-	To construct an interactive browser application to enhance collaboration between MEGP researchers and facilitate identification of genetic causes of rare disease.

hg37-aligned BAM files were analysed using a novel application of ExomeDepth package; CNVs were detected by comparing samples within matched sequencing batches.
To improve recognition of disease-causing CNVs, Regions-of-Homozygosity (ROHs) were calculated for all samples and CNVs were annotated with known disease associations using the Online Mendelian Inheritance in Man (OMIM) database of genotype-phenotype correlations.
There is limited publicly available CNV data for the Middle-Eastern population, with current databases heavily biased towards individuals of European descent.
To address these ethnic disparities, we built an in-house CNV frequency database unique to the Middle-Eastern population.
We employed a SQL relational database to combine our CNV findings with ROH data, OMIM phenotype annotations and the Middle-Eastern CNV frequency database.

We constructed a web application to combine the CNV findings with ROH, OMIM and frequency data.
The application, which utilises a server-side SQL database, allows users to select from various parameters, view a sample output and download the desired data in csv format.
We successfully identified two novel disease-causing homozygous deletions; a deletion in SEPN causing a congenital myopathy and a deletion in GALC causing a paediatric white-matter neurological disorder.

The initial identification of two CNVs shows this technique’s power in identifying disease-causing homozygous recessive deletions and demonstrates the research value of this cohort. We will continue to search for further findings in this dataset.
We aim to incorporate this technique into routine exome analysis to aid identification of disease-causing CNVs.

Available at: https://lpeacock.shinyapps.io/CNV_dev_version/

More information: https://www.genomemed.org/megp
