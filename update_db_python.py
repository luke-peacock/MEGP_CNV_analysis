#!/bin/python

#####	python import csvs #####

import sys
import getopt

#IMPORT the combined batch csv

from pathlib import Path
Path('/apittman-aramis/sgul/shares/aramis/Genetics_Centre_Bioinformatics/CNV_SGUL_Exomes/SQL/CNV_data.db').touch()

#create sqlite3 table
import sqlite3
conn = sqlite3.connect('/apittman-aramis/sgul/shares/aramis/Genetics_Centre_Bioinformatics/CNV_SGUL_Exomes/SQL/CNV_data.db')
c = conn.cursor()

############################################        full_combined          #############################################

#drop old table
c.execute("DROP TABLE fullcombined")

#table with following columns

create_temp = ["CREATE TABLE 'fullcombined' ('batch_ID' text, 'ID' text, 'start.p' numeric, 'end.p' numeric, 'type' text, 'nexons' numeric, 'start' numeric, 'end' numeric, 'chromosome' numeric, 'CNV_id' text, 'BF' numeric, 'reads.expected' numeric, 'reads.observed' numeric, 'reads.ratio' numeric, 'Conrad.hg19' text, 'exons.hg19' text)"]
create_table = ""
create_table = create_table.join(create_temp)

c.execute(create_table)

csv_comb_all_path = '/apittman-aramis/sgul/shares/aramis/Genetics_Centre_Bioinformatics/CNV_SGUL_Exomes/SQL/combined_all_csv_full_ID.csv'

#load CSV into sqlite table
import pandas as pd
#load data into Pandas dataframe
all_comb_table = pd.read_csv(csv_comb_all_path)     #, sep=',\s+', quoting=csv.QUOTE_ALL
#select wanted columns - to prevent errors
#all_comb_table_safe = all_comb_table.loc[:, 'batch_ID':'exons.hg19']
#write data to sql table
all_comb_table.to_sql('fullcombined', conn, if_exists='replace', index = False)

#rename column
rename = ["ALTER TABLE fullcombined RENAME COLUMN ID_CNV_id_type TO CNV_id;"]
#c.execute(rename)

############################################        roh_combined          #############################################
#drop old table
c.execute("DROP TABLE rohcombined")

#table with following columns

create_temp = ["CREATE TABLE 'rohcombined' ('ROH_chromosome' text, 'ROH_start' text, 'ROH_end' text, 'ROH_sample_ID' text, 'ROH_length' numeric, 'ROH_number_of_markers' numeric, 'ROH_quality_score' numeric, 'CNV_chromosome' text, 'CNV_start' text, 'CNV_end' text, 'CNV_ID' text, 'batch_ID' text, 'sample_ID' text, 'start_p' numeric, 'end_p' numeric, 'CNV_type' text, 'nexons' numeric, 'BF' numeric, 'reads_expected' numeric, 'reads_observed' numeric, 'reads_ratio' numeric, 'Conrad_hg19' text, 'exons_hg19' text)"]

create_table = ""
create_table = create_table.join(create_temp)

c.execute(create_table)

csv_roh_all_path = '/apittman-aramis/sgul/shares/aramis/Genetics_Centre_Bioinformatics/CNV_SGUL_Exomes/SQL/combined_roh_calls_2.csv'

#load CSV into sqlite table
import pandas as pd
#load data into Pandas dataframe
all_comb_table = pd.read_csv(csv_roh_all_path)
#write data to sql table
all_comb_table.to_sql('rohcombined', conn, if_exists='replace', index = False)

#rename column
rename = ["ALTER TABLE rohcombined RENAME COLUMN ROH_sample_ID_CNV_ID_CNV_type TO CNV_id;"]
#c.execute(rename)

############################################        OMIM_full          #############################################
#drop old table
c.execute("DROP TABLE OMIM_full")

#table with following columns

create_temp = ["CREATE TABLE 'OMIM_full' ('batch_ID' text, 'sample_ID' text, 'start.p' numeric, 'end.p' numeric, 'CNV_type' text, 'nexons' numeric, 'CNV_start' numeric, 'CNV_end' numeric, 'CNV_chromosome' text, 'CNV_ID' text, 'BF' numeric, 'reads.expected' numeric, 'reads.observed' numeric, 'reads.ratio' numeric, 'Conrad.hg19' text, 'exons.hg19' text, 'chrom_OMIM' text, 'OMIM_start' numeric, 'OMIM_end' numeric, 'OMIM_gene_name' text, 'OMIM_phenotype' text, 'OMIM_bp_overlap' numeric)"]

create_table = ""
create_table = create_table.join(create_temp)

c.execute(create_table)

csv_path = '/apittman-aramis/sgul/shares/aramis/Genetics_Centre_Bioinformatics/CNV_SGUL_Exomes/SQL/combined_OMIM_main.csv'

#load CSV into sqlite table
import pandas as pd
#load data into Pandas dataframe
comb_table = pd.read_csv(csv_path)
#write data to sql table
comb_table.to_sql('OMIM_full', conn, if_exists='replace', index = False)

############################################        OMIM_roh          #############################################
#drop old table
c.execute("DROP TABLE OMIM_roh")

#table with following columns

create_temp = ["CREATE TABLE 'OMIM_roh' ('batch_ID' text, 'sample_ID' text, 'start.p' numeric, 'end.p' numeric, 'CNV_type' text, 'ROH_chromosome' text, 'ROH_start' numeric, 'ROH_end' numeric, 'ROH_sample_ID' text, 'ROH_length' numeric, 'ROH_number_of_markers' numeric, 'ROH_quality_score' numeric, 'nexons' numeric, 'CNV_start' numeric, 'CNV_end' numeric, 'CNV_chromosome' text, 'CNV_ID' text, 'BF' numeric, 'reads.expected' numeric, 'reads.observed' numeric, 'reads.ratio' numeric, 'Conrad.hg19' text, 'exons.hg19' text, 'roh_bp_overlap' numeric, 'chrom_OMIM' text, 'OMIM_start' numeric, 'OMIM_end' numeric, 'OMIM_gene_name' text, 'OMIM_phenotype' text, 'OMIM_bp_overlap' numeric)"]

create_table = ""
create_table = create_table.join(create_temp)

c.execute(create_table)

csv_path = '/apittman-aramis/sgul/shares/aramis/Genetics_Centre_Bioinformatics/CNV_SGUL_Exomes/SQL/combined_OMIM_roh.csv'

#load CSV into sqlite table
import pandas as pd
#load data into Pandas dataframe
comb_table = pd.read_csv(csv_path)
#write data to sql table
comb_table.to_sql('OMIM_roh', conn, if_exists='replace', index = False)

############################################        import frequencies          #############################################

from pathlib import Path
Path('/apittman-aramis/sgul/shares/aramis/Genetics_Centre_Bioinformatics/CNV_SGUL_Exomes/SQL/CNV_data.db').touch()

#create sqlite3 table
import sqlite3
conn = sqlite3.connect('/apittman-aramis/sgul/shares/aramis/Genetics_Centre_Bioinformatics/CNV_SGUL_Exomes/SQL/CNV_data.db')
c = conn.cursor()

#drop old table
c.execute("DROP TABLE frequencies")

#table with following columns

create_temp = ["CREATE TABLE 'frequencies' ('CNV_ID' text, 'lowest_frequency_exon' text, 'number_CNVs_at_exon' text)"]
create_table = ""
create_table = create_table.join(create_temp)

c.execute(create_table)

csv_comb_all_path = '/apittman-aramis/sgul/shares/aramis/Genetics_Centre_Bioinformatics/CNV_SGUL_Exomes/frequency_CNVs/all_freq_combined_NAs.csv'

#load CSV into sqlite table
import pandas as pd
#load data into Pandas dataframe
freq_table = pd.read_csv(csv_comb_all_path)     #, sep=',\s+', quoting=csv.QUOTE_ALL
#select wanted columns - to prevent errors
#all_comb_table_safe = all_comb_table.loc[:, 'batch_ID':'exons.hg19']
#write data to sql table
freq_table.to_sql('frequencies', conn, if_exists='replace', index = False)

############
#close connection
conn.commit()
conn.close()

