#!/bin/python

import sys
import getopt
import csv

#IMPORT the combined batch csv

from pathlib import Path
Path('/apittman-aramis/sgul/shares/aramis/Genetics_Centre_Bioinformatics/CNV_SGUL_Exomes/SQL/CNV_data.db').touch()

#create sqlite3 table
import sqlite3
conn = sqlite3.connect('/apittman-aramis/sgul/shares/aramis/Genetics_Centre_Bioinformatics/CNV_SGUL_Exomes/SQL/CNV_data.db')
c = conn.cursor()

#drop old table
c.execute('DROP TABLE frequencies')

#table with following columns

create_temp = ["CREATE TABLE 'frequencies' ('CNV_ID' text, 'lowest_frequency' text, 'number_CNVs_at_lowest_exon' text, 'highest_frequency' text, 'number_CNVs_at_highest_exon' text, 'mean_frequency' text, 'mean_number_of_exons' text)"]
create_table = ''
create_table = create_table.join(create_temp)

c.execute(create_table)

csv_comb_all_path = '/apittman-aramis/sgul/shares/aramis/Genetics_Centre_Bioinformatics/CNV_SGUL_Exomes/frequency_CNVs/all_freq_combined.csv'

#load CSV into sqlite table
import pandas as pd
#load data into Pandas dataframe
freq_table = pd.read_csv(csv_comb_all_path)     #, sep=',\s+', quoting=csv.QUOTE_ALL
#select wanted columns - to prevent errors
#all_comb_table_safe = all_comb_table.loc[:, 'batch_ID':'exons.hg19']
#write data to sql table
freq_table.to_sql('frequencies', conn, if_exists='replace', index = False)

conn.commit()
conn.close()