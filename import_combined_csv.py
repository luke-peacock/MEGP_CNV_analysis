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

conn.commit()
conn.close()