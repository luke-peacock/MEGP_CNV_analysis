#!/bin/python

import sys
import getopt

#IMPORT the combined batch csv

from pathlib import Path
Path('/apittman-aramis/sgul/shares/aramis/Genetics_Centre_Bioinformatics/CNV_SGUL_Exomes/SQL/CNV_data.db').touch()

#create sqlite3 table
import sqlite3
conn = sqlite3.connect('/apittman-aramis/sgul/shares/aramis/Genetics_Centre_Bioinformatics/CNV_SGUL_Exomes/SQL/CNV_data.db')
c = conn.cursor()

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
