#!/bin/Rscript

library(readxl)
library(tidyverse)
library(DBI)
library(RSQLite)

path <- "/apittman-aramis/sgul/shares/aramis/Genetics_Centre_Bioinformatics/CNV_SGUL_Exomes/SGUL_GeneticsCentre_Exome_Sample_Inventory.xlsx"

inventory_df <- path %>%
  readxl::excel_sheets() %>%
  purrr::set_names() %>%
  purrr::map_df(~ read_excel(path = path, sheet = .x, range = NULL), .id = "sheet")
print(inventory_df, n = Inf)

inv_slim <- dplyr::select(inventory_df, "sheet", "PI", "projectID", "sampleID", "DRIVE", "BATCH_ID", "Run_Date")

conn <- dbConnect(RSQLite::SQLite(), '/apittman-aramis/sgul/shares/aramis/Genetics_Centre_Bioinformatics/CNV_SGUL_Exomes/SQL/CNV_data.db')

fields <- c(sheet="text", PI="text", projectID="text", sampleID="text", DRIVE="text", BATCH_ID="text", Run_Date="TEXT")

DBI::dbWriteTable(conn, "sample_info", inv_slim, overwrite = TRUE, field.types = fields)