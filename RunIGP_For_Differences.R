### Clear workspace ###
rm(list = ls())

### Required Packages ###
library(rjags)
library(R2jags)
library(fields)
library(tidyverse)
library(readxl)

##############################################################

### Load required R functions ###
Rfiles <- list.files(file.path(getwd(), "R"))
Rfiles <- Rfiles[grepl(".R", Rfiles)]
sapply(paste0("R/", Rfiles), source)

##############################################################

### Model settings ###
interval <- 50

# orig data is BP age scale 
lwr <- 4450
upr <- 900

excel_file <- "data/tests.xlsx"
sheet_names <- excel_sheets(excel_file)

##############################################################

### Loop over Excel sheets ###
for (sheet in sheet_names) {
  
  message("Running model for dataset: ", sheet)
  
  ## Read sheet
  df <- read_excel(excel_file, sheet = sheet)
  
  ## Write temporary CSV (required by dataprep)
  temp_csv <- file.path(tempdir(), paste0(sheet, ".csv"))
  write.csv(df, temp_csv, row.names = FALSE)
  
  ## Prepare data
  data.raw <- dataprep(
    data_path = temp_csv,
    dataname  = sheet,
    BP_age_scale = TRUE
  )
  
  ## Run model
  RunIGPModel(
    data.raw = data.raw,
    interval = interval,
    EIV = TRUE,
    fast = FALSE,
    upper = upr,
    lower = lwr,
    incl.errbounds = FALSE
  )
  
  ## Diagnostics
  get_diagnostics(data.raw = data.raw)
  
  ## Output results
  IGPResults(
    data.raw = data.raw,
    interval = interval,
    upper = upr,
    lower = lwr
  )
}

