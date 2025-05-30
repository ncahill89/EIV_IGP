###Clear workspace###
rm (list = ls( )) 

### Required Packages and Libraries ###
## uncomment the install.packages if these packages aren't already installed ###
#install.packages(c("rjags","R2jags","fields","tidyverse))
library(rjags)
library(R2jags)
library(fields)
library(tidyverse)
##############################################################

### Load the required R functions ###
Rfiles <- list.files(file.path(getwd(), "R"))
Rfiles <- Rfiles[grepl(".R", Rfiles)]
sapply(paste0("R/", Rfiles), source)

### An example of the correct format for the data is shown for New York
## Do not change the headings (RSL,RSLError,Age,AgeError)
## Errors must be 1-sigma
## If Age is given in the BP scale, change the BP_age_scale argument to TRUE in the dataprep() function below

### Prepare the data for the model ###
## User to input the data file location and the name of their dataset e.g., dataname = "New York"
## Plot of the raw data will be saved to fig folder
data.raw <- dataprep(data_path = "data/NYC.csv",
                   dataname = "New York",
                   BP_age_scale = TRUE)

## NOTE!!! If correcting for GIA change GIA = TRUE and provide a rate of GIA in mm/yr

### Run the model ###
## User to provide the interval for how often they want predictions from the model
interval = 50

## run
RunIGPModel(data.raw=data.raw,
            interval = interval)

## Check convergence for GP parameters
## Note: this won't work if fast = TRUE
get_diagnostics(data.raw=data.raw)

### Output the model figures and results ###
## Figures will output to fig folder
## .csv file will output to results folder
IGPResults(data.raw = data.raw,
           interval = interval)

