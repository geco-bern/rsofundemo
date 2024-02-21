# create validtion file used in p_model
library(dplyr)
library(tidyr)
library(ggplot2)
library(reshape2)
library(lubridate)
library(FluxDataKit)
library(rsofun)
library(here)
source("./R/read_meta_fdk.R")
source("./R/generation_validation_bysitename.R")

# insert all the files in the data-raw folder !!!NEED ALSO valid_years.csv

sites <- c("ES-Amo","FR-Pue")

csv_path <- here("data-raw//")
lsm_path <- here("data-raw//")

# creation of validation file, the start data are obtained with the repo https://github.com/geco-bern/FluxDataKit

validation = lapply(sites,function(x){create_validation(x,csv_path,lsm_path)})
validation <- dplyr::bind_rows(validation)

# save validation file, used in vignette
saveRDS(
  validation,
  paste0(here("data//"),"rsofun_validation_data.rds"),
  compress = "xz"
)
