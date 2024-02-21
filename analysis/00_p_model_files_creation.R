# create drivcers and validtion files used in p_model
library(dplyr)
library(tidyr)
library(ggplot2)
library(reshape2)
library(lubridate)
library(FluxDataKit)
library(ncdf4)
library(rsofun)
source("./R/read_meta_fdk.R")
source("./R/generation_validation_bysitename.R")

# select your own paths, the outpath doesn't need to exist since it will be created
# insert all the files in the raw_data folder

# xxx put external data into data-raw/

# xxx provide info for where the external data comes from (e.g., point to script in external repository, or DOI for data download, ...)

site <- "ES-Amo"

# use here::here("data-raw/")
csv_path <- paste0(getwd(),"/raw/")
lsm_path <- paste0(getwd(),"/raw/")
nc_path <- paste0(getwd(),"/raw/")
out_path <- paste0(getwd(),"/data/")

# running once per site since the results are saved
create_driver_validation(site,csv_path,lsm_path,nc_path,out_path)

# xxx split creation of forcing (driver) and validation
