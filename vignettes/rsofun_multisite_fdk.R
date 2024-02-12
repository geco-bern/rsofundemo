# xxx work with library(here), see AGDS book
# params <-
# list(output_dir = ".")

#' ---
#' title: "Generate input data for rsofun Phydro"
#' author: "Francesco Grossi"
#' date: "2024_02_09"
#' output: html_document
#' params:
#'   output_dir: "."
#' ---
#'
## ----setup, include=FALSE--------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

#'
## --------------------------------------------------------------------
library(tidyverse)
library(reshape2) # recommended to use dplyr (you've loaded tidyverse already)
#library(FluxDataKit)
library(lubridate)
library(ncdf4)

# xxx cannot read this file - it's not part of the repo.
# Either source from the same repo with source(here::here("R/function_name.R")),
# or load external functions from libraries. Here, I guess this is from the external
# library FluxDataKit. Load it.
library(FluxDataKit)
# source("read_meta_fdk.R")

## --------------------------------------------------------------------
# xxx this statement is only used if the script is run from the shell (terminal)
# with arguments provided from the shell
# args = commandArgs(trailingOnly=TRUE)

sites <- c("GF-Guy","BE-Vie")

# xxx highlight that this is user specific - data is external to this project
# root_data_dir <- "./fluxdatakit_oct3"
root_data_dir <- "/Users/benjaminstocker/data/FluxDataKit"

## --------------------------------------------------------------------

lsm_path = paste0(root_data_dir, "/FLUXDATAKIT_LSM/")
csv_path = paste0(root_data_dir, "/FLUXDATAKIT_FLUXNET/")
out_path = paste0(root_data_dir, "/phydro_drivers/")
# xxx just to clarify: we're not running phydro, phydro uses additional forcing time series which are not relevat for our modelling. our workflow should thus be simpler.

# xxx work with here::here()
dir.create(out_path, showWarnings = F)
# dir.create(paste0(out_path, "data_gen_figures/"), showWarnings = F)
dir.create(here::here("data_gen_figures"), showWarnings = F)

#figures_prefix = paste0(out_path,"/data_gen_figures/", site)

files_csv = list.files(csv_path)
files_lsm = list.files(lsm_path)


csv_acr <- sapply((substr(files_csv,5,10)), function(x) x[1])

file_csv <- files_csv[csv_acr %in% sites]
file_csv <- file_csv[grep("HH", file_csv)]

# get metadata
# --------------------------------------------------------
message("- reading Metadata for site")

tmp_fun <- function(x){
  read_meta_fdk(
    site = x,
    path = lsm_path,
    meta_data = T)}

meta <- lapply(sites,tmp_fun)
names(meta) <- sites
# get half-hourly data  --------------------------------------------------------
# message("- convert to FLUXNET standard CSV file")

message("- reading FLUXNET format halfhourly data")

backup <- getwd()
setwd(csv_path)
hhdf <- lapply(file_csv, read_csv)
setwd(backup)
names(hhdf) <- sites


# Add date and time columns to hhdf for easier further processing.
# ---------------------------------------------------------
message("- Add date and time columns")
message("- Add SW_OUT=NA if not present")


tmp_fun <- function(x){
  x <- x|>
    mutate(time = lubridate::as_datetime(as.character(TIMESTAMP_START), tz = "GMT", format="%Y%m%d%H%M")) |>
    mutate(date = lubridate::as_date(time))
  if (!("SW_OUT" %in% colnames(x))) {
    x$SW_OUT = NA}
  return(x)
}
hhdf <- lapply(hhdf,tmp_fun)

# Aggregate to daily 24-hr means  ----------------------------------------------------------
message("- downsampling FLUXNET format - 24 hr means")

tmp_fun <-  function(x){
  result <- hhdf[names(hhdf)==x][[1]] |>
    group_by(date) |>
    select(-TIMESTAMP_START, -TIMESTAMP_END) |>
    summarize_all(.funs = mean)
  return(result)
}

ddf_24hr_mean <- lapply(sites,tmp_fun)
names(ddf_24hr_mean) <- sites
# P1 check

# xxx always use here::here() to make sure the code runs irrespective of where the script file is saved
valid_years <- read.csv(here::here("vignettes/ancillary_data/valid_years_final.csv"), header = T, fileEncoding = "UTF-16")

tmp_fun <- function(x){
  ystart = valid_years %>% filter(Site==x) %>% pull(start_year)
  yend = valid_years %>% filter(Site==x) %>% pull(end_year)
  yspan <- data.frame(ystart=ystart,yend=yend,site = x)
  return(yspan)
}

year_span <- sapply(sites,tmp_fun)
year_span <- data.frame(t(year_span))
rownames(year_span) = NULL

# Aggregate around daily maximum ppfd for acclimating model
# ---------------------------------------------------------

tmp_fun <- function(x){
  ystart <-  year_span %>% filter(site==x) %>% pull(ystart)
  yend <-  year_span %>% filter(site==x) %>% pull(yend)

  result <- hhdf[names(hhdf)==x][[1]] |> filter(date >= as_date(paste0(floor((ystart[[1]]+yend[[1]])/2),"-06-01")) &
                                                  date <= as_date(paste0(floor((ystart[[1]]+yend[[1]])/2),"-06-03")) )
  return(result)
}


test_3day <- lapply(sites,tmp_fun)
names(test_3day) <- sites


aggregate_daily_3hr_maxima = function(df){
  # Get the time at which SW_IN is maximum
  maxppfd <- df |>
    filter(SW_IN_F_MDS == max(SW_IN_F_MDS))
  max_t <- maxppfd$time[1]

  # Select times that lie in 3-hr interval around max_t
  df_aroundmax <- df |>
    filter(time < (max_t + 1.5*3600) & time > (max_t - 1.5*3600) )

  # take mean of selected entries
  df_mean <- df_aroundmax |>
    select(-TIMESTAMP_START, -TIMESTAMP_END) |>
    summarize_all(.funs = mean)

  df_mean
}

# Test aggregation
# ----------------
# TO DO (not work)

tmp_fun <- function(x){
  result <- test3day[names(test3day)==x][[1]] |>
    mutate(across(ends_with("date"), as.Date, format = "%m/%d/%Y")) |>
    group_by(date)|>
    do(aggregate_daily_3hr_maxima(.)) |>
    ungroup()
  return(result)
}

test_3day_3hr <- lapply(test_3day, tmp_fun)
#test.3day.3hr = test.3day %>% group_by(date) %>% do(aggregate_daily_3hr_maxima(.)) %>% ungroup()

#P2 test

## --------------------------------------------------------------------
# Apply 3hr maxima aggregation to all data
# ----------------------------------------
message("- downsampling FLUXNET format - daily 3-hr means around max ppfd")
tmpfun <- function(x){
  result <- hhdf[names(hhdf)==x] [[1]]|>
    group_by(date) |>
    do(aggregate_daily_3hr_maxima(.)) |>
    ungroup()
  return(result)
}


ddf_3hr_maxima <- lapply(sites,tmp_fun)
names(ddf_3hr_maxima) <- sites

#'
## --------------------------------------------------------------------
aggregate_daily_daylength = function(df){
  # Get the time at which SW_IN > 0
  pos_ppfd <- df %>% filter(SW_IN_F_MDS > 10)
  # if SW_IN is unavailable in that year calc daylength based on NETRAD
  if (nrow(pos_ppfd) < 2){
    pos_ppfd <- df %>% filter(NETRAD > 25)
  }

  tmax <- max(pos_ppfd$time)
  tmin <- min(pos_ppfd$time)

  # Select times that lie in 3-hr interval around max_t
  df_aroundmax <- df %>% filter(time <= tmax &
                                time >= tmin )

  # take mean of selected entries
  df_mean <- df_aroundmax |>
    select(-TIMESTAMP_START, -TIMESTAMP_END) |>
    summarize_all(.funs = mean) |>
    mutate(daylength = difftime(tmax, tmin, units="hours") |> as.numeric())

  df_mean
}

# Test aggregation
# ----------------

tmp_fun <- function(x){
  result <- test3day[names(test3day)==x] [[1]]|>
    group_by(date) |>
    do(aggregate_daily_daylength(.))|>
    ungroup()
  return(result)
}

test3daydaylen <- lapply(sites,tmp_fun)

#P3 test
## --------------------------------------------------------------------
# Apply daytime mean aggregation to all data
# ------------------------------------------
message("- downsampling FLUXNET format - daytime means")

tmp_fun <- function(x){
  result <- hhdf[names(hhdf)==x][[1]]|>
    group_by(date) |>
    do(aggregate_daily_daylength(.)) |>
    ungroup()
  return(result)
}

ddf_daytime_mean <- lapply(sites,tmp_fun)
names(ddf_daytime_mean) <- sites

#P4 test
## --------------------------------------------------------------------
# Calculate daily tmax and tmin from hh data
# ------------------------------------------

tmp_fun <- function(x){
  result <- hhdf[names(hhdf)==x][[1]]|>
    group_by(date) |>
    summarize(
      tmax = max(TA_F_MDS),
      tmin = min(TA_F_MDS)
    )
  return(result)
}

tmaxmin <- lapply(sites,tmp_fun)
names(tmaxmin) <- sites

## --------------------------------------------------------------------
# Creating driver object  ------------------------------------------------------

# common part
message("- compiling drivers")
load(here::here("data/p_model_drivers.rda"))  # xxx added here::here()
nc <- nc_open(here::here("vignettes/ancillary_data/cwdx80.nc"))   # xxx added here::here()
lons <- ncvar_get(nc, "lon")
lats <- ncvar_get(nc, "lat")
S80 <- ncvar_get(nc, "cwdx80")
p_hydro_drivers <- p_model_drivers[-1]
kfFEC <- 2.04
n <- 1


# nested list creation
p_hydro_drivers_list <- list()
tmp_fun <- function(x){
  p_hydro_drivers_list <- append(p_hydro_drivers_list,p_hydro_drivers)
}

p_hydro_drivers_list <- lapply(sites,tmp_fun)
names(p_hydro_drivers_list) <- sites

# function 1
#big function until a df is created

tmp_fun <- function(x){

  # variable creation
  site_lon <- meta[names(meta)==x][[1]][[1]]
  site_lon <- site_lon$longitude
  site_lat <- meta[names(meta)==x][[1]][[1]]
  site_lat <- site_lat$latitude
  lonid <- which(lons > site_lon)[1]-1
  latid <- which(lats > site_lat)[1]-1
  S80_slice <- S80[(lonid-n):(lonid+n), (latid-n):(latid+n)]
  whc_site <- mean(as.numeric(S80_slice, na.rm=T))
  whc_site_sd <- sd(as.numeric(S80_slice, na.rm=T))

  # df selection
  result <- p_hydro_drivers_list[names(p_hydro_drivers_list)==x][[1]]
  # data creation to insert in site info
  lon <- meta[names(meta)==x][[1]][[1]]
  lon <- lon$longitude
  lat <- meta[names(meta)==x][[1]][[1]]
  lat <- lat$latitude
  elv <- meta[names(meta)==x][[1]][[1]]
  elv <- elv$elevation
  IGBP_veg_short <- meta[names(meta)==x][[1]][[1]]
  IGBP_veg_short <- IGBP_veg_short$IGBP_veg_short

  canpoy_data <-  meta[names(meta)==x][[1]][[1]]
  canpoy_data <- canpoy_data$canopy_height
  canopy_height <- ifelse(is.na(canpoy_data), yes = 20, canpoy_data)

  reference_data <-  meta[names(meta)==x][[1]][[1]]
  reference_data <- reference_data$reference_height
  reference_height <- ifelse(is.na(reference_data), yes = 22, reference_data)

  # tibble creation

  result$site_info[[1]] =
    tibble(
      lon=lon,
      lat=lat,
      elv = elv,
      canopy_height = canopy_height,
      reference_height = reference_height,
      whc = whc_site,
      whc_sd = whc_site_sd,
      IGBP_veg_short = IGBP_veg_short
    )
  return(result)
}

p_hydro_drivers_list<- lapply(sites,tmp_fun)
names(p_hydro_drivers_list) <- sites

# function 2

tmp_fun <- function(x){

  #data creation
  start_year <- valid_years %>% filter(Site==x) %>% pull(start_year)
  end_year <- valid_years %>% filter(Site==x) %>% pull(end_year)
  # tibble selection
  ddf_24hr_mean <-  ddf_24hr_mean[names(ddf_24hr_mean)==x][[1]]
  # tibble creation
  result <- p_hydro_drivers_list[names(p_hydro_drivers_list)==x][[1]]
  result$forcing <- ddf_24hr_mean |>
    dplyr::filter(lubridate::year(date) %in% start_year:end_year) |>
    dplyr::filter(!(lubridate::mday(date) == 29 & lubridate::month(date) == 2)) |>
    left_join(tmaxmin[names(tmaxmin)==x][[1]]) |>
    group_by(date) |>
    summarize(
      date = date,
      temp = TA_F_MDS,
      vpd = VPD_F_MDS * 100,
      ppfd = SW_IN_F_MDS * kfFEC * 1e-06,
      netrad = NETRAD,
      patm = PA_F * 1000,
      snow = 0,
      rain = P_F * 48 /(60 * 60 * 24),
      tmin = tmaxmin[names(tmaxmin)==sites[1]][[1]]$tmin, # TMIN_F_MDS,
      tmax = tmaxmin[names(tmaxmin)==sites[1]][[1]]$tmax, # TMAX_F_MDS,
      fapar = FPAR,
      co2 = CO2_F_MDS,
      ccov = 0
    ) |>
    list()
  return(result)
}

p_hydro_drivers_list<- lapply(sites,tmp_fun)
names(p_hydro_drivers_list) <- sites


# function 3
# TO DO (not work)

tmp_fun <- function(x){

  #data creation
  start_year <- valid_years %>% filter(Site==x) %>% pull(start_year)
  end_year <- valid_years %>% filter(Site==x) %>% pull(end_year)
  # tibble selection
  ddf_3hr_maxima <-  ddf_3hr_maxima[names(ddf_3hr_maxima)==x][[1]]
  # tibble creation

  result <- p_hydro_drivers_list[names(p_hydro_drivers_list)==x][[1]]
  result$forcing_acclimforcing_acclim <- ddf_3hr_maxima|>
    dplyr::filter(lubridate::year(date) %in% start_year:end_year) |>
    dplyr::filter(!(lubridate::mday(date) == 29 & lubridate::month(date) == 2)) |>
    left_join(tmaxmin[names(tmaxmin)==x][[1]]) |>
    group_by(date) |>
    summarize(
      date = date,
      time = time,
      temp = TA_F_MDS,
      vpd = VPD_F_MDS * 100,
      ppfd = SW_IN_F_MDS * kfFEC * 1e-06,
      netrad = NETRAD,
      patm = PA_F * 1000,
      snow = 0,
      rain = NA, # P_F * 48 / (60 * 60 * 24),
      tmin =  tmaxmin[names(tmaxmin)==sites[1]][[1]]$tmin,
      tmax = tmaxmin[names(tmaxmin)==sites[1]][[1]]$tmax,
      fapar = FPAR,
      co2 = CO2_F_MDS,
      ccov = 0
    ) |>
    list()
}

#I use this notation (and store multiple element) due to problem in function 3, but in this way I'm unable to save
p_hydro_drivers_list_2<- lapply(sites,tmp_fun)
names(p_hydro_drivers_list_2) <- sites


#function 4

tmp_fun <- function(x){

  #data creation
  start_year <- valid_years %>% filter(Site==x) %>% pull(start_year)
  end_year <- valid_years %>% filter(Site==x) %>% pull(end_year)
  # tibble selection
  ddf_daytime_mean <-  ddf_daytime_mean[names(ddf_daytime_mean)==x][[1]]
  # tibble creation
  result <- p_hydro_drivers_list[names(p_hydro_drivers_list)==x][[1]]
  result$forcing_daytime_mean <- ddf_daytime_mean |>
  dplyr::filter(lubridate::year(date) %in% start_year:end_year) |>
  dplyr::filter(!(lubridate::mday(date) == 29 & lubridate::month(date) == 2)) |>
  left_join(tmaxmin[names(tmaxmin)==x][[1]]) |>
  group_by(date) |>
  summarize(
    date = date,
    time = time,
    temp = TA_F_MDS,
    vpd = VPD_F_MDS * 100,
    ppfd = SW_IN_F_MDS * kfFEC * 1e-06,
    netrad = NETRAD,
    patm = PA_F * 1000,
    snow = 0,
    rain = NA,  # P_F * 48 / (60 * 60 * 24),
    tmin = tmaxmin[names(tmaxmin)==sites[1]][[1]]$tmin, # TMIN_F_MDS,
    tmax = tmaxmin[names(tmaxmin)==sites[1]][[1]]$tmax, # TMAX_F_MDS,
    fapar = FPAR,
    co2 = CO2_F_MDS,
    ccov = 0,
    daylength = daylength
  ) |>
  list()
  return(result)
}

p_hydro_drivers_list_3<- lapply(sites,tmp_fun)
names(p_hydro_drivers_list_3) <- sites

# function 5

tmp_fun <- function(x){
  #data creation
  start_year <- valid_years %>% filter(Site==x) %>% pull(start_year)
  end_year <- valid_years %>% filter(Site==x) %>% pull(end_year)
  # tibble selection
  hhdf <-  hhdf[names(hhdf)==x][[1]]

  result <- p_hydro_drivers_list[names(p_hydro_drivers_list)==x][[1]]
  result$forcing_halfhourly <- hhdf |>
    dplyr::filter(lubridate::year(date) %in% start_year:end_year) |>
    dplyr::filter(!(lubridate::mday(date) == 29 & lubridate::month(date) == 2)) |>
    group_by(time) |>
    summarize(
      date = date,
      time = time,
      temp = TA_F_MDS,
      vpd = VPD_F_MDS * 100,
      ppfd = SW_IN_F_MDS * kfFEC * 1e-06,
      netrad = NETRAD,
      patm = PA_F * 1000,
      snow = 0,
      rain = P_F * 48 / (60 * 60 * 24),
      fapar = FPAR,
      co2 = CO2_F_MDS,
      ccov = 0
    ) |>
    list()
  return(result)
}

p_hydro_drivers_list_4<- lapply(sites,tmp_fun)
names(p_hydro_drivers_list_4) <- sites

# write all drivers to file
# apply compression to minimize space
filn <- paste0(out_path, "/",site,"_p_hydro_drivers.rda")
message(paste0("- writing to file: ", filn))
save(p_hydro_drivers,
     file = filn
     )


# JJ Note: The gpp and latenth units here are different from the demo dataset supplied with rsofun. Here the units are matched to the output units from rsofun (see conversion below)
# Write validation data

load("./data/p_model_validation.rda")
p_hydro_validation <- p_model_validation[-1]
p_hydro_validation <- p_hydro_validation[[1]]


# nested list creation
p_hydro_validationList <- list()
tmp_fun <- function(x){
  p_hydro_validationList <- append(p_hydro_validationList,p_hydro_validation)
}

p_hydro_validationList <- lapply(sites,tmp_fun)
names(p_hydro_validationList) <- sites

# function 6

tmp_fun <- function(x){
  #data creation
  start_year <- valid_years %>% filter(Site==x) %>% pull(start_year)
  end_year <- valid_years %>% filter(Site==x) %>% pull(end_year)
  # tibble selection
  ddf_24hr_mean <-  ddf_24hr_mean[names(ddf_24hr_mean)==x][[1]]
  # tibble creation
  result <- p_hydro_validationList[names(p_hydro_validationList)==x][[1]]
  result$data <- ddf_24hr_mean |>
    dplyr::filter(lubridate::year(date) %in% start_year:end_year) |>
    dplyr::filter(!(lubridate::mday(date) == 29 & lubridate::month(date) == 2)) |>
    group_by(date) |>
    summarise(
      date = date,
      time = time,
      gpp = GPP_DT_VUT_REF,
      latenth = LE_F_MDS
    ) |>
    mutate(gpp = gpp*86400/1e6*12)  |>    # convert [umol m-2 s-1] to [gC m-2 day-1]
    mutate(latenth = latenth*86400) |>    # convert [W m-2] to [J m-2 day-1]
    list()
  return(result)
}

p_hydro_validationList_1 <- lapply(sites, tmp_fun)
names(p_hydro_validationList_1) <- sites

#function 7

tmp_fun <- function(x){
  #data creation
  start_year <- valid_years %>% filter(Site==x) %>% pull(start_year)
  end_year <- valid_years %>% filter(Site==x) %>% pull(end_year)
  # tibble selection
  hhdf <-  hhdf[names(hhdf)==x][[1]]
  # tibble creation
  result <- p_hydro_validationList[names(p_hydro_validationList)==x][[1]]
  result$data_hh <-
    hhdf |>
    dplyr::filter(lubridate::year(date) %in% start_year:end_year) |>
    dplyr::filter(!(lubridate::mday(date) == 29 & lubridate::month(date) == 2)) |>
    group_by(time) |>
    summarise(
      date = date,
      time = time,
      gpp = GPP_DT_VUT_REF,
      latenth = LE_F_MDS
    ) |>
    mutate(gpp = gpp*86400/1e6*12)  |>    # convert [umol m-2 s-1] to [gC m-2 day-1]
    mutate(latenth = latenth*86400) |>    # convert [W m-2] to [J m-2 day-1]
    list()
  return(result)
}

p_hydro_validationList_2 <- lapply(sites, tmp_fun)
names(p_hydro_validationList_2) <- sites


# fnuction 8
# xxx we don't need 3-hour maxima. That's specific for phydro.
tmp_fun <- function(x){
  #data creation
  start_year <- valid_years %>% filter(Site==x) %>% pull(start_year)
  end_year <- valid_years %>% filter(Site==x) %>% pull(end_year)
  # tibble selection
  ddf_3hr_maxima <-  ddf_3hr_maxima[names(ddf_3hr_maxima)==x][[1]]
  # tibble creation
  result <- p_hydro_validationList[names(p_hydro_validationList)==x][[1]]
  result$data_3hr_mean <-
    ddf_3hr_maxima |>
    dplyr::filter(lubridate::year(date) %in% start_year:end_year) |>
    dplyr::filter(!(lubridate::mday(date) == 29 & lubridate::month(date) == 2)) |>
    group_by(time) |>
    summarise(
      date = date,
      time = time,
      gpp = GPP_DT_VUT_REF,
      latenth = LE_F_MDS
    ) |>
    mutate(gpp = gpp*86400/1e6*12)  |>    # convert [umol m-2 s-1] to [gC m-2 day-1]
    mutate(latenth = latenth*86400) |>    # convert [W m-2] to [J m-2 day-1]
    list()
  return(result)
}

p_hydro_validationList_3 <- lapply(sites, tmp_fun)
names(p_hydro_validationList_3) <- sites

# function 9

tmp_fun <- function(x){
  #data creation
  start_year <- valid_years %>% filter(Site==x) %>% pull(start_year)
  end_year <- valid_years %>% filter(Site==x) %>% pull(end_year)
  # tibble selection
  ddf_daytime_mean <-  ddf_daytime_mean[names(ddf_daytime_mean)==x][[1]]
  # tibble creation
  result <- p_hydro_validationList[names(p_hydro_validationList)==x][[1]]
  result$data_daytime_mean <-
    ddf_daytime_mean |>
    dplyr::filter(lubridate::year(date) %in% start_year:end_year) |>
    dplyr::filter(!(lubridate::mday(date) == 29 & lubridate::month(date) == 2)) |>
    group_by(time) |>
    summarise(
      date = date,
      time = time,
      gpp = GPP_DT_VUT_REF,
      latenth = LE_F_MDS
    ) |>
    mutate(gpp = gpp*86400/1e6*12)  |>    # convert [umol m-2 s-1] to [gC m-2 day-1]
    mutate(latenth = latenth*86400) |>    # convert [W m-2] to [J m-2 day-1]
    list()
  return(result)
}

p_hydro_validationList_4 <- lapply(sites, tmp_fun)
names(p_hydro_validationList_4) <- sites

filn <- paste0(out_path,"/",site,"_p_hydro_validation.rda")
message(paste0("- writing to file: ", filn))
save(p_hydro_validation,
     file = filn
     )
# P5
