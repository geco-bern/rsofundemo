library(dplyr)
library(tidyr)
library(rsofun)
library(here)
library(ggplot2)
library(FluxDataKit)

# function used to create dataframes
get_annual_aet_pet <- function(df){
  df |>
    mutate(year = lubridate::year(date)) |>
    group_by(year) |>
    summarise(aet = sum(aet),
              pet = sum(pet)) |>
    ungroup() |>
    summarise(aet = mean(aet),
              pet = mean(pet))
}

get_annual_prec_cond <- function(df){
  df |>
    mutate(year = lubridate::year(date)) |>
    group_by(year) |>
    summarise(prec_cond = sum(prec_cond)) |>
    ungroup() |>
    summarise(prec_cond = mean(prec_cond))
}

source(here("analysis","monthly2daily_weather.R"))

# paramter for p model
params_modl <- list(
  kphio              = 0.04998,
  kphio_par_a        = 0.0,
  kphio_par_b        = 1.0,
  soilm_thetastar    = 0.6 * 240,
  soilm_betao        = 0.0,
  beta_unitcostratio = 146.0,
  rd_to_vcmax        = 0.014,
  tau_acclim         = 30.0,
  kc_jmax            = 0.41
)

test_driver <- readRDS(here("data-raw","rsofun_driver_data_v3.rds"))

# filter by land use (remove crop and wetland)
keep <- fdk_site_info|>
  filter(igbp_land_use != "CRO" & igbp_land_use != "WET")

test_driver <- test_driver[which(test_driver$sitename %in% keep$sitename),]

# filter by good quality le_corr
keep <- fdk_site_fullyearsequence |>
  filter(drop_lecorr != TRUE)

test_driver <- test_driver[which(test_driver$sitename %in% keep$sitename),]


# get rain per month

rain = cost_whc_driver |>
  unnest(forcing) |>
  select(sitename,date,rain) 

rain$month = substr(rain$date,1,7)

monthly_rain = rain |>
  group_by(sitename,month) |>
  summarise(rain_month = sum(rain),
            wet_days = sum(rain != 0))

# since every catchment has exactly N number of years, we can call the function
# monthly2daily_weather.R without problem

# for now, the same value as the true value 

# before run the monthly2daily_weather, I multiply the monthly rain by 10^5 otherwise get 0 in everything


monthly_rain$rain_month = monthly_rain$rain_month*  60 * 60 * 24



for(i in 1:1246){
  
  # get index in df
  start_row = 1 + 12*(i-1)
  end_row =  12*i
  
  # get info to change the data
  sitename = unique(monthly_rain[start_row:end_row,]$sitename)
  year = unique(substr(monthly_rain[start_row:end_row,]$month,1,4))
  
  # get results
  mval_prec = monthly_rain[start_row:end_row,]$rain_month
  mval_wet = monthly_rain[start_row:end_row,]$wet_days
  prdaily_random = matrix(runif(730),ncol=2)
  result = monthly2daily_weather(mval_prec, mval_wet, prdaily_random)
  result = result / (60 * 60 * 24)
  
  # insert results in driver data
  fitler_1 = test_driver[test_driver$sitename == sitename,]$forcing[[1]]
  fitler_1[substr(fitler_1$date,1,4) == year,]$rain = result
}

# we can use the p model as before

test_output <- rsofun::runread_pmodel_f(
  test_driver,
  par = params_modl
)

test_output <- test_output |>
  mutate(len = purrr::map_int(data, ~nrow(.))) |>
  filter(len != 1) |>
  select(-len)

#df creation

test_adf <- test_output |>
  mutate(test_adf = purrr::map(data, ~get_annual_aet_pet(.))) |>
  unnest(test_adf) |>
  select(sitename, aet,pet)

test_adf <- test_driver |>
  unnest(forcing) |>
  left_join(
    test_output |>
      unnest(data) |>
      select(-snow, -netrad, -fapar, gpp_pmodel = gpp),
    by = c("sitename", "date")
  ) |>
  mutate(prec = (rain + snow) * 60 * 60 * 24) |>
  mutate(prec_cond = prec + cond) |>
  group_by(sitename) |>
  nest() |>
  mutate(df = purrr::map(data, ~get_annual_prec_cond(.))) |>
  unnest(df) |>
  select(sitename, prec_cond) |>
  right_join(
    test_adf,
    by = "sitename"
  )

# BUdyko relationship

ggplot(test_adf) +
  geom_point(aes(x = pet/prec_cond, y = aet/prec_cond)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dotted") +
  geom_hline(yintercept = 1, linetype = "dotted") +
  ylim(0, NA) +
  xlim(0, NA) + 
  xlab("PET / prec ") +
  ylab("AET / prec ")
