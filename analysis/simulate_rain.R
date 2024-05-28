library(dplyr)
library(tidyr)
library(rsofun)
library(here)
library(ncdf4)
library(cwd)
library(FluxDataKit)

source(here("data-raw","monthly2daily_weather.R"))

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

cost_whc_driver <- test_driver <- readRDS(here("data-raw","filtered_driver.rds"))

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


