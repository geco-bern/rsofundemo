library(dplyr)
library(tidyr)
library(ggplot2)
library(rsofun)
library(here)
library(lubridate)

# We use `rsofun_driver_data_clean.rds`, provided on Zenodo
# ([Hufkens, 2022](https://doi.org/10.5281/zenodo.8403081)).
# Download that file and specify its local path. Note: adjust this path
# depending on where the file is located on your computer.
path_forcingfile <- "~/data/FluxDataKit/rsofun_driver_data_clean.rds"
driver <- readRDS(path_forcingfile)

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

# run the model
output <- rsofun::runread_pmodel_f(
  driver,
  par = params_modl
)

# some simulations failed due to missing values in the forcing. If that happens,
# the number of rows in the output data is 1.
output <- output |>
  mutate(len = purrr::map_int(data, ~nrow(.))) |>
  filter(len != 1) |>
  select(-len)

# aggregate outputs to annual totals, then average across years (aet and pet are in mm d-1, )
get_annual_totals <- function(df){
  df |>
    mutate(year = lubridate::year(date)) |>
    group_by(year) |>
    summarise(aet = sum(aet),
              pet = sum(pet)) |>
    ungroup() |>
    summarise(aet = mean(aet),
              pet = mean(pet))
}

output <- output |>
  mutate(adf = purrr::map(data, ~get_annual_totals(.))) |>
  unnest(adf)

# Bug: AET > PET for many sites - What is happening???
output |>
  ggplot(aes(pet, aet)) +
  geom_point() +
  geom_abline(slope = 1, linetype = "dotted")

# Look at time series of an individual site for which AET > PET
tmp <- output |>
  filter(sitename == "CZ-BK1") |>
  pull(data)

ggplot(aes(x = date),
       data = tmp[[1]] |>
         slice(1:365)) +
  geom_line(aes(y = aet), color = "royalblue") +
  geom_line(aes(y = pet), color = "tomato")
