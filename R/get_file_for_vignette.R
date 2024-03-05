library(dplyr)
library(tidyr)
library(rsofun)
library(here)
library(ncdf4)

# funuction used to create dataframes
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

# paramter for p model
params_modl <- list(
  kphio              = 0.04998,    # setup ORG in Stocker et al. 2020 GMD
  kphio_par_a        = 0.0,        # set to zero to disable temperature-dependence of kphio
  kphio_par_b        = 1.0,
  soilm_thetastar    = 0.6 * 240,  # to recover old setup with soil moisture stress
  soilm_betao        = 0.0,
  beta_unitcostratio = 146.0,
  rd_to_vcmax        = 0.014,      # value from Atkin et al. 2015 for C3 herbaceous
  tau_acclim         = 30.0,
  kc_jmax            = 0.41
)


# run p model
old_driver <- readRDS(here("data","rsofun_driver_data_clean.rds"))

old_output <- rsofun::runread_pmodel_f(
  old_driver,
  par = params_modl
)

# some simulations failed due to missing values in the forcing. If that happens,
# the number of rows in the output data is 1.
old_output <- old_output |>
  mutate(len = purrr::map_int(data, ~nrow(.))) |>
  filter(len != 1) |>
  select(-len)

# create dataframe
old_adf <- old_output |>
  mutate(old_adf = purrr::map(data, ~get_annual_aet_pet(.))) |>
  unnest(old_adf) |>
  select(sitename, aet, pet)


old_adf <- old_driver |>
  unnest(forcing) |>
  left_join(
    old_output |>
      unnest(data) |>
      select(-snow, -netrad, -fapar, gpp_pmodel = gpp),
    by = c("sitename", "date")
  ) |>
  mutate(prec = (rain + snow) * 60 * 60 * 24) |>
  mutate(prec_cond = prec + cond) |>
  group_by(sitename) |>
  nest() |>
  mutate(old_adf = purrr::map(data, ~get_annual_prec_cond(.))) |>
  unnest(old_adf) |>
  select(sitename, prec_cond) |>
  right_join(
    old_adf,
    by = "sitename"
  )

whc <-  old_output |>
  unnest(site_info) |>
  select(whc)

old_adf$whc <- whc[[1]]

# updating whc

nc_file <- nc_open(here("data","cwdx80.nc"))

whc = ncvar_get(nc_file, "cwdx80")
lons = ncvar_get(nc_file, "lon")
lats = ncvar_get(nc_file, "lat")

geo <- old_driver |>
  unnest(site_info) |>
  select(lon  , lat)

geo$sitename <- old_driver$sitename 

n <- 1 # parameter to select size of slice to average

new_whc <- lapply(geo$sitename, function(x){
  geo <- geo[geo$sitename == x,]
  lonid <- which(lons > geo$lon)[1]
  latid <- which(lats > geo$lat)[1]
  zroot_region <- whc[(lonid-n):(lonid+n), (latid-n):(latid+n)]
  whc_site <- mean(as.numeric(zroot_region, na.rm=T))
  return(whc_site)
  })

new_whc <- unlist(new_whc)

new_driver <- old_driver

for(i in 1:241){
  new_driver$site_info[i][[1]][4] <- new_whc[i]
}
rm(whc) # for memnory
```

```{r warning=FALSE}
# run p model
new_output <- rsofun::runread_pmodel_f(
  new_driver,
  par = params_modl
)

# some simulations failed due to missing values in the forcing. If that happens,
# the number of rows in the output data is 1.
new_output <- new_output |>
  mutate(len = purrr::map_int(data, ~nrow(.))) |>
  filter(len != 1) |>
  select(-len)

# create dataframe
new_adf <- new_output |>
  mutate(new_adf = purrr::map(data, ~get_annual_aet_pet(.))) |>
  unnest(new_adf) |>
  select(sitename, aet, pet)


new_adf <- new_driver |>
  unnest(forcing) |>
  left_join(
    new_output |>
      unnest(data) |>
      select(-snow, -netrad, -fapar, gpp_pmodel = gpp),
    by = c("sitename", "date")
  ) |>
  mutate(prec = (rain + snow) * 60 * 60 * 24) |>
  mutate(prec_cond = prec + cond) |>
  group_by(sitename) |>
  nest() |>
  mutate(new_adf = purrr::map(data, ~get_annual_prec_cond(.))) |>
  unnest(new_adf) |>
  select(sitename, prec_cond) |>
  right_join(
    new_adf,
    by = "sitename"
  )

whc <-  new_output |>
  unnest(site_info) |>
  select(whc)

new_adf$whc <- whc[[1]]

saveRDS(old_adf,here("data","budyko_old_whc.rds"))
saveRDS(new_adf,here("data","budyko_new_whc.rds"))
saveRDS(geo,here("data","site_coordinates.rds"))
