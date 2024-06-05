library(dplyr)
library(tidyr)
library(rsofun)
library(here)
library(ggplot2)
library(readr)

# In this function I calculated the precipitation as sum of the entire period
source(here("caravan_database","function_check.R"))


driver = readRDS(here("data-raw","rsofun_driver_data_v3.rds"))

# these files will be used to calculate the "observed AET"
df <- total_precipitation_comparison("US",0.05,"C:/Users/Lenovo/Desktop/berna/rsofun-master/data-raw/rsofun_driver_data_v3.rds",
                                     "C:/Users/Lenovo/Downloads/Nuova cartella/caravan/attributes/hysets/attributes_other_hysets.csv",
                                     "C:/Users/Lenovo/Downloads/Nuova cartella/caravan/timeseries/csv/hysets/")

# check precipitation and location 

ggplot(prova, aes(x = site_rain, y = location_rain,color=distance)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0, linetype = "dotted") +
  ggtitle("Total precipitaion in rsofun driver and caravan database")

world <- map_data("world")

US <- world[world$long < -70 & world$long > -130,]
US <- US[US$lat < 50 & US$lat > 25,]

ggplot(prova) +
  geom_map(data = US, map = US,aes(long, lat, map_id = region),
           color = "black", fill = "lightgray", size = 0.1) +
  geom_point(aes(x= tower_lon, y = tower_lat, color= "red")) +
  geom_point(aes(x= lon, y = lat, color= "blue"))



path_to_caravan_csv = "C:/Users/Lenovo/Downloads/Nuova cartella/caravan/timeseries/csv/hysets/"

obs_pred_AET = NULL


for(i in unique(df$location)){
  file = read.csv(paste0(path_to_caravan_csv,i,".csv"))
  message("file opened")
  
  AET = file |>
    mutate(date = substr(date,1,4)) |>
    group_by(date) |>
    summarise(AET = sum(total_precipitation_sum) - sum(streamflow, na.rm= TRUE))
  
  
  # each location is associated to multiple sites, so we need another cycle
  associated_sites <- df[df$location == i,]$sitename
  
  for(j in associated_sites){
    sub_driver = driver[driver$sitename == j,]
    date = sub_driver  |>
      unnest(forcing) |>
      select(date)|>
      mutate(date =  substr(date,1,4))
    
    sub_AET = AET[which(AET$date %in% unique(date$date)),]
    
    tmp = data.frame(sitename = j,
                     location = i,
                     observed_AET = mean(sub_AET$AET))
    
    
    obs_pred_AET = rbind(obs_pred_AET,tmp)
  }
}

# Now I'll import the output from P model

get_annual_aet <- function(df){
  df |>
    mutate(year = lubridate::year(date)) |>
    group_by(year) |>
    summarise(aet = sum(aet)) |>
    ungroup() |>
    summarise(aet = mean(aet))
}

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

# run P model variable WHC

output <- rsofun::runread_pmodel_f(
  driver,
  par = params_modl
)

output <- output |>
  mutate(len = purrr::map_int(data, ~nrow(.))) |>
  filter(len != 1) |>
  select(-len)

#df creation

variable_whc <- output |>
  mutate(df = purrr::map(data, ~get_annual_aet(.))) |>
  unnest(df) |>
  select(sitename, aet)  |>
  mutate(variable_aet = aet) |>
  select( - aet)

variable_whc =  variable_whc[which(variable_whc$sitename %in% obs_pred_AET$sitename),]
obs_pred_AET =  obs_pred_AET[which(obs_pred_AET$sitename %in% variable_whc$sitename),]

# change to previous WHC

whc_2m = read.csv(here("data-raw","whc_2m.csv"),sep=" ")

for(i in 1:dim(driver)[1]){
  driver$site_info[i][[1]][4] <- whc_2m$WHC[i]
}

# run P model costant WHC

output <- rsofun::runread_pmodel_f(
  driver,
  par = params_modl
)

output <- output |>
  mutate(len = purrr::map_int(data, ~nrow(.))) |>
  filter(len != 1) |>
  select(-len)

#df creation

costant_WHC <- output |>
  mutate(df = purrr::map(data, ~get_annual_aet(.))) |>
  unnest(df) |>
  select(sitename, aet)  |>
  mutate(costant_aet = aet) |>
  select( - aet)

costant_WHC =  costant_WHC[which(costant_WHC$sitename %in% obs_pred_AET$sitename),]

# join the results of P model

obs_pred_AET = left_join(obs_pred_AET,variable_whc, by="sitename")
obs_pred_AET = left_join(obs_pred_AET,costant_WHC, by="sitename")


ggplot(obs_pred_AET)+
  geom_point(aes(x = observed_AET, y = variable_aet, color= "variable")) +
  geom_point(aes(x = observed_AET, y = costant_aet, color= "costant")) +
  geom_abline(slope = 1, intercept = 0, linetype = "dotted") +
  xlab("AET from caravan") +
  ylab("AET from P model")

