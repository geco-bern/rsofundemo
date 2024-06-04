library(tidyverse)
library(readr)
library(here)
library(ggplot2)


# I'll call sites the sites from rsofun_driver_data_v3 and location the sites from caravan

# I check if the total precipitation in the closest location in caravan_US is similar to our total precipitation

# insert your path to file here
df_rsofun <- readRDS(here("data-raw","rsofun_driver_data_v3.rds"))
meta_data_1 <- read_csv("C:/Users/Lenovo/Downloads/Nuova cartella/caravan/attributes/camels/attributes_other_camels.csv")
meta_data_2 <- read_csv("C:/Users/Lenovo/Downloads/Nuova cartella/caravan/attributes/hysets/attributes_other_hysets.csv")

meta_data <-  meta_data_1 #rbind(meta_data_1,meta_data_2) to include also hysets dataset
rm(meta_data_1,meta_data_2)

# since are alphabetically ordered, the US sites are between 198 and 335

df_rsofun <- df_rsofun[198:335,]

coordinates <- df_rsofun|>
  unnest(site_info) |>
  select(sitename,lon,lat,whc)

# for each site in our driver data, I search for the closest site in caravan

df <- NULL

for(i in 1:dim(coordinates)[1]){
  filter <- sqrt((coordinates$lon[i] - meta_data$gauge_lon)^2 + 
                 (coordinates$lat[i] - meta_data$gauge_lat)^2)  

  
  tmp <- data.frame(sitename = coordinates$sitename[i],
                    tower_lon = coordinates$lon[i],
                    tower_lat = coordinates$lat[i],
                    location = meta_data$gauge_id[which(filter %in% min(filter))],
                    lon = meta_data$gauge_lon[which(filter %in% min(filter))],
                    lat = meta_data$gauge_lat[which(filter %in% min(filter))],
                    distance = min(filter))
  
  df = rbind(df,tmp)
}

# filter out far sites (more than 2 degree)

df <- df[df$distance < 2,]

# as check, I calculate the total precipitation in that site and compare to to our result

driver_data_precipitation =  df_rsofun|>
  unnest(forcing) |>
  select(sitename,date,rain) |>
  group_by(sitename) |>
  summarise(site_rain = sum(rain)*24*60*60)

driver_data_precipitation <- driver_data_precipitation[which(driver_data_precipitation$sitename %in% df$sitename),]


# load caravan data and check

file_to_open <- unique(df$location)

path_to_file <- "C:/Users/Lenovo/Downloads/Nuova cartella/caravan/timeseries/csv/camels/"


# to have the same dates

date = df_rsofun|>
  unnest(forcing) |>
  select(sitename,date)

df_rain = NULL

for(i in 1:length(file_to_open)){
  
  file = read.csv(paste0(path_to_file,file_to_open[i],".csv"))
  # each location is associated to multiple sites, so we need another cycle
  associated_sites <- df[df$location == file_to_open[i],]$sitename
  
  for(j in 1:length(associated_sites)){
    filter_date <- date[date$sitename == associated_sites[j],]$date
    
    rain <- sum(file[which(file$date %in% filter_date),38])
    
    tmp <- data.frame(location = file_to_open[i],
                      sitename = associated_sites[j],
                      location_rain = rain)
    
    df_rain = rbind(df_rain,tmp)
  }
}


df_rain <- left_join(df_rain,driver_data_precipitation, by="sitename")

ggplot(df_rain, aes(x = site_rain, y = location_rain)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0, linetype = "dotted") +
  ggtitle("Total precipitaion in rsofun driver and caravan database")
