total_precipitation_comparison <- function(state,distance_threshold,path_to_driver,path_to_metadata,path_to_caravan_csv){
  
  
  df_rsofun <- readRDS(path_to_driver)
  meta_data <- read_csv(path_to_metadata)
  message("-- metadata loaded")
  
  df_rsofun <- df_rsofun[grep(state,df_rsofun$sitename),]
  
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
  
  df <- df[df$distance < distance_threshold,]
  
  # as check, I calculate the total precipitation in that site and compare to to our result
  
  driver_data_precipitation =  df_rsofun|>
    unnest(forcing) |>
    select(sitename,date,rain) |>
    group_by(sitename) |>
    summarise(site_rain = sum(rain)*24*60*60)
  
  driver_data_precipitation <- driver_data_precipitation[which(driver_data_precipitation$sitename %in% df$sitename),]
  
  
  # load caravan data and check
  
  file_to_open <- unique(df$location)
 
  
  # to have the same dates
  
  date = df_rsofun|>
    unnest(forcing) |>
    select(sitename,date)
  
  df_rain = NULL
  
  message(paste0("-- starts to open caravan files, location = ",length(file_to_open)," sites = ",dim(df)[1]))
  
  for(i in 1:length(file_to_open)){
    
    file = read.csv(paste0(path_to_caravan_csv,file_to_open[i],".csv"))
    
    message("-- open file ",i," out of ",length(file_to_open))
    
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
  
  df_rain <- left_join(df_rain,df|>select(- location), by="sitename")
}



