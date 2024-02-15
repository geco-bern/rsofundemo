create_driver_validation <- function(site,csv_path,lsm_path,nc_path,out_path){
  
  dir.create(out_path, showWarnings = F)
  files_csv = list.files(csv_path)
  files_lsm = list.files(lsm_path)
  
  # Get filename for HH data for matching site
  file_csv = files_csv[intersect(grep(site, files_csv), 
                                 grep("HH", files_csv))]
  
  # get metadata
  # --------------------------------------------------------
  message("- reading Metadata for site")
  meta <- suppressWarnings(
    try(
      read_meta_fdk(
        site = site,
        path = lsm_path,
        meta_data = T
      )
    )
  )
  
  # get half-hourly data  --------------------------------------------------------
  # message("- convert to FLUXNET standard CSV file")
  hhdf <- suppressWarnings(
    try(
      fdk_convert_lsm(
        site = site,
        fluxnet_format = TRUE,
        path = lsm_path
      )
    )
  )
  
  if(inherits(hhdf, "try-error")){
    message("!!! conversion to FLUXNET failed  !!!")
    return(NULL)
  }
  
  
  # Add date and time columns to hhdf for easier further processing.
  # ---------------------------------------------------------
  hhdf =
    hhdf |>
    mutate(time = lubridate::as_datetime(as.character(TIMESTAMP_START), tz = "GMT", format="%Y%m%d%H%M")) |>
    mutate(date = lubridate::as_date(time))
  
  message("- Add SW_OUT=NA if not present")
  if (!("SW_OUT" %in% colnames(hhdf))) {
    hhdf$SW_OUT = NA
  }
  
  # Aggregate to daily 24-hr means  ----------------------------------------------------------
  message("- downsampling FLUXNET format - 24 hr means")
  ddf_24hr_mean <-
    try(
      hhdf |>
        group_by(date) |>
        select(-TIMESTAMP_START, -TIMESTAMP_END) |>
        summarize_all(.funs = mean)
    )
  
  valid_years = read.csv(paste0(csv_path,"/valid_years_final.csv"), header = T, fileEncoding = "UTF-16")
  
  
  # Get valid data years
  ystart = valid_years %>% filter(Site==site) %>% pull(start_year)
  yend = valid_years %>% filter(Site==site) %>% pull(end_year)
  
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
  
  message("- downsampling FLUXNET format - daytime means")
  ddf_daytime_mean <- hhdf |>
    group_by(date) |>
    do(aggregate_daily_daylength(.)) |>
    ungroup()
  
  
  tmaxmin <-
    hhdf |>
    group_by(date) |>
    summarize(
      tmax = max(TA_F_MDS),
      tmin = min(TA_F_MDS)
    )
  
  # Creating driver object  ------------------------------------------------------
  message("- compiling drivers")
  
  nc = nc_open(paste0(nc_path,"cwdx80.nc"))
  lons = ncvar_get(nc, "lon")
  lats = ncvar_get(nc, "lat")
  S80 = ncvar_get(nc, "cwdx80")
  
  site_lon = meta[[1]]$longitude
  site_lat = meta[[1]]$latitude
  
  lonid = which(lons > site_lon)[1]-1
  latid = which(lats > site_lat)[1]-1
  n = 1
  S80_slice = S80[(lonid-n):(lonid+n), (latid-n):(latid+n)]
  whc_site = mean(as.numeric(S80_slice, na.rm=T))
  whc_site_sd = sd(as.numeric(S80_slice, na.rm=T))
  
  p_hydro_drivers <- p_model_drivers
  p_hydro_drivers$sitename[[1]] <- site
  p_hydro_drivers$site_info[[1]] <-
    tibble(
      lon=meta[[1]]$longitude,
      lat=meta[[1]]$latitude,
      elv = meta[[1]]$elevation,
      #canopy_height = ifelse(is.na(meta[[1]]$canopy_height), yes = 20, meta[[1]]$canopy_height),
      #reference_height = ifelse(is.na(meta[[1]]$reference_height), yes = 22, meta[[1]]$reference_height),
      whc = whc_site
      #whc_sd = whc_site_sd,
      #IGBP_veg_short = meta[[1]]$IGBP_veg_short
    )
  kfFEC = 2.04
  
  start_year = ystart
  end_year = yend
  
  # for demo, use just a subset of years
  p_hydro_drivers$forcing <-
    ddf_24hr_mean |>
    dplyr::filter(lubridate::year(date) %in% start_year:end_year) |>
    dplyr::filter(!(lubridate::mday(date) == 29 & lubridate::month(date) == 2)) |>
    left_join(tmaxmin) |>
    group_by(date) |>
    summarize(
      date = date,
      temp = TA_F_MDS,
      vpd = VPD_F_MDS * 100,
      ppfd = SW_IN_F_MDS * kfFEC * 1e-06,
      netrad = NETRAD,
      patm = PA_F * 1000,
      snow = 0,
      rain = P_F * 48 /(60 * 60 * 24), # P_F [mm timestep-1] * 48 [timesteps day-1] / 86400 [secs day-1 ]
      tmin = tmin, # TMIN_F_MDS,
      tmax = tmax, # TMAX_F_MDS,
      fapar = FPAR,
      co2 = CO2_F_MDS,
      ccov = 0
    ) |>
    list()
  
  # write all drivers to file
  # apply compression to minimize space
  name <- p_hydro_drivers
    
  filn <- paste0(out_path,site,"_p_model_drivers.rda")
  message(paste0("- writing to file: ", filn))
  saveRDS(p_hydro_drivers,filn)
  
  # Write validation data
  message("- compiling validation")
  
  p_hydro_validation <- p_model_validation
  p_hydro_validation$sitename[[1]] = site
  
  p_hydro_validation$data <-
    ddf_24hr_mean |>
    dplyr::filter(lubridate::year(date) %in% start_year:end_year) |>
    dplyr::filter(!(lubridate::mday(date) == 29 & lubridate::month(date) == 2)) |>
    group_by(date) |>
    summarise(
      date = date,
      gpp = GPP_DT_VUT_REF,
      gpp_unc = 0
    ) |>                                          
    mutate(gpp = gpp*86400/1e6*12)  |>    # convert [umol m-2 s-1] to [gC m-2 s-1]
    list()
  
  # write all drivers to file
  # apply compression to minimize space
  filn <- paste0(out_path,site,"_p_model_validation.rda")
  message(paste0("- writing to file: ", filn))
  saveRDS(p_hydro_validation,filn)
}
