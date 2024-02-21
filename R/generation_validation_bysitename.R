create_validation <- function(site,csv_path,lsm_path){

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


  # Get valid data years
  valid_years = read.csv(paste0(csv_path,"/valid_years_final.csv"), header = T, fileEncoding = "UTF-16")

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

  # Write validation data
  message("- compiling validation")

  p_hydro_validation <- p_model_validation
  p_hydro_validation$sitename[[1]] = site

  p_hydro_validation$data <-
    ddf_24hr_mean |>
    dplyr::filter(lubridate::year(date) %in% ystart:yend) |>
    dplyr::filter(!(lubridate::mday(date) == 29 & lubridate::month(date) == 2)) |>
    group_by(date) |>
    summarise(
      date = date,
      gpp = GPP_DT_VUT_REF,
      gpp_unc = GPP_DT_VUT_SE
    ) |>
    mutate(gpp = gpp*86400/1e6*12)  |> # convert [umol m-2 s-1] to [gC m-2 s-1]
    mutate(gpp_unc = gpp_unc*86400/1e6*12)  |>  # convert [umol m-2 s-1] to [gC m-2 s-1]
    list()
 return(p_hydro_validation)
}
