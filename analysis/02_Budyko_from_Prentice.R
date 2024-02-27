library(here)
source(here("R","const.R"))
source(here("R","data.R"))
source(here("R","evap.R"))
source(here("R","splash.R"))
driver <- readRDS(here("data","rsofun_driver_data_clean.rds"))

create_data <- function(df_driver) {
    my_data <- list()
    my_data$file_name <- df_driver$sitename
    tmp <- df_driver[[4]][[1]]
    my_data$sf <- 1- tmp$ccov 
    my_data$tair <- tmp$temp
    my_data$pn <- tmp$rain * (60*24*24)
    my_data$num_lines <- nrow
    my_data$year <- as.numeric(substr(tmp$date[1],1,4))
    tmp <- df_driver[[3]][[1]]
    my_data$lat_deg <- tmp$lat
    my_data$elv_m <- tmp$elv

    return (my_data)}

# run singe site
sub_df <- driver[1,]

nrow = dim(sub_df[[4]][[1]])[1]
daily_totals <- matrix(data=rep(0,  nrow*9), nrow=nrow, ncol=9)
daily_totals <- as.data.frame(daily_totals)
names(daily_totals) <- c("ho",   # daily solar irradiation, J/m2
                         "hn",   # daily net radiation, J/m2
                         "qn",   # daily PPFD, mol/m2
                         "cn",   # daily condensation, mm
                         "wn",   # daily soil moisture, mm
                         "ro",   # daily runoff, mm
                         "eq_n", # daily equilibrium ET, mm
                         "ep_n", # daily potential ET, mm
                         "ea_n") # daily actual ET, mm


y <- 2013 # to fix, need mutliple year as input
ny <- 365

my_data <- create_data(sub_df)
# Spin up the soil moisture content
daily_totals <- spin_up(my_data, daily_totals)

# Run one day:
daily_vals <- run_one_day(my_data$lat_deg,
                            my_data$elv_m,
                            172,
                            my_data$year,
                            145.0,0.5,17.3,10.0)

for (i in rep(1:365,times=dim(prova[[4]][[1]])[1]/365)) {
    n <- i
    idx <- (n - 1)
    if (idx < 1) {
        idx <- ny}
    daily_vals <- run_one_day(my_data$lat_deg,
                              my_data$elv_m,
                              n,
                              my_data$year,
                              daily_totals$wn[idx],
                              my_data$sf[n],
                              my_data$tair[n],
                              my_data$pn[n])

    # Update daily values:
    daily_totals$wn[n] <- daily_vals$wn
    daily_totals$ro[n] <- daily_vals$ro

    # Save daily results:
    daily_totals$ho[n] <- daily_vals$ho
    daily_totals$hn[n] <- daily_vals$hn
    daily_totals$qn[n] <- daily_vals$ppfd
    daily_totals$cn[n] <- daily_vals$cond
    daily_totals$eq_n[n] <- daily_vals$eet
    daily_totals$ep_n[n] <- daily_vals$pet
    daily_totals$ea_n[n] <- daily_vals$aet
}

# should be less than one
mean(daily_totals$ea_n)/mean(my_data$pn)


# run for multiple sites

AET <- c()
PET <- c()
sites <- c()

# run for different site
for(i in 1:100){
  sub_df <- driver[i,]
  
  nrow = dim(sub_df[[4]][[1]])[1]
  daily_totals <- matrix(data=rep(0,  nrow*9), nrow=nrow, ncol=9)
  daily_totals <- as.data.frame(daily_totals)
  names(daily_totals) <- c("ho",   # daily solar irradiation, J/m2
                           "hn",   # daily net radiation, J/m2
                           "qn",   # daily PPFD, mol/m2
                           "cn",   # daily condensation, mm
                           "wn",   # daily soil moisture, mm
                           "ro",   # daily runoff, mm
                           "eq_n", # daily equilibrium ET, mm
                           "ep_n", # daily potential ET, mm
                           "ea_n") # daily actual ET, mm
  
  
  y <- 2013 # to fix, need mutliple year as input
  ny <- 365
  
  my_data <- create_data(sub_df)
  # Spin up the soil moisture content
  daily_totals <- spin_up(my_data, daily_totals)
  
  # Run one day:
  daily_vals <- run_one_day(my_data$lat_deg,
                            my_data$elv_m,
                            172,
                            my_data$year,
                            145.0,0.5,17.3,10.0)
  
  for (i in rep(1:365,times=dim(prova[[4]][[1]])[1]/365)) {
    n <- i
    idx <- (n - 1)
    if (idx < 1) {
      idx <- ny}
    daily_vals <- run_one_day(my_data$lat_deg,
                              my_data$elv_m,
                              n,
                              my_data$year,
                              daily_totals$wn[idx],
                              my_data$sf[n],
                              my_data$tair[n],
                              my_data$pn[n])
    
    # Update daily values:
    daily_totals$wn[n] <- daily_vals$wn
    daily_totals$ro[n] <- daily_vals$ro
    
    # Save daily results:
    daily_totals$ho[n] <- daily_vals$ho
    daily_totals$hn[n] <- daily_vals$hn
    daily_totals$qn[n] <- daily_vals$ppfd
    daily_totals$cn[n] <- daily_vals$cond
    daily_totals$eq_n[n] <- daily_vals$eet
    daily_totals$ep_n[n] <- daily_vals$pet
    daily_totals$ea_n[n] <- daily_vals$aet
    
  }
  
  # AET should be less than one
  
  AET <- append(AET,mean(daily_totals$ea_n)/mean(my_data$pn))
  PET <- append(AET,mean(daily_totals$ep_n)/mean(my_data$pn))
}
