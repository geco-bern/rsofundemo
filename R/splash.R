# R version 3.2.3 (2015-12-10) -- "Wooden Christmas-Tree"
#
# splash.R
#
# VERSION: 1.0
# LAST UPDATED: 2016-02-19
#
# ~~~~~~~~
# license:
# ~~~~~~~~
# Copyright (C) 2016 Prentice Lab
#
# This file is part of the SPLASH model.
#
# SPLASH is free software: you can redistribute it and/or modify it under
# the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 2.1 of the License, or
# (at your option) any later version.
#
# SPLASH is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with SPLASH.  If not, see <http://www.gnu.org/licenses/>.
#
# ~~~~~~~~~
# citation:
# ~~~~~~~~~
# T. W. Davis, I. C. Prentice, B. D. Stocker, R. J. Whitley, H. Wang, B. J.
# Evans, A. V. Gallego-Sala, M. T. Sykes, and W. Cramer, Simple process-led
# algorithms for simulating habitats (SPLASH): Robust indices of radiation
# evapo-transpiration and plant-available moisture, Geoscientific Model
# Development, 2016 (in progress)
#
# ~~~~~~~~~~~~
# description:
# ~~~~~~~~~~~~
# This script contains functions for running SPLASH for point-based data, i.e.:
#   spin_up(list mdat, list dtot)
#   quick_run(double lat, double elv, double n, double y, double wn, double sf,
#             double tc, double pn)
#   run_one_day(double lat, double elv, double n, double y, double wn,
#               double sf, double tc, double pn)
#
# ~~~~~~~~~~
# changelog:
# ~~~~~~~~~~
# - added fix to daily soil moisture when n and ny are 365 [15.01.27]
#
#### IMPORT SOURCES ##########################################################
source(here("R","const.R"))
source(here("R","evap.R"))


#### DEFINE FUNCTIONS ########################################################
# ************************************************************************
# Name:     spin_up
# Inputs:   - list, meteorological data (mdat)
#               $num_lines ..... double, length of meteorol. variable lists
#               $lat_deg ....... double latitude (degrees)
#               $elv_m ......... double, elevation (m)
#               $year .......... double, year
#               $sf ............ list, fraction of sunshine hours
#               $tair .......... list, mean daily air temperature (deg. C)
#               $pn ............ list, precipitation (mm/d)
#           - list, daily totals (dtot)
#               $wm ............ list, daily soil moisture (mm)
# Returns:  list, daily totals
# Features: Updates the soil moisture in daily totals until equilibrium
# Depends:  quick_run
# ************************************************************************
spin_up <- function(mdat, dtot) {
    # Run one year:
    for (i in seq(from=1, to=mdat$num_lines, by=1)) {
        if (i == 1) {
            wn <- dtot$wn[mdat$num_lines]
        } else {
            wn <- dtot$wn[i - 1]
        }
        my_vals <- quick_run(mdat$lat_deg, mdat$elv_m, i, mdat$year, wn,
                             mdat$sf[i], mdat$tair[i], mdat$pn[i])
        dtot$wn[i] <- my_vals$sm
    }

    # Calculate the change:
    start_sm <- dtot$wn[1]
    end_vals <- quick_run(mdat$lat_deg, mdat$elv_m, 1, mdat$year,
                          dtot$wn[mdat$num_lines], mdat$sf[1], mdat$tair[1],
                          mdat$pn[1])
    diff_sm <- abs(end_vals$sm - start_sm)

    # Equilibrate:
    spin_count <- 1
    while (diff_sm > 1.0) {
        for (i in seq(from=1, to=mdat$num_lines, by=1)) {
            if (i == 1) {
                wn <- dtot$wn[mdat$num_lines]
            } else {
                wn <- dtot$wn[i - 1]
            }
            my_vals <- quick_run(mdat$lat_deg, mdat$elv_m, i, mdat$year, wn,
                                 mdat$sf[i], mdat$tair[i], mdat$pn[i])
            dtot$wn[i] <- my_vals$sm
        }

        # Calculate the change:
        start_sm <- dtot$wn[1]
        end_vals <- quick_run(mdat$lat_deg, mdat$elv_m, 1, mdat$year,
                              dtot$wn[mdat$num_lines], mdat$sf[1],
                              mdat$tair[1], mdat$pn[1])
        diff_sm <- abs(end_vals$sm - start_sm)
        spin_count <- spin_count + 1
    }
    cat(paste("Spun", spin_count, "years\n"))
    return(dtot)
}


# ************************************************************************
# Name:     quick_run
# Inputs:   - double, latitude, deg (lat)
#           - double, elevation, m (elv)
#           - double, day of year (n)
#           - double, year (y)
#           - double, daily soil moisture content, mm (wn)
#           - double, daily fraction of bright sunshine (sf)
#           - double, daily air temperature, deg C (tc)
#           - double, daily precipitation, mm (pn)
# Returns:  list
#             $sm - soil moisture, mm
#             $ro - runoff, mm
# Features: Returns daily soil moisture and runoff
# Depends:  evap
# ************************************************************************
quick_run <- function(lat, elv, n, y, wn, sf, tc, pn) {
    # Calculate evaporative supply (mm/hr)
    sw <- kCw*wn/kWm

    # Compute daily radiation and evaporations values:
    ET <- calc_daily_evap(lat, n, elv, y, sf, tc, sw)

    # Update daily soil moisture:
    sm <- wn + pn + ET$cond_mm - ET$aet_mm

    if (sm > kWm) {
        # Bucket is full:
        # - set soil moisture to capacity
        # - add remaining water to runoff
        ro <- sm - kWm
        sm <- kWm
    } else if (sm < 0) {
        # Bucket is empty:
        # - set runoff and soil moisture equal to zero
        ro <- 0
        sm <- 0
    } else {
        ro <- 0
    }

    rval <- list()
    rval$sm <- sm
    rval$ro <- ro
    return(rval)
}


# ************************************************************************
# Name:     run_one_day
# Inputs:   - double, latitude, deg (lat)
#           - double, elevation, m (elv)
#           - double, day of year (n)
#           - double, year (y)
#           - double, daily soil moisture content, mm (wn)
#           - double, daily fraction of bright sunshine (sf)
#           - double, daily air temperature, deg C (tc)
#           - double, daily precipitation, mm (pn)
# Returns:  list
#             $ho - daily solar irradiation, J/m2
#             $hn - daily net radiation, J/m2
#             $ppfd - daily PPFD, mol/m2
#             $cond - daily condensation water, mm
#             $eet - daily equilibrium ET, mm
#             $pet - daily potential ET, mm
#             $aet - daily actual ET, mm
#             $wn - daily soil moisture, mm
#             $ro - daily runoff, mm
# Features: Runs SPLASH at a single location for one day.
# Depends:  evap
# ************************************************************************
run_one_day <- function(lat, elv, n, y, wn, sf, tc, pn) {
    # Return values
    rvals <- list()

    # Calculate evaporative supply (mm/hr)
    sw <- kCw*wn/kWm

    # Compute daily radiation and evaporations values:
    ET <- calc_daily_evap(lat, n, elv, y, sf, tc, sw)
    rvals$ho <- ET$ra_j.m2
    rvals$hn <- ET$rn_j.m2
    rvals$ppfd <- ET$ppfd_mol.m2
    rvals$cond <- ET$cond_mm
    rvals$eet <- ET$eet_mm
    rvals$pet <- ET$pet_mm
    rvals$aet <- ET$aet_mm

    # Update daily soil moisture:
    sm <- wn + pn + ET$cond_mm - ET$aet_mm

    print( paste( "in run_one_day: pn =", pn))

    if (sm > kWm) {
        # Bucket is full:
        # - set soil moisture to capacity
        # - add remaining water to runoff
        ro <- sm - kWm
        sm <- kWm
    } else if (sm < 0) {
        # Bucket is empty:
        # - reduce actual ET by discrepancy amount
        # - set runoff and soil moisture equal to zero
        rvals$aet <- rvals$aet + sm
        ro <- 0
        sm <- 0
    } else {
        ro <- 0
    }

    rvals$wn <- sm
    rvals$ro <- ro
    return(rvals)
}
