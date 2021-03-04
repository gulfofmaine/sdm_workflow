library(gmRi)
library(tidyverse)
library(raster)


# set data paths
trawlFolder <- shared.path("unix", "RES_Data", "NMFS_trawl")

climFolder <- shared.path("unix", "RES_Data", "SODA")

trawlData <- "NEFSC_BTS_02182021.RData"

climData <- "SODA_Temp_Red_bottomLayer.nc"

# load data

trawlPath <- paste0(trawlFolder, trawlData)

climPath <- paste0(climFolder, climData)

load(trawlPath)

climData <- raster::stack(climPath)

# Extract time series from lat lons
getTrawlPoints <- function(climData, survey){
  trawlLocs <- survey$survdat %>% dplyr::select(LAT, LON) %>% unique() %>% tibble()
  
  xy <- cbind(trawlLocs$LON, trawlLocs$LAT)
  
  climTimeSeries <- raster::extract(climData, xy, df=TRUE, cellnumbers=TRUE)
  
  locTSdf <- climTimeSeries %>% 
    pivot_longer(cols = c(-ID, -cells), names_to="Date", values_to="temp") %>% 
    mutate(yr = as.numeric(str_sub(Date, 2, 5)),
           mon = as.numeric(str_sub(Date, 7,8)))
  
  return(locTSdf)
}


test <- ncdf4::nc_open("/Users/mdzaugis/Box/RES_Data/CMIP6/BottomSal/StGrid/stGrid_so_FIO-ESM-2-0_r1i1p1f1_historical.nc")

times <- ncdf4::ncvar_get(test, "time")
