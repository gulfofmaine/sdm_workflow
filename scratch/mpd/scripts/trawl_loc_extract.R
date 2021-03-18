library(gmRi)
library(tidyverse)
library(raster)


# set data paths
trawlFolder <- shared.path("unix", "RES_Data", "NMFS_trawl")

climFolder <- shared.path("unix", "RES_Data", "CMIP6")

trawlData <- "NEFSC_BTS_02182021.RData"

climData <- "BiasCorrected/surf_temp_OISST_bias_corrected_mean.grd"

# load data

trawlPath <- paste0(trawlFolder, trawlData)

climPath <- paste0(climFolder, climData)

load(trawlPath)

climData <- raster::stack(climPath)

# Extract time series from lat lons
getTrawlPoints <- function(climData, survey, time){
  
  time <- rlang::ensym(time)
  trawlLocs <- survey$survdat %>% dplyr::select(LAT, LON) %>% unique() %>% tibble()
  
  xy <- cbind(gmRi::make360(trawlLocs$LON), trawlLocs$LAT)
  
  climTimeSeries <- raster::extract(climData, xy, df=TRUE) %>% 
    bind_cols("LON" = trawlLocs$LON, "LAT" = trawlLocs$LAT)
  
  # if you are looking to speed up this process - the pivot_longer is the step that takes a long time
  
  locTSdf <- climTimeSeries %>% 
    pivot_longer(cols = c(-ID, -LON, -LAT), names_to="Date", values_to="temp") %>% 
    mutate(yr = as.numeric(str_sub(Date, 2, 5)),
           mon = as.numeric(str_sub(Date, 7,8)),
           qt = case_when(mon %in% c(1,2,3) ~ 'winter',
                          mon %in% c(4,5,6) ~ 'spring',
                          mon %in% c(7,8,9) ~ 'summer',
                          mon %in% c(10,11,12) ~ 'fall')) %>%
    group_by(!!time, yr, ID, LAT, LON) %>% 
    summarise(temp = mean(temp), 
              Date = paste(!!time, yr, sep = "-"),
              .groups = "drop") %>% dplyr::select(-!!time, -yr) %>% 
    pivot_wider(names_from = Date, values_from = temp)
    
  return(locTSdf)
}


trawlLocsTemp <- getTrawlPoints(climData, survey, time="qt")

head(trawlLocsTemp)
