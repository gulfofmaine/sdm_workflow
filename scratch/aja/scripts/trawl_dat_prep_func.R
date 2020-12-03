trawl_dat_prep<- function(survdat_path = "~/Volumes/Shared/Research/MillsLab/NMFS Trawl Data/Survdat_Nye_allseason.RData", out_path = "./Data/") {
  ## Details
  # This function reads in an .RData file of NOAA NEFSC bottom trawl survey data and then makes some modifications to produce a data file ready for fitting with SDM, such that for every tow there is a record of every species (presence, log(biomass+1) given presence) occurrence. The function saves the tidy data set as "model_dat.rds."
  
  # Args:
  # survdat_path = Path to survey data file 
  # out_path = Path to save processsed rds data file
  
  # Returns: RDS Data file, which is also saved in folder specified by out_path
  
  ## Start function
  # Preliminaries -----------------------------------------------------------
  # Load libraries, using package_check to download if it is missing
  libraries_needed<- c("tidyverse", "gmRi")
  library_check(libraries_needed)
  
  # For debugging
  if(FALSE){
    survdat_path = paste(shared.path(os.use = os_use, group = "root", folder = "RES Data/NMFS_trawl/"), "Survdat_Nye_allseason.Rdata", sep = "")
    out_path = here::here("/scratch/aja/data/")
  }
  
  # Load in the data as survdat
  load(survdat_path)
  
  # Data cleaning -----------------------------------------------------------
  dat<- survdat %>%
    # We won't need all of the columns available in the raw survey data
    dplyr::select(ID, EST_YEAR, EST_MONTH, SEASON, STRATUM, DECDEG_BEGLAT, DECDEG_BEGLON, COMNAME, CATCHSEX, BIOMASS, AVGDEPTH, ABUNDANCE, LENGTH, NUMLEN) %>%
    # Filtering out some strata which are in Canadian waters, as well as those more recently part of the NOAA SEFSC seasonal surveys
    filter(STRATUM >= 01010 & STRATUM <= 01760) %>%
    filter(STRATUM != 1310 & STRATUM != 1320 & STRATUM != 1330 & STRATUM != 1350 &
             STRATUM != 1410 & STRATUM != 1420 & STRATUM != 1490) %>%
    filter(SEASON == "SPRING" | SEASON == "FALL") %>%
    # Some work with biomass as there are occasionally issues with biomass and abundance 
    mutate(BIOMASS = ifelse(is.na(BIOMASS) == TRUE & ABUNDANCE > 0, 0.01, BIOMASS)) %>% 
    mutate(BIOMASS = ifelse(BIOMASS == 0 & ABUNDANCE > 0, 0.01, BIOMASS)) %>% 
    mutate(ABUNDANCE = ifelse(is.na(ABUNDANCE) == TRUE & BIOMASS > 0, 1, ABUNDANCE)) %>% 
    mutate(ABUNDANCE = ifelse(ABUNDANCE == 0 & BIOMASS > 0, 1, ABUNDANCE)) %>%   
    filter(!is.na(BIOMASS),
           !is.na(ABUNDANCE)) %>%
    # Keep distinct rows, ignoring duplicated rows for catch sex of different lengths
    distinct(ID, EST_YEAR, EST_MONTH, SEASON, STRATUM, DECDEG_BEGLAT, DECDEG_BEGLON, COMNAME, CATCHSEX, .keep_all = TRUE) %>%
    # Now, aggregate to get our total abundance and biomass by tow-species
    group_by(., ID, EST_YEAR, EST_MONTH, SEASON, STRATUM, DECDEG_BEGLAT, DECDEG_BEGLON, COMNAME, AVGDEPTH) %>%
    summarize(., "SUM_BIOMASS" = sum(BIOMASS, na.rm = TRUE),
              "SUM_ABUNDANCE" = sum(ABUNDANCE, na.rm = TRUE))
  
  # Now, we have unique tow-species-biomass-abundance records for tows were a species was caught. However, we also need to absences, so we need biomass and abundance of zero for each species when they were not caught at a specific tow.
  # First, we create a null dataset for merging with all species and all trawls
  null_dat<- expand.grid("ID" = unique(dat$ID), "COMNAME" = unique(dat$COMNAME)) %>%
    mutate(., "COMNAME" = as.character(COMNAME),
           "SUM_BIOMASS" = rep(0, nrow(.)),
           "SUM_ABUNDANCE" = rep(0, nrow(.)))
  
  # Add in "tow" information based on ID...
  tow_info<- dat %>%
    ungroup() %>%
    select(., ID, EST_YEAR, EST_MONTH, SEASON, STRATUM, DECDEG_BEGLAT, DECDEG_BEGLON, AVGDEPTH) %>%
    distinct() %>% 
    filter(., ID %in% null_dat$ID)
  null_dat_join<- null_dat %>%
    left_join(., tow_info, by = c("ID" = "ID")) %>%
    select(., ID, EST_YEAR, EST_MONTH, SEASON, STRATUM, DECDEG_BEGLAT, DECDEG_BEGLON, COMNAME, AVGDEPTH, SUM_BIOMASS, SUM_ABUNDANCE) %>%
    group_by(ID, EST_YEAR, EST_MONTH, SEASON, STRATUM, DECDEG_BEGLAT, DECDEG_BEGLON, COMNAME)
  
  # Keep the rows from null_dat that are NOT in the original dataset -- these would be the rows we want to add for the imputed absences
  null_dat_add<- anti_join(null_dat_join, dat, by = c("ID" = "ID", "COMNAME" = "COMNAME"))
  
  # Full join with the presence data
  dat<- dat %>%
    full_join(., null_dat_add)
  
  # One final addition, adding in a PRESENCE column and a LOG.BIOMASS column
  dat$PRESENCE<- ifelse(dat$SUM_BIOMASS > 0, 1, dat$SUM_BIOMASS) # Create presence/absence vector based on sum biomass
  dat$LOG_BIOMASS<- log(dat$SUM_BIOMASS+1)
  
  # Return and write out processed file -----------------------------------------------------------
  return(dat)
  saveRDS(dat, file = paste(out_path, "model_dat.rds", sep = ""))
  
  #########
  ## End
  #########
}
