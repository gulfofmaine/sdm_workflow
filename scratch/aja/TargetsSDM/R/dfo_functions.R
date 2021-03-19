####  Common Resources  ####


####  Functions  ####
####

## Dynamic loading functions -- each takes in a directory path and reads in files from it. 
dfo_GSINF_load<- function(dfo_raw_dir){
  # Load raw data
  rdata_files_load<- list.files(dfo_raw_dir, pattern = ".GSINF", full.names = TRUE)
  lapply(rdata_files_load, load, environment(), verbose = FALSE)
  
  # New name
  dfo_GSINF<- GSINF
  return(dfo_GSINF)
}

dfo_GSMISSIONS_load<- function(dfo_raw_dir){
  # Load raw data
  rdata_files_load<- list.files(dfo_raw_dir, pattern = ".GSMISSIONS", full.names = TRUE)
  lapply(rdata_files_load, load, environment(), verbose = FALSE)
  
  # New name
  dfo_GSMISSIONS<- GSMISSIONS
  return(dfo_GSMISSIONS)
}

dfo_GSCAT_load<- function(dfo_raw_dir){
  # Load raw data
  rdata_files_load<- list.files(dfo_raw_dir, pattern = ".GSCAT", full.names = TRUE)
  lapply(rdata_files_load, load, environment(), verbose = FALSE)
  
  # New name
  dfo_GSCAT<- GSCAT
  return(dfo_GSCAT)
}

## Cleaning and processing functions
#' @title DFO tow data
#' 
#' @description This function uses the `dfo_GSINF` and `dfo_GSMISSIONS` databases to create a dataframe that has one row for every unique tow, with column names that match those readily available from the NOAA NMFS surveys.
#' 
#' @param dfo_GSINF
#' @param dfo_GSMISSIONS
#' @param out_dir Directory to save the cleaned dataframe as an .rds file
#' 
#' @return A datafame with information of all unique DFO tows and matching column names to NMFS data. This file is also saved in out_dir.
#' 
#' @export
dfo_get_tows<-function(dfo_GSINF, dfo_GSMISSIONS, out_dir){
  
  # For debugging
  if(FALSE){
    dfo_GSINF = dfo_GSINF_load(here::here("scratch/aja/targets_flow/data/dfo/raw"))
    dfo_GSMISSIONS = dfo_GSMISSIONS_load(here::here("scratch/aja/targets_flow/data/dfo/raw"))
    out_dir = here::here("scratch/aja/targets_flow/data/dfo/clean")
  }
  
  # Join GSING with GSMISSIONS, create a unique ID and get DATE info
  temp_tows<- dfo_GSINF %>% 
    dplyr::left_join(dfo_GSMISSIONS, by = "MISSION") %>% # join with missions df
    mutate(ID = paste(MISSION, SETNO,sep="")) %>% # create a unique ID
    dplyr::filter(!is.na(SDATE)) %>% # only 2 NA SDATEs 
    tidyr::separate(SDATE, into = c("DATE", "EST_TIME"), sep = " ", remove = F) %>%
    filter(!is.na(LATITUDE) & !is.na(LONGITUDE))
  
  # Add US-equivalent date columns based on SDATE
  temp_tows$DATE<- lubridate::as_date(temp_tows$DATE)
  temp_tows$EST_YEAR<- lubridate::year(temp_tows$SDATE)
  temp_tows$EST_MONTH<- lubridate::month(temp_tows$SDATE)
  
  # Not all "SEASON" data is present, so classify into correct season based on month, where SPRING = Jan-Apr, SUMMER = May-Aug, FALL = Sept-Dec
  temp_tows$SEASON = NA
  temp_tows$SEASON[temp_tows$EST_MONTH <=4] = "SPRING"
  temp_tows$SEASON[temp_tows$EST_MONTH >= 5 & temp_tows$EST_MONTH <= 8] = "SUMMER"
  temp_tows$SEASON[temp_tows$EST_MONTH >= 9] = "FALL"
  
  # Create SVVESSEL column
  temp_tows$SVVESSEL<-as.factor(temp_tows$VESEL)

  # Rename and select columns to make US-equivalent
  # in rename function, first name is new name, second name is old name
  # select function = selects variables to keep and puts them in a US-equivalent order
  dfo_tows_out<- temp_tows %>% 
    dplyr::rename("DECDEG_BEGLAT" = "LATITUDE", "DECDEG_BEGLON" = "LONGITUDE") %>% 
    dplyr::select(ID, DATE, EST_YEAR, SEASON, SVVESSEL, DECDEG_BEGLAT, DECDEG_BEGLON)
  
  # Return and save it
  saveRDS(dfo_tows_out, file = paste(out_dir, "dfo_tows.rds", sep = "/"))
  return(dfo_tows_out)
  
  ## End function
}

#' @title DFO "tidy" occupancy
#' 
#' @description This function uses the `dfo_GSCAT` information to create a "tidy" occupancy dataset, where each row represents a unique tow - species - occupancy record for every species listed in `species_table`. In creating this "tidy" dataset, we impute "absence" observations, such that there is a record for every species in `species_table` at every unique tow.  
#' 
#' @param dfo_GSCAT 
#' @param dfo_tows
#' @param species_table = Dataframe for the species of interest, including their DFO and NMFS species codes
#' @param out_dir = Directory to save the processed dataframe as an .rds file
#' 
#' @return A "tidy" dataframe of species occurrences. This file is also saved in out_dir.
#' 
#' @export
dfo_make_tidy_occu<-function(dfo_GSCAT, dfo_tows, species_table, out_dir){
  
  # For debugging
  if(FALSE){
    dfo_GSCAT = dfo_GSCAT_load(here::here("scratch/aja/targets_flow/data/dfo/raw"))
    dfo_tows = dfo_GSCAT_load(here::here("scratch/aja/targets_flow/data/dfo/raw"))
    species_table = species_read_csv(here::here("scratch/aja/targets_flow/data/supporting"))
    out_dir = here::here("scratch/aja/targets_flow/data/dfo/clean")
  }
  
  # Create a long dataframe containing biomass and abundance data for all ID/species from GSCAT; this should be a presence only dataset
  presence_data<- dfo_GSCAT %>%
    mutate(ID = paste(MISSION, SETNO, sep="")) %>% # create a unique ID
    rename(., DFO_SPEC = SPEC) %>% # renaming to keep DFO/NMFS codes clear
    dplyr::filter(DFO_SPEC %in% c(species_table$DFO_SPEC)) %>%  # keep only species in species_table
    mutate(BIOMASS = ifelse(TOTWGT == 0 & TOTNO > 0, 0.001, TOTWGT)) %>%  # if TOTNO/ABUNDANCE > 0 but TOTWGT is 0, make TOTWGT non-zero (set equal to 0.001)
    rename("ABUNDANCE" = "TOTNO") %>%   # rename TOTNO to ABUNDANCE
    mutate(PRESENCE = ifelse(ABUNDANCE > 0, 1, 0)) %>%  # PRESENCE = 1 if ABUNDANCE >=1, PRESENCE = 0 if ABUNDANCE = 0   
    dplyr::select(ID, DFO_SPEC, PRESENCE, BIOMASS, ABUNDANCE) # keep only cols of interest                           
  
  # Create a dataframe of all possible survey ID/species combinations
  all_ID_SPEC_possibilities<- tibble::tibble(ID = rep(dfo_tows$ID, length(unique(species_table$DFO_SPEC)))) %>% 
    dplyr::arrange(ID) %>% 
    mutate(DFO_SPEC = rep(unique(species_table$DFO_SPEC), length(dfo_tows$ID))) %>%
    left_join(., species_table, by = "DFO_SPEC")
  
  # Create full presence absence dataset
  dfo_tidy_occu<- all_ID_SPEC_possibilities %>% 
    dplyr::left_join(presence_data, by = c("ID", "DFO_SPEC")) %>%                           
    # populate "possibilities" dataset with presence data 
    mutate(PRESENCE = ifelse(is.na(PRESENCE) == T, 0, PRESENCE)) %>%      
    mutate(BIOMASS = ifelse(is.na(BIOMASS) == T, 0.000, BIOMASS)) %>%    
    mutate(ABUNDANCE = ifelse(is.na(ABUNDANCE) == T, 0, ABUNDANCE)) %>% 
    dplyr::select(ID, NMFS_SVSPP, DFO_SPEC, PRESENCE, BIOMASS, ABUNDANCE) # keep only cols of interest       
  
  # Return and save it
  saveRDS(dfo_tidy_occu, file = paste(out_dir, "dfo_tidy_occu.rds", sep = "/"))
  return(dfo_tidy_occu)
 
  ## End function
  
}
