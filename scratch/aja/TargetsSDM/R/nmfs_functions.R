####  Common Resources  ####


####  Functions  ####
####
## Dynamic loading functions -- each takes in a directory path and reads in files from it. 
nmfs_load<- function(nmfs_raw_dir){
  # Read and return data
  nmfs_raw<- read.csv(paste(nmfs_raw_dir, "march_2021_survdat_allcols.csv", sep = "/"))
  return(nmfs_raw)
}

## Cleaning and processing functions
#' @title NMFS tow data
#' 
#' @description This function reads in the nmfs_raw data (originally loaded as `survdat` and then renamed in the loading `nmfs_load` to `nmfs_raw`) and then creates a dataframe that has one row for every unique tow
#' 
#' @param nmfs_raw Raw NMFS data, loaded with `nmfs_load`
#' @param out_dir Directory to save the cleaned dataframe as an .rds file
#' 
#' @return A datafame with information of all unique NMFS tows. This file is also saved in out_dir.
#' 
#' @export
nmfs_get_tows<-function(nmfs_raw, out_dir){
  
  # For debugging
  if(FALSE){
    nmfs_raw = nmfs_load(here::here("scratch/aja/targets_flow/data/nmfs/raw"))
    out_dir = here::here("scratch/aja/targets_flow/data/nmfs/clean")
  }
  
  # Get unique tows
  temp_tows<- nmfs_raw %>%
    dplyr::mutate(date = lubridate::date(est_towdate)) %>% 
    dplyr::distinct(id, date, est_year, season, svvessel, decdeg_beglat, decdeg_beglon) %>%
    filter(!is.na(decdeg_beglat) & !is.na(decdeg_beglon)) %>% # Remove NA locations
    mutate(., id = as.character(id))
  
  # Renaming...revist this
  nmfs_tows_out<- temp_tows
  colnames(nmfs_tows_out)<- toupper(colnames(nmfs_tows_out))

  # Return and save it
  saveRDS(nmfs_tows_out, file = paste(out_dir, "nmfs_tows.rds", sep = "/"))
  return(nmfs_tows_out)
}


#' @title NMFS "tidy" occupancy
#' 
#' @description This function uses the raw NMFS data (`nmfs_raw`) to create a "tidy" occupancy dataset, where each row represents a unique tow - species - occupancy record for every species listed in `species_table`. In creating this "tidy" dataset, we impute "absence" observations, such that there is a record for every species in `species_table` at every unique tow.  
#' 
#' @param nmfs_raw = Raw NMFS data, loaded with `nmfs_load`
#' @param nmfs_tows = 
#' @param species_table = Dataframe for the species of interest, including their DFO and NMFS species codes
#' @param out_dir = Directory to save the processed dataframe as an .rds file
#' 
#' @return A "tidy" dataframe of species occurrences. This file is also saved in out_dir.
#' 
#' @export
nmfs_make_tidy_occu<-function(nmfs_raw, nmfs_tows, species_table, out_dir){
  
  # For debugging
  if(FALSE){
    nmfs_raw = nmfs_load(here::here("scratch/aja/targets_flow/data/nmfs/raw"))
    nmfs_tows = nmfs_load(here::here("scratch/aja/targets_flow/data/nmfs/raw"))
    species_table = species_read_csv(here::here("scratch/aja/targets_flow/data/supporting"))
    out_dir = here::here("scratch/aja/targets_flow/data/nmfs/clean")
  }
  
  # Renaming columns...
  colnames(nmfs_raw)<- toupper(colnames(nmfs_raw))
  colnames(nmfs_tows)<- toupper(colnames(nmfs_tows))
  
  # Create a long dataframe containing biomass and abundance data for all ID/species; this should be a presence only dataset
  presence_data<- nmfs_raw %>% 
    dplyr::mutate(SVSPP = as.double(SVSPP)) %>% 
    dplyr::filter(SVSPP %in% c(species_table$NMFS_SVSPP)) %>% #keep only Shackell species
    dplyr::group_by(ID, SVSPP) %>% 
    dplyr::summarise(ABUNDANCE = sum(ABUND_ADJ), BIOMASS = sum(BIOM_ADJ)) %>% 
    dplyr::mutate(PRESENCE = ifelse(ABUNDANCE > 0, 1, 0)) %>% #should all be 1s
    #presence = 1 if abundance >=1, presence = 0 if abundance = 0
    dplyr::select(ID, SVSPP, PRESENCE, BIOMASS, ABUNDANCE) %>%
    rename(., "NMFS_SVSPP" = SVSPP) %>%
    ungroup() %>%
    mutate(., ID = as.character(ID))
    
  # Create a dataframe of all possible survey ID/species combinations
  all_ID_SPEC_possibilities<- tibble::tibble(ID = rep(nmfs_tows$ID, length(unique(species_table$NMFS_SVSPP)))) %>% 
    dplyr::arrange(ID) %>% 
    mutate(NMFS_SVSPP = rep(unique(species_table$NMFS_SVSPP), length(nmfs_tows$ID))) %>%
    left_join(., species_table, by = "NMFS_SVSPP") %>%
    mutate(., ID = as.character(ID))
  
  # Create full presence absence dataset
  nmfs_tidy_occu<- all_ID_SPEC_possibilities %>% 
    dplyr::left_join(presence_data, by = c("ID", "NMFS_SVSPP")) %>%                           
    #populate "possibilities" dataset with presence data                       
    mutate(PRESENCE = ifelse(is.na(PRESENCE) == T, 0, PRESENCE)) %>%     
    mutate(BIOMASS = ifelse(is.na(BIOMASS) == T, 0.000, BIOMASS)) %>%  
    mutate(ABUNDANCE = ifelse(is.na(ABUNDANCE) == T, 0, ABUNDANCE)) %>%  
    dplyr::select(ID, NMFS_SVSPP, DFO_SPEC, PRESENCE, BIOMASS, ABUNDANCE) %>% #keep only cols of interest
    mutate(., ID = as.character(ID))

  # Return and save it
  saveRDS(nmfs_tidy_occu, file = paste(out_dir, "nmfs_tidy_occu.rds", sep = "/"))
  return(nmfs_tidy_occu)
  
  ## End function
}

