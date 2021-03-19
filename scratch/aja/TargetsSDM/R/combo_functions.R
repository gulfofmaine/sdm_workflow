####  Common Resources  ####


####  Functions  ####
####
## Dynamic loading functions -- each takes in a directory path and reads in files from it. 
species_read_csv<- function(species_table_dir){
  # Read it
  species_table<- read.csv(paste(species_table_dir, "species_table.csv", sep = "/"))
  
  # Return it
  return(species_table)
}

## Cleaning and processing functions
#' @title Bind NMFS and DFO tow data 
#' 
#' @description A short function that will row bind a NMFS tows dataframe and DFO tows dataframe to create one dataframe that has all tows combined, with one row per unique tow.
#'
#' @param nmfs_tows = NMFS tows dataframe, created by `nmfs_get_tows()`
#' @param dfo_tows = DFO tows dataframe, created by `dfo_get_tows()`
#' @param out_dir = Directory to save the combined dataframe as an .rds file
#' 
#' @return A datafame with information of all unique tows. This file is also saved in out_dir. 
#' 
#' @export

bind_nmfs_dfo_tows<- function(nmfs_tows, dfo_tows, out_dir){
  
  # For debugging
  if(FALSE){
    nmfs_tows = nmfs_tows_out
    dfo_tows = dfo_tows_out
    out_dir = here::here("scratch/aja/targets_flow/data/combined")
  }
  
  # Load in both datasets, add survey column
  nmfs_tows<- nmfs_tows %>%
    mutate(., "SURVEY" = rep("NMFS", nrow(.)))
  dfo_tows<- dfo_tows %>%
    mutate(., "SURVEY" = rep("DFO", nrow(.)))
  
  # Bind them together
  all_tows<- bind_rows(nmfs_tows, dfo_tows) %>%
    arrange(., DATE)
  
  # Return and save
  saveRDS(all_tows, file = paste(out_dir, "all_tows.rds", sep = "/"))
  return(all_tows)
}

#' @title Bind NMFS and DFO "tidy" occupancy data
#' 
#' @description A short function that will row bind a NMFS tidy occupancy dataframe with a DFO tidy occupancy dataframe to create one data frame that has all combined observation, where each row is a unique tow - species - occupancy record for all species in `species_table`.
#'
#' @param nmfs_tidy_occu = NMFS tidy occupancy dataframe, created by `nmfs_make_tidy_occu()`
#' @param dfo_tidy_occu = DFO tidy occupancy dataframe, created by `dfo_make_tidy_occu()`
#' @param out_dir = Directory to save the combined dataframe as an .rds file
#' 
#' @return A datafame with information of all occupancy records. This file is also saved in out_dir. 
#' 
#' @export

bind_nmfs_dfo_tidy_occu<- function(nmfs_tidy_occu, dfo_tidy_occu, out_dir){
  
  # For debugging
  if(FALSE){
    nmfs_tidy_occu = nmfs_tidy_occu
    dfo_tidy_occu = dfo_tidy_occu
    out_dir = here::here("scratch/aja/targets_flow/data/combined")
  }
  
  # Load in both datasets, add survey column
  nmfs_tidy_occu<- nmfs_tidy_occu %>%
    mutate(., "SURVEY" = rep("NMFS", nrow(.)))
  dfo_tidy_occu<- dfo_tidy_occu %>%
    mutate(., "SURVEY" = rep("DFO", nrow(.)))
  
  # Bind them together
  all_tidy_occu<- bind_rows(nmfs_tidy_occu, dfo_tidy_occu)
  
  # Return and save
  saveRDS(all_tidy_occu, file = paste(out_dir, "all_tidy_occu.rds", sep = "/"))
  return(all_tidy_occu)
}

#' @title Make a "tidy" model dataframe
#' 
#' @description This function brings together information on the unique tows, which will usually include measures of effort as well as habitat covariates, with the occupancy data to create a "tidy" model dataframe. The "tidy" model dataframe is set up so that each row is a unique tow - species - occupancy record with additional information for habitat covariates that we might include in a species distribution model. 
#'
#' @param all_tows = A dataframe with information on all unique tows, including any potential habitat covariates describing characteristics of these locations
#' @param all_tidy_occu = A dataframe with all occupancy records, created by `bind_nmfs_dfo_tidy_occu`
#' @param out_dir = Directory to save the tidy model dataframe as an .rds file
#' 
#' @return A tidy model datafame with all the information (tows, habitat covariates, species occurrences) needed to fit a species distribution model. This file is also saved in out_dir. 
#' 
#' @export
make_tidy_mod_data<- function(all_tows, all_tidy_occu, out_dir){
  
  # For debugging
  if(FALSE){
    all_tows = all_tows
    all_tidy_occu = all_tidy_occu 
    out_dir = here::here("scratch/aja/targets_flow/data/combined")
  }
  
  # Need to join up the tow info with the tidy occu
  tidy_mod_data<- all_tidy_occu %>%
    left_join(., all_tows, by = c("ID", "SURVEY"))
  
  # Keep only what we need..
  tidy_mod_data_out<- tidy_mod_data %>%
    dplyr::select(., ID, EST_DATE, EST_YEAR, SEASON, SURVEY, SVVESSEL, DECDEG_BEGLAT, DECDEG_BEGLON, NMFS_SVSPP, DFO_SPEC, PRESENCE, BIOMASS, ABUNDANCE) 
  
  # Return and save
  saveRDS(tidy_mod_data_out, file = paste(out_dir, "tidy_mod_data.rds", sep = "/"))
  return(tidy_mod_data_out)
}

# Function to convert points to sf object
points_to_sf<- function(points){
  sf_out<- st_as_sf(points, coords = c("DECDEG_BEGLON", "DECDEG_BEGLAT"), crs = 4326, remove = FALSE)
  return(sf_out)
}