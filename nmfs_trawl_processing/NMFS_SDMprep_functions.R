## Functions for processing DFO bottom trawl surveys for SDM workflow

# Setup ---------------------------------------------------------
# Detect the operating system
os.use<- .Platform$OS.type

# Set path to shared folders

#########CHANGE THIS COMPUTER NAME IF PC
computer.name<- "lcarlson" # Needed for PC users

shared.path<- switch(os.use, 
                     "unix" = paste("~/Box/", user.name, sep = ""),
                     "windows" = paste("C:/Users/", computer.name, "/Box/", sep = ""))






# Helper Function ---------------------------------------------------------
library_check<- function(libraries) {
  ## Details
  # This function will check for and then either load or install any libraries needed to run subsequent functions
  
  # Args:
  # libraries = Vector of required library names
  
  # Returns: NA, just downloads/loads libraries into current work space.
  
  ## Start function
  lapply(libraries, FUN = function(x) {
    if (!require(x, character.only = TRUE)) {
      install.packages(x, dependencies = TRUE)
      library(x, character.only = TRUE)
    }
  })
  ## End function
}










# Demographic Data Function ---------------------------------------------------------
create_alltows_nmfs<-function(datapath){
  
  # Details
  
  # This function reads in preprocessed NMFS bottom trawl survey data. The function creates a dataframe that contains demographic data about all unique tow IDs
  
  # Args:
  # datapath = path to pre-processed NMFS datafile (output file will be in the same location)
  
  
  # Returns:
  # RDS Data file, which is also saved in folder specified by datapath
  
  
  ## Start function
  # Preliminaries ----------------------
  
  libraries.needed<- c("tidyverse", "lubridate")
  library_check(libraries.needed)
  
  
  # Load in data
  
  NMFS_input<-read_csv(paste(datapath,"march_2021_survdat_allcols.csv",sep = "")) 
  
  NMFS_alltows<-NMFS_input %>% 
    dplyr::mutate(date = lubridate::date(est_towdate)) %>% 
    dplyr::distinct(id,date,est_year,season,svvessel,decdeg_beglat,decdeg_beglon)
  
  # write data file
  
  saveRDS(NMFS_alltows, file = paste(datapath, "NMFS_alltows.dat", sep = "/"))
  
  ## End function
  
}


#example function call: create_alltows_nmfs(datapath = paste(shared.path,"RES_Data/NMFS_trawl/processed_data/",sep=""))












# Presence/Absence, Biomass, Abundance Function ---------------------------------------------------------
create_tidyoccurence_nmfs<-function(datapath){
  
  # Details
  
  # This function reads in preprocessed NMFS bottom trawl survey data as well as the processed demographic "alltows" dataset (created with create_alltows_nmfs) and a supporting dataset (species_conversion). The function then creates a long dataset containing presence/absence for all possible combinations of survey ID and species; plus associated biomass and abundance data
  
  # Args:
  # datapath = path to pre-processed NMFS datafile (output file will be in the same location)
  
  
  # Returns:
  # RDS Data file, which is also saved in folder specified by datapath
  
  
  ## Start function
  # Preliminaries ----------------------------
  
  libraries.needed<- c("tidyverse", "lubridate")
  library_check(libraries.needed)
  
  
  # Load in data
  
  NMFS_input<-read_csv(file=paste(datapath,"march_2021_survdat_allcols.csv",sep = "")) 
  NMFS_alltows<-readRDS(file=paste(datapath, "NMFS_alltows.dat", sep = "/"))
  species_conversion<-read_csv(file=paste(shared.path,"Mills Lab/Projects/DFO_survey_data/supporting_data/species_naming_conversion.csv",sep = ""))
  
  
  # create a long dataframe containing biomass and abundance data for all ID/species; this should be a presence only dataset
  
  nmfs_presence_data<-NMFS_input %>%   
    dplyr::mutate(svspp = as.double(svspp)) %>% 
    dplyr::filter(svspp %in% c(species_conversion$SVSPP)) %>%                      #keep only Shackell species
    dplyr::group_by(id, svspp) %>% 
    dplyr::summarise(abundance = sum(abund_adj), biomass = sum(biom_adj)) %>% 
    dplyr::mutate(presence = ifelse(abundance > 0, 1, 0)) %>%                      #should all be 1s
    #presence = 1 if abundance >=1, presence = 0 if abundance = 0
    dplyr::select(id, svspp ,presence, biomass, abundance)      
  
  
  # create a dataframe of all possible survey ID/species (svspp) combinations
  
  all_id_svspp_possibilities<-tibble::tibble(id = rep(NMFS_alltows$id,length(unique(species_conversion$SVSPP)))) %>% 
    dplyr::arrange(id) %>% 
    dplyr::mutate(svspp = rep(unique(species_conversion$SVSPP),length(unique(NMFS_alltows$id))))
  
  
  # create full presence absence dataset
  
  NMFS_tidyoccurence<-all_id_svspp_possibilities %>% 
    dplyr::left_join(nmfs_presence_data, by = c("id","svspp")) %>%                           
    #populate "possibilities" dataset with presence data 
    dplyr::mutate(presence = ifelse(is.na(presence) == T, 0,presence)) %>%      #populate with 0s
    dplyr::mutate(biomass = ifelse(is.na(biomass) == T, 0.0000,biomass)) %>%     #populate with 0s
    dplyr::mutate(abundance = ifelse(is.na(abundance) == T, 0,abundance)) %>%   #populate with 0s
    dplyr::select(id,svspp,presence,biomass,abundance)   
  
  
  # write data file
  
  saveRDS(NMFS_tidyoccurence, file = paste(datapath, "NMFS_tidyoccurence.dat", sep = ""))
  
  ## End function
  
}


#example function call:   create_tidyoccurence_nmfs(datapath = paste(shared.path,"RES_Data/NMFS_trawl/processed_data/",sep=""))



