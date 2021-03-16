

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
create_demographic<-function(datapath,outpath){
  
  # Details
  
  # This function reads in DFO bottom trawl survey data from the Mar.datawrangling package (download using devtools). The function then merges the GSINF and GSMISSIONS dataframes and adds/renames/reorders columns to mirror NMFS survey data formatting
  
  # Args:
  # datapath = path to DFO datafile
  # outpath = path to save processed dataset
  
  
  # Returns:
  # RDS Data file, which is also saved in folder specified by outpath
  
  
  ## Start function
  # Preliminaries ----------------------
  
  libraries.needed<- c("tidyverse", "lubridate", "Mar.datawrangling")
  library_check(libraries.needed)
  
  
  # Load in data
  
  Mar.datawrangling::get_data(db='rv', data.dir = datapath)
  
  
  # will need GSINF dataframe from DFO bottom trawl survey (Maritimes package)
  # will need GSMISSIONS dataframe from DFO bottom trawl survey (Maritimes package)
  
  
  temp_demographic<-GSINF %>% 
    dplyr::left_join(GSMISSIONS, by = "MISSION") %>%                                #join with missions df
    dplyr::mutate(ID = paste(MISSION,SETNO,sep="")) %>%                             #create a unique ID
    dplyr::filter(!is.na(SDATE)) %>%                                                #only 2 NA SDATEs
    tidyr::separate(SDATE,into = c("DATE","EST_TIME"), sep = " ", remove=F)         #break SDATE: date and time
  
  
  # add US-equivalent date columns based on SDATE
  
  temp_demographic$DATE<- lubridate::as_date(temp_demographic$DATE)
  temp_demographic$EST_YEAR<- lubridate::year(temp_demographic$SDATE)
  temp_demographic$EST_MONTH<- lubridate::month(temp_demographic$SDATE)
  
  
  # not all "SEASON" data is present, so classify into correct season based on month
  # SPRING = Jan-Apr
  # SUMMER = May-Aug
  # FALL = Sept-Dec
  
  temp_demographic$SEASON = NA
  temp_demographic$SEASON[temp_demographic$EST_MONTH <=4] = "SPRING"
  temp_demographic$SEASON[temp_demographic$EST_MONTH >= 5 & temp_demographic$EST_MONTH <= 8] = "SUMMER"
  temp_demographic$SEASON[temp_demographic$EST_MONTH >= 9] = "FALL"
  
  
  # create SVVESSEL column
  
  temp_demographic$SVVESSEL<-as.factor(temp_demographic$VESEL)
  
  
  # rename and select columns to make US-equivalent
  # in rename function, first name is new name, second name is old name
  # select function = selects variables to keep and puts them in a US-equivalent order
  
  DFO_demographic<-temp_demographic %>% 
    dplyr::filter(EST_YEAR >= 1982) %>%                                   #data pre-1982 used different methods
    dplyr::rename("TOWDUR" = "DUR", "RPM" = "SPEED",  
                  "DECDEG_BEGLAT" = "LATITUDE", "DECDEG_BEGLONG" = "LONGITUDE") %>% 
    dplyr::select(ID, DATE, EST_YEAR, SEASON, SVVESSEL, TOWDUR, RPM, DECDEG_BEGLAT, DECDEG_BEGLONG)
  
  
  # write data file
  
  saveRDS(DFO_demographic, file = paste(outpath, "DFO_demographic.dat", sep = "/"))
  
  ## End function
  
}














# Presence/Absence, Biomass, Abundance Function ---------------------------------------------------------
create_presence_absence<-function(datapath,demopath,specpath,outpath){
  
  # Details
  
  # This function reads in DFO bottom trawl survey data from the Mar.datawrangling package (download using devtools) as well as the processed demographic dataset (created with create_demographic) and a supporting dataset (species_conversion). The function then creates a long dataset containing presence/absence for all possible combinations of survey ID and species; plus associated biomass and abundance data
  
  # Args:
  # datapath = path to DFO datafile
  # demopath = path to processed demographic dataset
  # specpath = path to supporting species dataset
  # outpath = path to save processed dataset
  
  # Returns:
  # RDS Data file, which is also saved in folder specified by outpath
  
  
  ## Start function
  # Preliminaries ----------------------------
  
  libraries.needed<- c("tidyverse", "lubridate", "Mar.datawrangling")
  library_check(libraries.needed)
  
  
  # Load in data
  
  Mar.datawrangling::get_data(db='rv', data.dir = datapath)
  DFO_demographic<-readRDS(file=paste(demopath,"DFO_demographic.dat",sep = "/"))
  species_conversion<-read_csv(file=paste(specpath,"species_naming_conversion.csv",sep = "/"))
  
  
  # create a long dataframe containing biomass and abundance data for all ID/species; this should be a presence only dataset
  
  presence_data<-GSCAT %>% 
    dplyr::mutate(ID = paste(MISSION,SETNO,sep="")) %>%                          #create a unique ID
    dplyr::filter(ID %in% c(DFO_demographic$ID)) %>%          
    #keep only IDs in demographic data (which is a complete list of surveys post 1982)
    dplyr::filter(SPEC %in% c(species_conversion$SPEC)) %>%                      #keep only Shackell species
    dplyr::mutate(BIOMASS = ifelse(TOTWGT == 0 & TOTNO > 0,0.001,TOTWGT)) %>%    
    #if TOTNO/ABUNDANCE > 0 but TOTWGT is 0, make TOTWGT non-zero (set equal to 0.001)
    dplyr::rename("ABUNDANCE" = "TOTNO") %>%                                     #rename TOTNO to ABUNDANCE
    dplyr::mutate(PRESENCE = ifelse(ABUNDANCE > 0, 1, 0)) %>%                     
    #PRESENCE = 1 if ABUNDANCE >=1, PRESENCE = 0 if ABUNDANCE = 0
    dplyr::select(ID,SPEC,PRESENCE, BIOMASS, ABUNDANCE)                          #keep only cols of interest
  
  
  # create a dataframe of all possible survey ID/species combinations
  
  all_ID_SPEC_possibilities<-tibble::tibble(ID = rep(DFO_demographic$ID,length(unique(species_conversion$SPEC)))) %>% 
    dplyr::arrange(ID) %>% 
    dplyr::mutate(SPEC = rep(unique(species_conversion$SPEC),length(unique(DFO_demographic$ID)))) 
  
  
  # create full presence absence dataset
  
  DFO_presabs<-all_ID_SPEC_possibilities %>% 
    dplyr::left_join(presence_data, by = c("ID","SPEC")) %>%                           
    #populate "possibilities" dataset with presence data 
    dplyr::left_join(species_conversion, by = "SPEC") %>%                       #add SVSPP (NMFS species ID)
    dplyr::mutate(PRESENCE = ifelse(is.na(PRESENCE) == T, 0,PRESENCE)) %>%      #populate with 0s
    dplyr::mutate(BIOMASS = ifelse(is.na(BIOMASS) == T, 0.000,BIOMASS)) %>%     #populate with 0s
    dplyr::mutate(ABUNDANCE = ifelse(is.na(ABUNDANCE) == T, 0,ABUNDANCE)) %>%   #populate with 0s
    dplyr::select(ID,SVSPP,PRESENCE,BIOMASS,ABUNDANCE)                          #keep only cols of interest
  
  # write data file
  
  saveRDS(DFO_presabs, file = paste(outpath, "DFO_presence_absence.dat", sep = "/"))
  
  ## End function
  
}
