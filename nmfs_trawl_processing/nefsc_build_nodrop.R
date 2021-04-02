

#### Survdat Data Prep No Column Drop  ####
#' @title  Load survdat file with standard data filters, keep all columns
#' 
#'
#' @description Processing function to prepare survdat data for size spectra analyses. 
#' Options to select various survdat pulls, or provide your own available.
#' 
#'
#' @param survdat optional candidate dataframe in the R environment to run through size spectra build.
#' @param survdat_source String indicating which survdat file to load from box
#'
#' @return Returns a dataframe filtered and tidy-ed for size spectrum analysis.
#' @export
#'
#' @examples
survdat_prep_nodrop <- function(survdat = NULL, survdat_source = "2020"){
  
  ####  Resource Paths  
  box_paths   <- research_access_paths(os.use = "unix")
  mills_path  <- box_paths$mills
  res_path    <- box_paths$res
  caccel_path <- paste0(mills_path, "Projects/NSF_CAccel/Data/")
  
  
  ####  Import supplemental data  ####
  
  # # Load Survey Stratum Area Info - old and in nautical miles (with errors)
  # stratum_area <- read_csv(str_c(caccel_path, "strata area.csv"), col_types = cols())
  
  # Load fresh stratum areas
  stratum_area <- read_csv(str_c(caccel_path, "strata_areas_km2.csv"), col_types = cols())  %>% 
    mutate(stratum = as.character(stratum))
  
  # EPU are info: source slucey_survdat_functions.R and {ecodata}
  epu_areas <- read_csv(str_c(caccel_path, "EPU_areas_km2.csv"), col_types = cols())
  
  
  
  
  ####  Import SURVDAT Data  ####
  
  # Testing:
  #survdat_source <- "2016"        ; survdat <- NULL
  #survdat_source <- "2019"        ; survdat <- NULL
  #survdat_source <- "2020"        ; survdat <- NULL
  #survdat_source <- "2021"        ; survdat <- NULL
  #survdat_source <- "bigelow"     ; survdat <- NULL
  #survdat_source <- "most recent" ; survdat <- NULL
  
  # convenience change to make it lowercase
  survdat_source <- tolower(survdat_source)
  
  
  # Build Paths to survdat for standard options
  survdat_path <- switch(EXPR = survdat_source,
                         "2016"        = paste0(mills_path, "Projects/WARMEM/Old survey data/Survdat_Nye2016.RData"),
                         "2019"        = paste0(res_path,   "NMFS_trawl/Survdat_Nye_allseason.RData"),
                         "2020"        = paste0(res_path,   "NMFS_trawl/Survdat_Nye_Aug 2020.RData"),
                         "2021"        = paste0(res_path,   "NMFS_trawl/survdat_slucey_01152021.RData"),
                         "bigelow"     = paste0(res_path,   "NMFS_trawl/survdat_Bigelow_slucey_01152021.RData"),
                         "most recent" = paste0(res_path, "NMFS_trawl/NEFSC_BTS_all_seasons_03032021.RData") )
  
  
  
  # If providing a starting point for survdat pass it in:
  if(is.null(survdat) == FALSE){ 
    trawldat <- survdat %>% clean_names() 
  } else if(is.null(survdat) == TRUE){
    
    # If not then load using the correct path
    load(survdat_path)
    
    
    # Bigelow data doesn't load in as "survdat"
    if(survdat_source == "bigelow"){
      survdat <- survdat.big
      rm(survdat.big)}
    
    # Most recent pulls load a list containing survdat
    if(survdat_source == "most recent"){
      survdat <- survey$survdat }
    
    # clean names up for convenience
    trawldat <- survdat %>% clean_names() 
  }
  
  # remove survdat once the data is in
  rm(survdat)
  
  
  
  ####_ 1.  Special Steps for Different SURVDAT versions  ####
  
  ####__ a. Missing data flags  ####
  
  # Flags for missing columns that need to be merged in or built
  has_comname  <- "comname" %in% names(trawldat)
  has_id_col   <- "id" %in% names(trawldat)
  has_towdate  <- "est_towdate" %in% names(trawldat)
  has_month    <- "est_month" %in% names(trawldat)
  
  # Flags for renaming or subsetting the data due to presence/absence of columns
  has_year      <- "est_year" %in% names(trawldat)
  has_catchsex  <- "catchsex" %in% names(trawldat)
  has_decdeg    <- "decdeg_beglat" %in% names(trawldat)
  has_avg_depth <- "avgdepth" %in% names(trawldat)
  
  
  
  
  
  
  
  ####__ b. Missing comname  ####
  
  # Use SVSPP to get common names for species
  if(has_comname == FALSE){
    message("no comnames found, merging records in with NMFS_trawl/spp_keys/sppclass.csv")
    # Load sppclass codes and common names
    spp_classes <- read_csv(
      paste0(res_path, "NMFS_trawl/spp_keys/sppclass.csv"), 
      col_types = cols()) %>%
      clean_names() %>%
      mutate(comname         = str_to_lower(common_name),
             scientific_name = str_to_lower(scientific_name)) %>%
      distinct(svspp, comname, scientific_name)
    
    
    # Add the common names over and format for rest of build
    trawldat <- mutate(trawldat, svspp = str_pad(svspp, 3, "left", "0")) %>% 
      left_join(spp_classes, by = "svspp") 
    
  }
  
  
  ####__ c. Missing ID  ####
  if(has_id_col == FALSE) {
    message("creating station id from cruise-station-stratum fields")
    trawldat <- trawldat %>%
      
      # Build ID column
      mutate(cruise6 = str_pad(cruise6, 6, "left", "0"),
             station = str_pad(station, 3, "left", "0"),
             stratum = str_pad(stratum, 4, "left", "0"),
             id = str_c(cruise6, station, stratum))}
  
  
  ####__ d. Field renaming  ####
  
  # Rename select columns for consistency
  if(has_year == FALSE)      {
    message("renaming year column to est_year")
    trawldat <- rename(trawldat, est_year = year) }
  if(has_decdeg == FALSE) {
    message("renaming lat column to decdeg_beglat")
    trawldat <- rename(trawldat, decdeg_beglat = lat) }
  if(has_decdeg == FALSE) {
    message("renaming lon column to decdeg_beglon")
    trawldat <- rename(trawldat, decdeg_beglon = lon) }
  if(has_avg_depth == FALSE)      {
    message("renaming depth column to avgdepth")
    trawldat <- rename(trawldat, avgdepth = depth) }
  
  
  
  ####____ d. build date structure for quick grab of date components
  if(has_towdate == TRUE) {
    message("building month/day columns from est_towdate")
    trawldat <- mutate(trawldat,
                       est_month = str_sub(est_towdate, 6,7),
                       est_month = as.numeric(est_month),
                       est_day   = str_sub(est_towdate, -2, -1),
                       est_day   = as.numeric(est_day))}
  
  
  
  ####_ 2. Column Changes  ####
  trawldat <- trawldat %>% 
    mutate(
      
      # Text Formatting 
      comname = tolower(comname),
      id      = format(id, scientific = FALSE),
      svspp   = as.character(svspp),
      svspp   = str_pad(svspp, 3, "left", "0"),
      
      # Stratum number, 
      # exclude leading and trailing codes for inshore/offshore, 
      # used for matching to stratum areas
      strat_num = str_sub(stratum, 2, 3)) %>%  
    
    # Replace NA's where there is some biomass/abundance
    mutate(biom_adj  = ifelse(biomass == 0 & abundance > 0, 0.0001, biomass), 
           .after = biomass) %>% 
    mutate(abund_adj = ifelse(abundance == 0 & biomass > 0, 1, abundance), 
           .after = abundance)
  
  
  
  

  
  ####_ 3. Row Filtering  ####
  
  # Things filtered:
  # 1. Strata
  # 2. Seasons
  # 3. Year limits
  # 4. Vessels
  # 5. Species Exclusion
  trawldat <- trawldat %>%
    filter(
      
      # Eliminate Canadian Strata and Strata No longer in Use 
      stratum >= 01010,
      stratum <= 01760,
      stratum != 1310,
      stratum != 1320,
      stratum != 1330,
      stratum != 1350,
      stratum != 1410,
      stratum != 1420,
      stratum != 1490,
      
      # Filter to just Spring and Fall
      season %in% c("SPRING", "FALL"),
      
      # Only the Albatross and Henry Bigelow
      #svvessel %in% c("AL", "HB"),
      est_year >= 1970,
      est_year < 2020,
      
      # Drop NA Biomass and Abundance Records
      !is.na(biomass),
      !is.na(abundance),
      
      # Exclude the Skrimps
      svspp %not in% c(285:299, 305, 306, 307, 316, 323, 910:915, 955:961),
      
      # Exclude the unidentified fish
      svspp %not in% c(0, 978, 979, 980, 998)
    )
  
  
  
  
  
  
  
  
  ####_ 4. Spatial Filtering  ####
  
  # This section assigns EPU's via overlay using Sean Lucey's code,
  # And also merges with stratum area information,
  # these are used to relate catch/effort to physical areas in km squared
  
  #### EPU assignment for survdat stations - function from Sean Luceys "RSurvey" repo
  
  # Area stratification code from Sean Lucey's Repo, renamed to signify what it is
  source(paste0(res_path, "NMFS_trawl/slucey_functions/slucey_survdat_functions.R"))
  
  
  # EPU shapefiles can be loaded from ecodata package
  epu_sf <- ecodata::epu_sf
  epu_sp <- suppressWarnings(as_Spatial(epu_sf))
  
  
  # Rename columns to match expected formats for Sean Lucey's poststrat()
  trawldat <- trawldat %>% 
    rename(CRUISE6 = cruise6,
           STATION = station,
           STRATUM = stratum,
           LAT     = decdeg_beglat,
           LON     = decdeg_beglon)
  
  # Post stratify the station positions using EPU polygons
  trawldat <- trawldat %>% 
    poststrat(survdat = ., stratum = epu_sp, strata.col = "EPU") %>% 
    rename(epu           = newstrata,
           cruise6       = CRUISE6,
           station       = STATION,
           stratum       = STRATUM,
           decdeg_beglat = LAT,
           decdeg_beglon = LON)
  
  
  # Stratum Area Key for which stratum correspond to larger regions we use
  strata_key <- list(
    "Georges Bank"          = as.character(13:23),
    "Gulf of Maine"         = as.character(24:40),
    "Southern New England"  = str_pad(as.character(1:12), width = 2, pad = "0", side = "left"),
    "Mid-Atlantic Bight"    = as.character(61:76))
  
  
  # Add the labels to the data
  trawldat <- trawldat %>%
    mutate(
      survey_area =  case_when(
        strat_num %in% strata_key$`Georges Bank`         ~ "GB",
        strat_num %in% strata_key$`Gulf of Maine`        ~ "GoM",
        strat_num %in% strata_key$`Southern New England` ~ "SNE",
        strat_num %in% strata_key$`Mid-Atlantic Bight`   ~ "MAB",
        TRUE                                             ~ "stratum not in key"))
  
  
  # Optional, Use strata_select to pull the strata we want individually
  strata_select <- c(strata_key$`Georges Bank`, strata_key$`Gulf of Maine`,
                     strata_key$`Southern New England`, strata_key$`Mid-Atlantic Bight`)
  
  
  # Filtering with strata_select
  trawldat <- trawldat %>% filter(strat_num %in% strata_select) %>% 
    mutate(stratum = as.character(stratum))
  
  
  
  ####_ 5. Stratum Area/Effort Ratios  ####
  # Stratum area ratio is the ratio between the area of the select 
  # stratum to the total area of all stratum sampled that year
  
  # Join to the files containing area of each stratum, epu
  trawldat <- trawldat %>% 
    left_join(stratum_area, by = "stratum") %>% 
    left_join(epu_areas, by = "epu") %>% 
    arrange(trawldat, id)
  
  
  
  
  # Get Total stratum area of all strata sampled in a year
  # (excludes ones we do not care about via left join)
  total_stratum_areas <- trawldat %>% 
    group_by(est_year) %>% 
    distinct(stratum, .keep_all = T) %>%  
    summarise(tot_s_area =  sum(s_area_km2, na.rm = T),
              .groups = "keep") %>% ungroup()
  
  # Get total area of EPU's sampled in a year
  total_epu_areas <- trawldat %>% 
    group_by(est_year) %>% 
    distinct(epu, .keep_all = T) %>%  
    summarise(tot_epu_area =  sum(epu_area_km2, na.rm = T), 
              .groups = "keep") %>% ungroup()
  
  
  # Calculate strata area relative to total area i.e. stratio or stratum weights
  trawldat <- trawldat %>% 
    left_join(total_stratum_areas, by = "est_year") %>% 
    left_join(total_epu_areas, by = "est_year") %>% 
    mutate(st_ratio   = s_area_km2 / tot_s_area,
           epu_ratio  = epu_area_km2 / tot_epu_area) 
  
  # We have to total areas, now we want effort in each
  # Number of unique tows per stratum
  yr_strat_effort <- trawldat %>% 
    group_by(est_year, stratum) %>% 
    summarise(strat_ntows = n_distinct(id), 
              .groups = "keep") %>% ungroup()
  
  # Number of unique tows per EPU
  yr_epu_effort <-  trawldat %>% 
    group_by(est_year, epu) %>% 
    summarise(epu_ntows = n_distinct(id), 
              .groups = "keep") %>% ungroup()
  
  
  # Add those yearly effort counts back for later 
  #(area stratified abundance)
  trawldat <- trawldat %>% 
    left_join(yr_strat_effort, by = c("est_year", "stratum")) %>% 
    left_join(yr_epu_effort, by = c("est_year", "epu"))
  
  
  
  
  
  ####_ 6. Adjusting NumLength to Match Abundance  ####
  # Sometimes there are more/less measured than initially tallied* in abundance
  
  
  # If catchsex is not a column then total abundance is assumed pooled
  if(has_catchsex == TRUE){
    abundance_groups <- c("id", "comname", "catchsex")
  } else {
    message("catchsex column not found, ignoring sex for numlen adjustments")
    abundance_groups <- c("id", "comname")}
  
  
  
  # Get the abundance value for each sex 
  # arrived at by summing across each length
  abundance_check <- trawldat %>%
    group_by(!!!syms(abundance_groups), abundance) %>%
    summarise(
      abund_actual = sum(numlen),               
      n_len_class  = n_distinct(length),
      .groups      = "keep") %>% ungroup()
  
  # Get the ratio between the original abundance column 
  # and the sum of numlen we just grabbed
  conv_factor <- trawldat %>% 
    #distinct(id, comname, catchsex, length, abund_adj) %>% 
    distinct(!!!syms(abundance_groups), length, abund_adj) %>% 
    #inner_join(abundance_check, by = c("id", "comname", "catchsex")) %>% 
    inner_join(abundance_check, by = abundance_groups) %>% 
    mutate(convers = abund_adj / abund_actual)
  
  
  
  # Merge back and convert the numlen field
  # original numlen * conversion factor = numlength adjusted
  survdat_processed <- trawldat %>%
    #left_join(conv_factor, by = c("id", "comname", "catchsex", "length", "abundance", "abund_adj")) %>%
    left_join(conv_factor, by = c(abundance_groups, "length", "abundance", "abund_adj")) %>%
    mutate(numlen_adj = numlen * convers, .after = numlen) %>% 
    select(-c(abund_actual, convers))
  
  # remove conversion factor from environment
  rm(abundance_check, conv_factor, strata_key, strata_select, stratum_area, epu_areas, epu_sf, epu_sp)
  
  
  
  
  
  ####_ 7. Distinct Station & Species Length Info   ####
  
  # For each station we need unique combinations of
  # station_id, species, catchsex, length, adjusted_numlen
  # to capture what and how many of each length fish is caught
  
  # Record of unique station catches: 
  # One row for every species * sex * length, combination in the data
  trawl_lens <- survdat_processed %>% 
    filter(is.na(length) == FALSE,
           is.na(numlen) == FALSE,
           numlen_adj > 0) 
  
  
  # Do we want to just keep all the station info here as well?
  # question to answer is whether any other columns repeat,
  # or if these are the only ones
  trawl_spectra <- trawl_lens %>%
    distinct(id, svspp, comname, catchsex, abundance, n_len_class,
             length, numlen, numlen_adj, biom_adj, .keep_all = TRUE)
  
  
  
  
  
  # Reorder the main columns
  
  
  
  
  
  # Return the dataframe
  # Row for each length class of every species caught
  return(trawl_spectra)
  
}