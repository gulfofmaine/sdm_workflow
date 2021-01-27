####  A. Kemberling  ####
# Trawl Cleaning Scripts

####  Packages  ####
library(devtools)
library(gmRi)
library(tidyverse)
library(sf)


# Load resource paths
box_paths <- research_access_paths(os.use = "unix")
mills_path <- shared.path("unix", "mills", "")

####  Load Cleanup Functions  ####
source_url("https://raw.githubusercontent.com/adamkemberling/nefsc_trawl/master/R/01_nefsc_ss_build.R")

# Just copy and paste this way to see the most current version... Obviously not ideal
# Really just didn't want to develop it in three different locations
survdat_prep


####_______________________________####


#### SURVDAT Cleanup  ####
# copied: 1/14/2021
function(survdat = NULL, survdat_source = "2020"){
  
  ####  Resource Paths  
  box_paths   <- research_access_paths(os.use = "unix")
  mills_path  <- box_paths$mills
  res_path    <- box_paths$res
  caccel_path <- paste0(mills_path, "Projects/NSF_CAccel/Data/")
  
  
  ####  Load supplemental data  ####
  
  
  
  # # Load Survey Stratum Area Info - old and in nautical miles (with errors)
  # stratum_area <- read_csv(str_c(caccel_path, "strata area.csv"), col_types = cols())
  
  # Load fresh stratum areas
  stratum_area <- read_csv(str_c(caccel_path, "strata_areas_km2.csv"), 
                           col_types = cols())  %>% 
    mutate(stratum = as.character(stratum))
  
  # EPU are info: source slucey_survdat_functions.R and {ecodata}
  epu_areas <- read_csv(str_c(caccel_path, "EPU_areas_km2.csv"), col_types = cols())
  
  
  
  
  ####  Load SURVDAT Data  ####
  
  #convenience change
  survdat_source <- tolower(survdat_source)
  
  # Testing:
  #survdat_source <- "2016"
  #survdat_source <- "2019"
  #survdat_source <- "2020"
  #survdat_source <- "2021"
  #survdat_source <- "bigelow"
  
  
  # Build Paths to survdat for standard options
  survdat_path <- switch(EXPR = survdat_source,
                         "2016"    = paste0(mills_path, "Projects/WARMEM/Old survey data/Survdat_Nye2016.RData"),
                         "2019"    = paste0(res_path,   "NMFS_trawl/Survdat_Nye_allseason.RData"),
                         "2020"    = paste0(res_path,   "NMFS_trawl/Survdat_Nye_Aug 2020.RData"),
                         "2021"    = paste0(res_path,   "NMFS_trawl/survdat_slucey_01152021.RData"),
                         "bigelow" = paste0(res_path,   "NMFS_trawl/survdat_Bigelow_slucey_01152021.RData")
  )
  
  
  
  # If providing a starting point for survdat pass it in:
  if(is.null(survdat) == FALSE){ 
    trawldat <- survdat %>% clean_names()
    
  } else if(is.null(survdat) == TRUE){
    load(survdat_path)
    trawldat <- survdat %>% clean_names()  
    
  }
  
  # remove survdat
  rm(survdat)
  
  
  
  ####  Special Steps for Different SURVDATs  ####
  
  # Use SVSPP to get common names
  if(survdat_source %in% c("2021", "bigelow")){
    
    # Load sppclass codes and common names
    spp_classes <- read_csv(
      paste0(res_path, "NMFS_trawl/spp_keys/sppclass.csv"), 
      col_types = cols()) %>%
      clean_names() %>%
      mutate(common_name     = str_to_lower(common_name),
             scientific_name = str_to_lower(scientific_name)) %>%
      distinct(svspp, comname = common_name, scientific_name)
    
    
    # Add the common names over and format for rest of build
    trawldat <- mutate(trawldat, svspp = str_pad(svspp, 3, "left", "0")) %>% 
      left_join(spp_classes, by = "svspp") %>%
      drop_na(comname) %>%
      
      # Build ID column
      mutate(cruise6 = str_pad(cruise6, 6, "left", "0"),
             station = str_pad(station, 3, "left", "0"),
             stratum = str_pad(stratum, 4, "left", "0"),
             id = str_c(cruise6, station, stratum)) %>%
      
      # Re-order columns
      select(id, year, station, stratum, svvessel, season, 
             lat, lon, depth, surftemp, surfsalin, bottemp, 
             botsalin, svspp, comname, scientific_name, everything()) %>% 
      
      # Rename to match other survdat files
      rename(est_year = year, 
             decdeg_beglat = lat, 
             decdeg_beglon = lon, 
             avgdepth = depth)
    # }
    # 
    # 
    # # seperate tow data
    # if(survdat_source %in% c("2021", "bigelow")){
    trawldat <- mutate(trawldat,
                       est_month = str_sub(est_towdate, 6,7),
                       est_month = as.numeric(est_month),
                       est_day   = str_sub(est_towdate, -2, -1),
                       est_day = as.numeric(est_day),)
  }
  
  
  
  ####__ 1. Column Changes  ####
  trawldat <- trawldat %>% 
    mutate(
      # Text Formatting 
      comname = tolower(comname),
      id      = format(id, scientific = FALSE),
      svspp   = as.character(svspp),
      svspp   = str_pad(svspp, 3, "left", "0"),
      
      # Biomass and abundance NA substitutions
      biomass   = ifelse(is.na(biomass) == TRUE & abundance > 0, 0.0001, biomass),
      abundance = ifelse(is.na(abundance) == TRUE & biomass > 0, 1, abundance),
      
      # Sratum number, excluding leading and trailing codes for inshore/offshore, for matching
      strat_num = str_sub(stratum, 2, 3)) %>%  
    
    # these are redundant, but leaving for now
    mutate(biom_adj  = ifelse(biomass == 0 & abundance > 0, 0.0001, biomass), .after = biomass) %>% 
    mutate(abund_adj = ifelse(abundance == 0 & biomass > 0, 1, abundance), .after = abundance)
  
  
  
  
  
  ####__ 2. Column Selections ####
  
  # currently a light-weight group of columns, 
  # leaves behind CTD and shipboard instrument details.
  # Favors the larger categorical group metadata
  
  # 2016 colnames lack month and day columns
  short_list <- c(
    "id", "cruise6", "station", "est_year", "svvessel", 
    "season", "stratum", "strat_num", "decdeg_beglat", 
    "decdeg_beglon", "avgdepth", "svspp", "comname", 
    "catchsex", "length", "numlen", "abundance", "abund_adj", 
    "biomass", "biom_adj")
  
  # These are the columns we want from the more inclusive pulls (has month and day...)
  long_list <- c(short_list, "est_month", "est_day")
  
  
  # Toggle for different survdat resources
  if(survdat_source %in% c("2016")){
    important_cols <- syms(short_list)
  } else {
    important_cols <- syms(long_list)
  }
  
  
  # select the appropriate columns
  trawldat <- select(trawldat, !!!important_cols)
  
  
  
  
  ####__ 3. Row Filtering  ####
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
      svvessel %in% c("AL", "HB"),
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
  
  
  
  
  
  
  
  
  ####__ 4. Spatial Filtering  ####
  
  
  #### EPU assignment for survdat stations - function from Sean Luceys "RSurvey" repo
  source(here("R/kathy_ss_code/slucey_survdat_functions.R"))
  epu_sf <- ecodata::epu_sf
  epu_sp <- suppressWarnings(as_Spatial(epu_sf))
  
  
  # assign EPU labels
  trawldat <- trawldat %>% 
    rename(CRUISE6 = cruise6,
           STATION = station,
           STRATUM = stratum,
           LAT     = decdeg_beglat,
           LON     = decdeg_beglon) %>% 
    poststrat(survdat = ., stratum = epu_sp, strata.col = "EPU") %>% 
    rename(epu           = newstrata,
           cruise6       = CRUISE6,
           station       = STATION,
           stratum       = STRATUM,
           decdeg_beglat = LAT,
           decdeg_beglon = LON)
  
  
  # Stratum Key for filtering specific areas
  strata_key <- list(
    "Georges Bank"          = as.character(13:23),
    "Gulf of Maine"         = as.character(24:40),
    "Southern New England"  = str_pad(as.character(1:12), width = 2, pad = "0", side = "left"),
    "Mid-Atlantic Bight"    = as.character(61:76))
  
  # Add labels to the data
  trawldat <- trawldat %>%
    mutate(
      survey_area =  case_when(
        strat_num %in% strata_key$`Georges Bank`         ~ "GB",
        strat_num %in% strata_key$`Gulf of Maine`        ~ "GoM",
        strat_num %in% strata_key$`Southern New England` ~ "SNE",
        strat_num %in% strata_key$`Mid-Atlantic Bight`   ~ "MAB",
        TRUE                                             ~ "not found"))
  
  
  # Optional, Use strata_select to pull the strata we want individually
  strata_select <- c(strata_key$`Georges Bank`, strata_key$`Gulf of Maine`,
                     strata_key$`Southern New England`, strata_key$`Mid-Atlantic Bight`)
  
  
  # Filtering with strata_select
  trawldat <- trawldat %>% filter(strat_num %in% strata_select) %>% 
    mutate(stratum = as.character(stratum))
  
  
  
  
  ####__ 5. Stratum Area/Effort Ratios  ####
  # Stratum area ratio is the ratio between the area of the select 
  # stratum to the total area of all stratum sampled that year
  
  # Join to the files containing area of each stratum, epu
  trawldat <- trawldat %>% 
    left_join(stratum_area, by = "stratum") %>% 
    left_join(epu_areas, by = "epu") %>% 
    arrange(trawldat, id)
  
  
  
  
  # Get Total stratum area of all strata 
  # (excludes ones we do not care about via left join)
  total_stratum_areas <- trawldat %>% 
    group_by(est_year) %>% 
    distinct(stratum, .keep_all = T) %>%  
    summarise(tot_s_area =  sum(s_area_km2, na.rm = T),
              .groups = "keep") %>% 
    ungroup()
  
  total_epu_areas <- trawldat %>% 
    group_by(est_year) %>% 
    distinct(epu, .keep_all = T) %>%  
    summarise(tot_epu_area =  sum(epu_area_km2, na.rm = T), 
              .groups = "keep") %>% 
    ungroup
  
  
  # Calculate strata area relative to total area i.e. stratio or stratum weights
  trawldat <- trawldat %>% 
    left_join(total_stratum_areas, by = "est_year") %>% 
    left_join(total_epu_areas, by = "est_year") %>% 
    mutate(st_ratio   = s_area_km2 / tot_s_area,
           epu_ratio  = epu_area_km2 / tot_epu_area) 
  
  
  # Get number of unique tows per stratum
  yr_strat_effort <- trawldat %>% 
    group_by(est_year, stratum) %>% 
    summarise(strat_ntows = n_distinct(id), 
              .groups = "keep") %>% 
    ungroup()
  
  
  yr_epu_effort <-  trawldat %>% 
    group_by(est_year, epu) %>% 
    summarise(epu_ntows = n_distinct(id), 
              .groups = "keep") %>% 
    ungroup()
  
  
  # Add those yearly effort counts back for later
  trawldat <- trawldat %>% 
    left_join(yr_strat_effort, by = c("est_year", "stratum")) %>% 
    left_join(yr_epu_effort, by = c("est_year", "epu"))
  
  
  
  
  
  ####__ 6. Adjusted NumLength  ####
  # Sometimes there are more/less measured than initially tallied*
  
  # Get the abundance value for each sex arrived at by summing across each length
  abundance_check <- trawldat %>%
    group_by(id, comname, catchsex, abundance) %>%
    summarise(
      abund_actual = sum(numlen),               
      n_len_class = n_distinct(length),
      .groups     = "keep") %>% 
    ungroup()
  
  # Get the ratio between the abundance column and the sum of numlen
  conv_factor <- trawldat %>% 
    distinct(id, comname, catchsex, length, abund_adj) %>% 
    inner_join(abundance_check) %>% 
    mutate(convers = abund_adj / abund_actual)
  
  
  
  # Merge back and convert the numlen field
  survdat_processed <- trawldat %>%
    left_join(conv_factor) %>%
    mutate(numlen_adj = numlen * convers, .after = numlen) %>% 
    select(-c(abund_actual, convers))
  
  # remove conversion factor from environment
  rm(abundance_check, conv_factor, strata_key, strata_select, stratum_area, epu_areas, epu_sf, epu_sp)
  
  
  ####__ 7. Distinct Station & Species Length Info   ####
  
  # For each station we need unique combinations of
  # station_id, species, catchsex, length, adjusted_numlen
  # Record of unique station catches: # rows for each species * sex * length
  trawl_lens <- survdat_processed %>% 
    filter(is.na(length) == FALSE,
           is.na(numlen) == FALSE,
           numlen_adj > 0) %>% 
    # Columns that uniquely identify a station and the different catches
    distinct(id, svspp, comname, catchsex, abundance, n_len_class, 
             length,  numlen, numlen_adj, biom_adj)
  
  
  
  # Pull distinct records of the stations themselves and metadata that match
  # these are susceptible to upstream changes
  
  if(survdat_source %in% c("2016")) {
    station_cols <- syms(c(
      # Tow Identification details
      "id", "est_year", "svvessel", "season", 
      # physical location details
      "decdeg_beglat", "decdeg_beglon", "avgdepth", 
      # NMFS/NEFSC Survey Stratum
      "stratum", "strat_num", "s_area_km2", "st_ratio", "strat_ntows", "tot_s_area",
      # Aggregate Regions/EPU's
      "survey_area", "epu", "epu_area_km2", "epu_ratio", "epu_ntows", "tot_epu_area"
    ))
  } else{
    station_cols <- syms(c(
      # Tow Identification details
      "id", "est_year", "est_month", "est_day", "svvessel", "season", 
      # physical location details
      "decdeg_beglat", "decdeg_beglon", "avgdepth", 
      # NMFS/NEFSC Survey Stratum
      "stratum", "strat_num", "s_area_km2", "st_ratio", "strat_ntows", "tot_s_area",
      # Aggregate Regions/EPU's
      "survey_area", "epu", "epu_area_km2", "epu_ratio", "epu_ntows", "tot_epu_area"
    ))
  }
  
  # Pull distinct
  trawl_stations <- survdat_processed %>% 
    select(!!!station_cols) %>% 
    distinct() 
  
  
  # recombine with the distinct station info
  trawl_spectra <- trawl_stations %>% 
    left_join(trawl_lens, by = "id")
  
  
  # Return the dataframe
  # Row for each length class of every species caught
  return(trawl_spectra)
  
}


#### Test it  ####
survdat_clean <- survdat_prep(survdat_source = "2020")



####_______________________________####
####  Add Length weight Info  ####
add_lw_info <- function(survdat_clean, cutoff = FALSE){
  
  
  ####__ 1. Match Species with Growth Coefficients  ####
  
  # This table is a combined table of wigley and fishbase L-W coefficients
  lw_combined <- read_csv(here::here("nmfs_trawl_processing/data/biomass_key_combined.csv"), col_types = cols()) %>% 
    mutate(svspp = str_pad(svspp, 3, "left", "0"))
  
  
  # Do a priority pass with the filter(lw_combined, source == "wigley)
  # merge on comname, season, and catchsex
  wigley_coefficients <- filter(lw_combined, source == "wigley") %>% 
    select(source, season, svspp, comname, scientific_name, spec_class, 
           hare_group, fishery, catchsex, a, b, ln_a)
  
  
  # Do a second pass with the filter(lw_combined, source == "fishbase")
  # merge on common names only
  fishbase_coefficients <- filter(lw_combined, source == "fishbase") %>% 
    select(source, -svspp, comname, scientific_name, spec_class, 
           hare_group, fishery, a, b, ln_a)  
  
  
  # First Pass - Wigley
  # Join just by svspp to account for name changes
  pass_1 <- survdat_clean %>% 
    select(-comname) %>% 
    inner_join(wigley_coefficients)
  
  # Want to pick up stragglers here
  # testing approaches to not lose the rest of these :
  # length(unique(wigley_coefficients$comname))
  # length(unique(pass_1$comname))
  survdat_clean %>% 
    filter(svspp %not in% pass_1$svspp) %>% 
    mutate(comname = ifelse(comname == "windowpane", "windowpane flounder", comname)) %>% 
    inner_join(select(wigley_coefficients, -svspp)) %>% 
    distinct(comname, svspp)
  
  # windowpane is taken care of already, so double counting is occurring
  
  
  
  # Second Pass - Fishbase, for the stragglers if any
  # currently has potential for double matching in event of name changes
  pass_2 <- survdat_clean %>% 
    filter(comname %not in% wigley_coefficients$comname,
           svspp %not in% pass_1$svspp) %>% 
    inner_join(fishbase_coefficients) 
  
  # common names of fishes coming from fishbase
  sort(unique(pass_2$comname))
  
  
  # Join them with bind rows (implicitly drops things that don't have growth coefs)
  trawl_weights <- bind_rows(pass_1, pass_2) %>% 
    arrange(est_year, season) %>% 
    mutate(
      b             = as.numeric(b),
      a             = as.numeric(a),
      a             = ifelse(is.na(a) & !is.na(ln_a), exp(ln_a), a),
      ln_a          = ifelse(is.na(ln_a), log(a), ln_a),  # log of a used if ln_a is isn't already there (some fish just had ln_a reported)
      llen          = log(length),
      ind_log_wt    = ln_a + (b * llen),
      ind_weight_kg = exp(ind_log_wt),                    # weight of an individual in size class
      sum_weight_kg = ind_weight_kg * numlen_adj) %>%     # Individual weight * adjusted numlen
    drop_na(ind_weight_kg) %>% 
    select(-ind_log_wt, -llen)
  
  
  # clean up environment
  rm(pass_1, pass_2, fishbase_coefficients, wigley_coefficients)
  
  
  
  
  ####__ 2. Use Coefficients to Re-calculate Biomass  ####
  
  # calculate total biomass again using weights from key 
  # make a key for the length weight coefficient sources
  survdat_weights <- trawl_weights %>%  
    arrange(est_year, season, comname, length) %>% 
    mutate(lw_group = str_c(comname, season, catchsex)) 
  
  
  
  
  
  
  ####__ 3. Drop comnames that don't align well with BIOMASS
  
  # these species were dropped at 50% mismatch threshold
  # code: 02_survdat_stratification_validation
  
  cutoff_50 <- c(
    "acadian redfish"          , "alewife"                  , "american plaice"         ,
    "american shad"            , "atlantic angel shark"     , "atlantic cod"            ,
    "atlantic croaker"         , "atlantic halibut"         , "atlantic mackerel"       ,
    "atlantic sharpnose shark" , "atlantic spadefish"       , "atlantic sturgeon"       ,
    "atlantic thread herring"  , "atlantic torpedo"         , "atlantic wolffish"       ,
    "blackbelly rosefish"      , "blueback herring"         , "bluefish"                ,
    "bluntnose stingray"       , "buckler dory"             , "bullnose ray"            ,
    "butterfish"               , "chain dogfish"            , "cownose ray"             ,
    "cunner"                   , "cusk"                     , "fawn cusk-eel"           ,
    "greater amberjack"        , "haddock"                  , "longhorn sculpin"        ,
    "northern kingfish"        , "ocean pout"               , "offshore hake"           ,
    "pollock"                  , "rosette skate"            , "roughtail stingray"      ,
    "round herring"            , "sand tiger"               , "sandbar shark"           ,
    "sea raven"                , "smooth butterfly ray"     , "smooth dogfish"          ,
    "southern kingfish"        , "spanish mackerel"         , "spanish sardine"         ,
    "spiny butterfly ray"      , "spiny dogfish"            , "spot"                    ,
    "striped bass"             , "tautog"                   , "thorny skate"            ,
    "weakfish"                 , "white hake"               , "windowpane flounder"     ,
    "winter flounder"          , "winter skate"             , "witch flounder"          ,
    "yellowtail flounder"
  )
  
  if(cutoff == TRUE){
    survdat_weights <- survdat_weights %>% filter(comname %in% cutoff_50)
  }
  
  
  
  
  return(survdat_weights)
}


####  Test it  ####
survdat_lw <- add_lw_info(survdat_clean = survdat_clean, cutoff = F)


####_______________________________####
####  Add Area Stratification  ####
add_area_stratification <- function(survdat_weights, include_epu = F){ 
  
  # Constants:
  # average area covered by an albatross standard tow in km2
  alb_tow_km2 <- 0.0384 
  
  # catchability coefficient, ideally should shift for different species guilds
  q <- 1                
  
  
  # Derived Estimates:
  survdat_weights <- survdat_weights  %>% 
    mutate(
      # Abundance per tow
      # abundance / ntows for the year within that strata/epu
      abund_tow_s   = numlen_adj / strat_ntows,    
      
      # Biomass is repeated across length classes at each station by species
      # the number of length classes is tallied where the conversion factor is done
      biom_per_lclass = (biom_adj / n_len_class),
      
      # Mean biomass/tow for the "biomass" column
      biom_tow_s      = biom_per_lclass / strat_ntows,
      
      # Stratified mean abundance, weighted by the stratum areas
      wt_abund_s     = abund_tow_s * st_ratio, 
      
      # Stratified mean BIOMASS
      wt_biom_s   = biom_tow_s * st_ratio,
      
      # convert from catch rate by area swept to total catch for entire stratum
      # So catch/tow times the total area, divided by how many tows would cover that area
      expanded_abund_s   = round((wt_abund_s * tot_s_area / alb_tow_km2) / q),
      
      # Total BIOMASS from the weighted biomass
      expanded_biom_s   = round((wt_biom_s * tot_s_area / alb_tow_km2) / q), 
      
      # Total lw-biomass from Projected abundances: Biomass = abundance * lw_weight
      expanded_lwbio_s    = sum_weight_kg * expanded_abund_s
      
    ) 
  
  ####  Optional - weighted by EPU areas
  if(include_epu == TRUE){
    survdat_weights <- survdat_weights  %>% 
      mutate(
        # Abundance per tow
        abund_tow_epu = numlen_adj / epu_ntows,       
        # Mean biomass/tow for the BIOMASS column
        biom_tow_epu    = biom_per_lclass / epu_ntows,
        # Stratified mean abundances, weighted by the stratum areas
        wt_abund_epu   = abund_tow_epu * epu_ratio,
        # Stratified mean BIOMASS
        wt_biom_epu = biom_tow_epu * epu_ratio,
        # Total catch for entire stratum
        expanded_abund_epu = round((wt_abund_epu * tot_epu_area/ alb_tow_km2) / q),
        # Total Biomass from the weighted biomass
        expanded_biom_epu = round((wt_biom_epu * tot_epu_area/ alb_tow_km2) / q),
        # LW Biomass from Expanded abundances: Biomass = abundance * lw_weight
        expanded_lwbio_epu  = sum_weight_kg * expanded_abund_epu)} 
  
  
  # Remove instances where there were fish that were measured but not weighed
  survdat_weights <- survdat_weights %>% 
    filter(expanded_lwbio_s != -Inf)
  
  return(survdat_weights)
}


#### Test it  ####
survdat_stratified <- add_area_stratification(survdat_weights = survdat_lw, include_epu = F)
