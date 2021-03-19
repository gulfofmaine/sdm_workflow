

#' @title Make VAST dataset
#' 
#' @description 
#'
#' @param tidy_mod_data
#' @param year_min
#' @param year_max
#' @param out_dir
#' 
#' @return  
#' 
#' @export

process_vast_seasonal_data<- function(tidy_mod_data, year_min, year_max, out_dir){
  
  # For debugging
  if(FALSE){
    tidy_mod_data = readRDS(here::here("scratch/aja/targets_flow/data/combined/tidy_mod_data.rds"))
  }
  
  # Some work on the time span and seasons
  data_temp<- tidy_mod_data %>%
    filter(., EST_YEAR >= year_min & EST_YEAR <= year_max) %>%
    filter(., SEASON %in% c("SPRING", "FALL")) %>%
    mutate(., "VAST_SEASON" = case_when(
      SEASON == "SPRING" & SURVEY == "DFO" ~ "DFO",
      SEASON == "SPRING" & SURVEY == "NMFS" ~ "SPRING",
      SEASON == "FALL" & SURVEY == "NMFS" ~ "FALL"))
  data_temp$VAST_SEASON<- ifelse(data_temp$SEASON == "SPRING" & data_temp$SURVEY == "DFO", "DFO", data_temp$SEASON)
  
  # Set of years and seasons. The DFO spring survey usually occurs before the NOAA NEFSC spring survey, so ordering accordingly.
  year_set<- sort(unique(data_temp$EST_YEAR))
  season_set<- c("DFO", "SPRING", "FALL")
  
  # Create a grid with all unique combinations of seasons and years and then combine these into one "year_season" variable
  yearseason_grid<- expand.grid("SEASON" = season_set, "EST_YEAR" = year_set)
  yearseason_levels<- apply(yearseason_grid[, 2:1], MARGIN = 1, FUN = paste, collapse = "_")
  yearseason_labels<- round(yearseason_grid$EST_YEAR + (as.numeric(factor(yearseason_grid$SEASON, levels = season_set))-1)/length(season_set), digits = 1)
  
  # Similar process, but for the observations
  yearseason_i<- apply(data_temp[, c("EST_YEAR", "SEASON")], MARGIN = 1, FUN = paste, collapse = "_")
  yearseason_i<- factor(yearseason_i, levels = yearseason_levels)
  
  # Add the year_season factor column to our sampling_data data set
  data_temp$YEAR_SEASON<- yearseason_i
  data_temp$VAST_SEASON = factor(data_temp$VAST_SEASON, levels = season_set)
  
  # Make dummy data for all year_seasons to estimate gaps in sampling if needed
  dummy_data<- data.frame(ID = sample(data_temp$ID, size = 1), Year = yearseason_grid[,'EST_YEAR'], Season = yearseason_grid[,'SEASON'], Year_Season = yearseason_levels, Lat = mean(data_temp$Lat, na.rm = TRUE), Lon = mean(data_temp$Lon, na.rm = TRUE), Biomass = 0, Swept = 0.05, PredTF = TRUE)
  
  # Combine with original dataset
  vast_data_out<- rbind(cbind(data_temp, PredTF = FALSE), dummy_data)

  # Save and return it
  saveRDS(vast_cov_out, file = paste(out_dir, "vast_cov.rds", sep = "/"))
  return(vast_cov_out)
}

#' @title Make VAST sample dataset
#' 
#' @description A short function that processes a sample dataset to work with VAST function and, if needed, extraction functions to get other covariate values before being used to fit VAST model
#'
#' @param all_sample_dir
#' @param year_min
#' @param year_max
#' @param out_dir
#' 
#' @return Combined sample data frame 
#' 
#' @export

vast_make_sample_data<- function(all_sample_dir, year_min, year_max, nmfs_species_code, out_dir){
  
  # For debugging
  if(FALSE){
    all_sample_dir = here::here("scratch/aja/targets_flow/data/combined")
    year_min = 1985
    year_max = 2015
    nmfs_species_code = 105
    out_dir = here::here("scratch/aja/targets_flow/data/dfo/combined")
  }
  
  # Bring in the data
  all_samp<- readRDS(here::here(all_sample_dir, "all_sample.rds"))
  
  # Some work on the time span and seasons
  samp_sub<- all_samp %>%
    filter(., EST_YEAR >= year_min & EST_YEAR <= year_max) %>%
    filter(., SEASON %in% c("SPRING", "FALL"))
  samp_sub$SEASON<- ifelse(samp_sub$SEASON == "SPRING" & samp_sub$SURVEY == "DFO", "DFO", samp_sub$SEASON)
  
  # Set of years and seasons. The DFO spring survey usually occurs before the NOAA NEFSC spring survey, so ordering accordingly.
  year_set<- sort(unique(samp_sub$EST_YEAR))
  season_set<- c("DFO", "SPRING", "FALL")
  
  # Create a grid with all unique combinations of seasons and years and then combine these into one "year_season" variable
  yearseason_grid<- expand.grid("SEASON" = season_set, "EST_YEAR" = year_set)
  yearseason_levels<- apply(yearseason_grid[, 2:1], MARGIN = 1, FUN = paste, collapse = "_")
  yearseason_labels<- round(yearseason_grid$EST_YEAR + (as.numeric(factor(yearseason_grid$SEASON, levels = season_set))-1)/length(season_set), digits = 1)
  
  # Similar process, but for the observations
  yearseason_i<- apply(samp_sub[, c("EST_YEAR", "SEASON")], MARGIN = 1, FUN = paste, collapse = "_")
  yearseason_i<- factor(yearseason_i, levels = yearseason_levels)
  
  # Add the year_season factor column to our sampling_data data set
  samp_sub$YEAR_SEASON<- yearseason_i
  samp_sub$SEASON = factor(samp_sub$SEASON, levels = season_set)
  
  # Subset and renaming...
  samp_sub2<- samp_sub %>%
    filter(., NMFS_SVSPP == nmfs_species_code) %>%
    select(., ID, EST_YEAR, SEASON, YEAR_SEASON, DECDEG_BEGLON, DECDEG_BEGLAT, BIOMASS) %>%
    rename(., c(Year = EST_YEAR, Season = SEASON, Year_Season = YEAR_SEASON, Lat = DECDEG_BEGLAT, Lon = DECDEG_BEGLON, Biomass = BIOMASS))
  samp_sub2$Swept = rep(0.05, nrow(samp_sub2))
  
  # Make dummy data for all year_seasons to estimate gaps in sampling if needed
  dummy_data<- data.frame(ID = sample(samp_sub2$ID, size = 1), Year = yearseason_grid[,'EST_YEAR'], Season = yearseason_grid[,'SEASON'], Year_Season = yearseason_levels, Lat = mean(samp_sub2$Lat, na.rm = TRUE), Lon = mean(samp_sub2$Lon, na.rm = TRUE), Biomass = 0, Swept = 0.05, PredTF = TRUE)
  
  # Combine with sampling data
  vast_samp_out<- rbind(cbind(samp_sub2, PredTF = FALSE), dummy_data)
  
  # Final step
  vast_samp_out<- vast_samp_out %>%
    mutate(., Year = as.numeric(Year_Season)-1) %>%
    select(., Year, Lat, Lon, Biomass, Swept, PredTF)
  
  # Save and return it
  saveRDS(vast_samp_out, file = paste(out_dir, "vast_sample.rds", sep = "/"))
  return(vast_samp_out)
}

#' @title Read in sample dataset 
#' 
#' @description A short function to read in sample data given a file path to csv
#'
#' @param samp_dat_path File path to sample data csv
#' 
#' @return sample data frame 
#' 
#' @export

read_samp_dat_csv<- function(samp_dat_path){
  
  # For debugging
  if(FALSE){
    samp_dat_path = here::here("scratch/aja/targets_flow/data", "sample_data.csv")
  }
  
  # Read in sample data from file path
  samp_dat<- read.csv(samp_dat_path)
  
  # Factor work
  samp_dat$Season<- factor(samp_dat$Season, levels = c("DFO", "SPRING", "FALL"))
 
  # Return it
  return(samp_dat)
}

#' @title Read in covariate dataset 
#' 
#' @description A short function to read in covariate data given a file path to csv
#'
#' @param cov_dat_path File path to sample data csv
#' 
#' @return covariate data frame
#' 
#' @export

read_cov_dat_csv<- function(cov_dat_path){
  
  # For debugging
  if(FALSE){
    cov_dat_path = here::here("scratch/aja/targets_flow/data", "covariate_data.csv")
  }
  
  # Read in covariate data from file path
  cov_dat<- read.csv(cov_dat_path)
  
  # Factor work
  cov_dat$Season<- factor(cov_dat$Season, levels = c("DFO", "SPRING", "FALL"))
  cov_dat$Year_Cov<- factor(cov_dat$Year_Cov, levels = unique(cov_dat$Year_Cov))
 
  # Return it
  return(cov_dat)
}

#' @title Read in shapefile
#' 
#' @description A short function to read in a shapefile given a file path
#'
#' @param file_path File path to geospatial vector polygon file with .shp extension, specifying the location and shape of the area of interest.
#' @param factor_vars Names of factor columns that should be checked and converted if necessary

#' 
#' @return SF poylgon 
#' 
#' @export

read_polyshape<- function(polyshape_path){
  
  # For debugging
  if(FALSE){
    polyshape_path = "~/Box/RES_Data/Shapefiles/NELME_regions/NELME_sf.shp"
  }
  
  # Read in polygon shapefile from file_path
  shapefile<- st_read(polyshape_path)
  
  # Return it
  return(shapefile)
}

####
#' @title Make VAST extrapolation grid settings from a shapefile
#' 
#' @description Create a list of with information defining the extrapolation grid and used by subsequent VAST functions, leveraging code here: https://github.com/James-Thorson-NOAA/VAST/wiki/Creating-an-extrapolation-grid. 
#'
#' @param shapefile A geospatial vector sf polygon file, specifying the location and shape of the area of interest.
#' @param cell_size The size of grid in meters (since working in UTM). This will control the resolution of the extrapolation grid.
#'
#' @return Tagged list containing extrapolation grid settings needed to fit a VAST model of species occurrence.
#' 
#' @export

vast_make_extrap_grid<- function(shapefile, cell_size){
  
  # For debugging
  if(FALSE){
    shapefile = shapefile
    cell_size = 25000
  }
  
  # Transform crs of shapefile to common WGS84 lon/lat format.
  shape_wgs84<- st_transform(shapefile, crs = "+proj=longlat +lat_0=90 +lon_0=180 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0 ")
  
  # Get UTM zone
  lon<- sum(st_bbox(shape_wgs84)[c(1,3)])/2
  utm_zone<- floor((lon + 180)/6)+1
  
  # Transform to the UTM zone
  crs_utm<- st_crs(paste0("+proj=utm +zone=",utm_zone," +ellps=WGS84 +datum=WGS84 +units=m +no_defs "))
  shape_utm<- st_transform(shape_wgs84, crs = crs_utm)
  
  # Make extrapolation grid with sf
  region_grid<- st_as_sf(st_make_grid(shape_utm, cellsize = cell_size, what = "centers"), crs = crs_utm) 
  
  # Now get only the points that fall within the shape polygon
  points_keep<- data.frame("pt_row" = seq(from = 1, to = nrow(region_grid), by = 1), "in_out" = st_intersects(region_grid, shape_utm, sparse = FALSE))            
  region_grid<- region_grid %>%
    mutate(., "in_poly" = st_intersects(region_grid, shape_utm, sparse = FALSE)) %>%
    filter(., in_poly == TRUE)
  
  # Finally convert back to WGS84 lon/lat, as that is what VAST expects.
  extrap_grid<- region_grid %>%
    st_transform(., crs = "+proj=longlat +lat_0=90 +lon_0=180 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0 ") %>%
    mutate(., "Lon" = st_coordinates(.)[,1],
           "Lat" = st_coordinates(.)[,2]) %>%
    st_drop_geometry() %>%
    dplyr::select(., Lon, Lat) %>%
    mutate(., Area_km2=((cell_size/1000^2)))
    
  # Return it
  return(extrap_grid)
}

####
#' @title Make VAST model settings
#' 
#' @description Create a list of model settings needed to fit a VAST model for species occurrence, largely copied from VAST::make_settings
#'
#' @param extrap_grid User created extrapolation grid from vast_make_extrap_grid.
#' @param FieldConfig A vector defining the number of spatial (Omega) and spatio-temporal (Epsilon) factors to include in the model for each of the linear predictors. For each factor, possible values range from 0 (which effectively turns off a given factor), to the number of categories being modeled. If FieldConfig < number of categories, VAST estimates common factors and then loading matrices. 
#' @param RhoConfig A vector defining the temporal structure of intercepts (Beta) and spatio-temporal (Epsilon) variation for each of the linear predictors. See `VAST::make_data` for options.
#' @param bias.correct Logical boolean determining if Epsilon bias-correction should be done. 
#'
#' @return Tagged list containing settings needed to fit a VAST model of species occurrence.
#' 
#' @export

vast_make_settings <- function(extrap_grid, FieldConfig, RhoConfig, bias.correct){
  
  # For debugging
  if(FALSE){
    extrap_grid = extrap_grid
    FieldConfig = c("Omega1" = 1, "Epsilon1" = 1, "Omega2" = 1, "Epsilon2" = 1)
    RhoConfig = c("Beta1" = 3, "Beta2" = 3, "Epsilon1" = 4, "Epsilon2" = 4)
    bias.correct = FALSE
  }
  
  # Get number of vertices in the mesh, which is based on our input_grid
  n_x_use<- nrow(extrap_grid)
  
  # Run FishStatsUtils::make_settings
  settings_out<- make_settings(n_x = n_x_use, Region = "User", purpose = "index2", FieldConfig = FieldConfig, RhoConfig = RhoConfig, ObsModel = c(1, 1), bias.correct = bias.correct, knot_method = "grid")
  
  # Return it
  return(settings_out)
}

####
#' @title Make VAST covariate effect objects
#' 
#' @description Create covariate effects for both linear predictors
#'
#' @param X1_coveff_vec A vector specifying the covariate effects for first linear predictor. 
#' @param X2_coveff_vec A vector specifying the covariate effects for second linear predictor. 
#'
#' @return A list with covariate effects for the first linear predictor (first list slot) and second linear predictor (second list slot). 
#'
#' @export

vast_make_coveff<- function(X1_coveff_vec, X2_coveff_vec){
  
  # For debugging
  if(FALSE){
   X1_coveff_vec = c(2, 3, 3, 2, rep(3, 32))
   X2_coveff_vec = c(2, 3, 3, 2, rep(3, 32))
  }
  
  # Combine into a list and name it 
  coveff_out<- list("X1config_cp" = matrix(X1_coveff_vec, nrow = 1), "X2config_cp" = matrix(X2_coveff_vec, nrow = 1))
  
  # Return it
  return(coveff_out)
}

####
#' @title Build VAST SDM
#' 
#' @description Build VAST species distribution model, without running it. This can be helpful to check settings before running `vast_fit_sdm`. Additionally, it can be helpful for making subsequent modifications, particularly to mapping.
#'
#' @param settings A tagged list with the settings for the model, created with `vast_make_settings`.
#' @param input_grid An extrapolation grid, created with `vast_make_extrap_grid`.
#' @param samp_dat A data frame with the biomass sample data for each species at each tow.
#' @param cov_dat A data frame with the covariate data for each tow.
#' @param X1_formula A formula for the first linear predictor.
#' @param X2_formula A formula for the second linear predictor.
#' @param X_contrasts A tagged list specifying the contrasts to use for factor covariates in the model.
#' @param Xconfig_list A tagged list specifying the covariate effects for first and second linear predictors.
#'
#' @return A VAST `fit_model` object, with the inputs and built TMB object components.
#'
#' @export

vast_build_sdm <- function(settings, extrap_grid, samp_dat, cov_dat, X1_formula, X2_formula, X_contrasts, Xconfig_list){
  
  # For debugging
  if(FALSE){
    settings = settings_out
    extrap_grid = extrap_grid
    samp_dat = samp_dat
    cov_dat = cov_dat
    X1_formula = ~ Season + Year_Cov
    X2_formula = ~ Season + Year_Cov
    X_contrasts = list(Season = contrasts(cov_dat$Season, contrasts = FALSE), Year_Cov = contrasts(cov_dat$Year_Cov, contrasts = FALSE))
    Xconfig_list = coveff_out
  }
  
  # Check names
  samp_dat_names<- c("Lat", "Lon", "Year", "Biomass", "Swept", "PredTF")
  if(!(all(samp_dat_names %in% names(samp_dat)))){
    stop(paste("Check names in sample data. Must include:", paste0(samp_dat_names, collapse = ","), sep = " "))
  }
  
  cov_dat_names1<- unlist(str_extract_all(X1_formula, boundary("word"))[[2]])
  cov_dat_names2<- unlist(str_extract_all(X2_formula, boundary("word"))[[2]])
  cov_dat_names_all<- unique(c(cov_dat_names1, cov_dat_names2))
  if(!(all(cov_dat_names_all %in% names(cov_dat)))){
    stop(paste("Check names in covariate data. Must include", paste0(cov_dat_names_all, collapse = ","), sep = " "))
  }
  
  if(!(all(c("X1config_cp", "X2config_cp") %in% names(Xconfig_list)))){
    stop(paste("Check names of Xconfig_list. Must be", paste0(c("X1config_cp", "X2config_cp"), collapse = ","), sep = ""))
  }
  
  # Run VAST::fit_model with correct info and settings
  vast_build_out<- fit_model("settings" = settings, input_grid = as.matrix(extrap_grid, ncol = 3), "Lat_i" = samp_dat[, 'Lat'], "Lon_i" = samp_dat[, 'Lon'], "t_i" = samp_dat[, 'Year'], "c_i" = rep(0, nrow(samp_dat)), "b_i" = samp_dat[, 'Biomass'], "a_i" = samp_dat[, 'Swept'], "PredTF_i" = samp_dat[, 'PredTF'], "X1config_cp" = Xconfig_list[['X1config_cp']], "X2config_cp" = Xconfig_list[['X2config_cp']], "covariate_data" = cov_dat, "X1_formula" = X1_formula, "X2_formula" = X2_formula, X_contrasts = X_contrasts, "newtonsteps" = 1, "getsd" = TRUE, "getReportCovariance" = TRUE, "run_model" = FALSE, "test_fit" = FALSE,  "Use_REML" = FALSE, "getJointPrecision" = FALSE)
  
  # Return it
  return(vast_build_out)
}

####
#' @title Adjust VAST SDM
#' 
#' @description Make adjustments to VAST SDM and the model returned in `vast_build_sdm`. This can either be the exact same as the one built using `vast_build_sdm`, or it can update that model with adjustments provided in a tagged list. 
#'
#' @param vast_build A VAST `fit_model` object.
#' @param adjustments Either NULL (default) or a tagged list identifying adjustments that should be made to the vast_build `fit_model` object. If NULL, the identical model defined by the `vast_build` is run and fitted.
#'
#' @return A VAST fit_model object, with the inputs and built TMB object components.
#' 
#' @export

vast_make_adjustments <- function(vast_build, adjustments = NULL){
  
  # For debugging
  if(FALSE){
    vast_build = vast_build_out
    adjustments = list("FieldConfig" = c("Omega1" = 1, "Epsilon1" = 0, "Omega2" = 1, "Epsilon2" = 0))
  }
  
  # If no adjustments are needed, just need to pull information from vast_build and then set "run_model" to TRUE
  if(is.null(adjustments)){
    vast_build_adjust_out<- fit_model("settings" = vast_build$settings, input_grid = vast_build$input_args$data_args_input$input_grid, "Lat_i" = vast_build$data_frame[, 'Lat_i'], "Lon_i" = vast_build$data_frame[, 'Lon_i'], "t_i" = vast_build$data_frame[, 't_i'], "c_iz" = vast_build$data_frame[, 'c_iz'], "b_i" = vast_build$data_frame[, 'b_i'], "a_i" = vast_build$data_frame[, 'a_i'], "PredTF_i" = vast_build$data_list[['PredTF_i']], "X1config_cp" = vast_build$input_args$data_args_input[['X1config_cp']], "X2config_cp" = vast_build$input_args$data_args_input[['X2config_cp']], "covariate_data" = vast_build$input_args$data_args_input$covariate_dat, "X1_formula" = vast_build$input_args$data_args_input$X1_formula, "X2_formula" = vast_build$input_args$data_args_input$X2_formula, X_contrasts = vast_build$input_args$data_args_input$X_contrasts, "newtonsteps" = 1, "getsd" = TRUE, "getReportCovariance" = TRUE, "run_model" = FALSE, "test_fit" = FALSE,  "Use_REML" = FALSE, "getJointPrecision" = FALSE)
  }
  
  # If there are adjustments, need to make those and then re run model. 
  if(!is.null(adjustments)){
    # Check names -- trying to think of what the possible adjustment flags would be in the named list
    adjust_names<- c("FieldConfig", "RhoConfig", "X1_formula", "X2_formula", "X1config_cp", "X2config_cp", "X_contrasts", "log_sigmaXi1_cp", "log_sigmaXi2_cp")
    
    if(!(all(names(adjustments) %in% adjust_names))){
      stop(paste("Check names in adjustment list. Must be one of", paste0(adjust_names, collapse = ","), sep = " "))
    }
    
    # First options are going to be in the settings bit..
    if(any(names(adjustments) %in% c("FieldConfig", "RhoConfig"))){
      # Get just the settings adjustments
      settings_adjusts<- names(adjustments)[which(names(adjustments) %in% names(vast_build$settings))]
      
      for(i in seq_along(settings_adjusts)){
        setting_adjust_i<- settings_adjusts[i]
        vast_build$settings[[{{setting_adjust_i}}]]<- adjustments[[{{setting_adjust_i}}]]
      }
    }
    
    # A lot of stuff is going to be in the `vast_build$input_args$data_args_input` object
    if(any(names(adjustments) %in% names(vast_build$input_args$data_args_input))){
      
      # Get just the data args adjustments
      data_adjusts<- names(adjustments)[which(names(adjustments) %in% names(vast_build$input_args$data_args_input))]
      
      for(i in seq_along(data_adjusts)){
        data_adjust_i<- data_adjusts[i]
        vast_build$input_args$data_args_input[[{{data_adjust_i}}]]<- adjustments[[{{data_adjust_i}}]]
      }
    }
    
    # Only other adjustment (for now) is Map. 
    if(any(names(adjustments) %in% c("log_sigmaXi1_cp", "log_sigmaXi2_cp"))){
      # Get the original, which we can then edit...
      map_adjust_out<- vast_build$tmb_list$Map
      
      # Get just the map adjustment names
      map_adjusts<- names(adjustments)[which(names(adjustments) %in% names(vast_build$tmb_list$Map))]
      
      # Loop over them
      for(i in seq_along(map_adjusts)){
        map_adjust_i<- map_adjusts[i]
        map_adjust_out[[{{map_adjust_i}}]]<- adjustments[[{{map_adjust_i}}]]
      }
    }
    
    # Now, re-build and fit model. This is slightly different if we have changed map or not...
    if(any(names(adjustments) %in% c("log_sigmaXi1_cp", "log_sigmaXi2_cp"))){
      # Adding Map argument
      vast_build_adjust_out<- fit_model("settings" = vast_build$settings, input_grid = vast_build$input_args$data_args_input$input_grid, "Lat_i" = vast_build$data_frame[, 'Lat_i'], "Lon_i" = vast_build$data_frame[, 'Lon_i'], "t_i" = vast_build$data_frame[, 't_i'], "c_iz" = vast_build$data_frame[, 'c_iz'], "b_i" = vast_build$data_frame[, 'b_i'], "a_i" = vast_build$data_frame[, 'a_i'], "PredTF_i" = vast_build$data_list[['PredTF_i']], "X1config_cp" = vast_build$input_args$data_args_input[['X1config_cp']], "X2config_cp" = vast_build$input_args$data_args_input[['X2config_cp']], "covariate_data" = vast_build$input_args$data_args_input$covariate_dat, "X1_formula" = vast_build$input_args$data_args_input$X1_formula, "X2_formula" = vast_build$input_args$data_args_input$X2_formula, X_contrasts = vast_build$input_args$data_args_input$X_contrasts, Map = map_adjust_out, "newtonsteps" = 1, "getsd" = TRUE, "getReportCovariance" = TRUE, "run_model" = FALSE, "test_fit" = FALSE,  "Use_REML" = FALSE, "getJointPrecision" = FALSE)
    } else {
      # No need for Map argument, just build and fit
      vast_build_adjust_out<- fit_model("settings" = vast_build$settings, input_grid = vast_build$input_args$data_args_input$input_grid, "Lat_i" = vast_build$data_frame[, 'Lat_i'], "Lon_i" = vast_build$data_frame[, 'Lon_i'], "t_i" = vast_build$data_frame[, 't_i'], "c_iz" = vast_build$data_frame[, 'c_iz'], "b_i" = vast_build$data_frame[, 'b_i'], "a_i" = vast_build$data_frame[, 'a_i'], "PredTF_i" = vast_build$data_list[['PredTF_i']], "X1config_cp" = vast_build$input_args$data_args_input[['X1config_cp']], "X2config_cp" = vast_build$input_args$data_args_input[['X2config_cp']], "covariate_data" = vast_build$input_args$data_args_input$covariate_dat, "X1_formula" = vast_build$input_args$data_args_input$X1_formula, "X2_formula" = vast_build$input_args$data_args_input$X2_formula, X_contrasts = vast_build$input_args$data_args_input$X_contrasts, "newtonsteps" = 1, "getsd" = TRUE, "getReportCovariance" = TRUE, "run_model" = FALSE, "test_fit" = FALSE,  "Use_REML" = FALSE, "getJointPrecision" = FALSE)
    }
  }
  
  # Return it
  return(vast_build_adjust_out)
}

#' @title Fit VAST SDM
#' 
#' @description Fit VAST SDM
#'
#' @param vast_build A VAST `fit_model` object.
#'
#' @return A VAST fit_model object, with the inputs and and outputs, including parameter estimates, extrapolation gid info, spatial list info, data info, and TMB info.
#'
#' @export

vast_fit_sdm <- function(vast_build_adjust){
  
  # For debugging
  if(FALSE){
   vast_build_adjust = vast_build_adjust_out
  }
  
  # Build and fit model
  vast_fit_out<- fit_model("settings" = vast_build_adjust$settings, input_grid = vast_build_adjust$input_args$data_args_input$input_grid, "Lat_i" = vast_build_adjust$data_frame[, 'Lat_i'], "Lon_i" = vast_build_adjust$data_frame[, 'Lon_i'], "t_i" = vast_build_adjust$data_frame[, 't_i'], "c_iz" = vast_build_adjust$data_frame[, 'c_iz'], "b_i" = vast_build_adjust$data_frame[, 'b_i'], "a_i" = vast_build_adjust$data_frame[, 'a_i'], "PredTF_i" = vast_build_adjust$data_list[['PredTF_i']], "X1config_cp" = vast_build_adjust$input_args$data_args_input[['X1config_cp']], "X2config_cp" = vast_build_adjust$input_args$data_args_input[['X2config_cp']], "covariate_data" = vast_build_adjust$input_args$data_args_input$covariate_dat, "X1_formula" = vast_build_adjust$input_args$data_args_input$X1_formula, "X2_formula" = vast_build_adjust$input_args$data_args_input$X2_formula, "X_contrasts" = vast_build_adjust$input_args$data_args_input$X_contrasts, "Map" = vast_build_adjust$tmb_list$Map, "newtonsteps" = 1, "getsd" = TRUE, "getReportCovariance" = TRUE, "run_model" = TRUE, "test_fit" = FALSE,  "Use_REML" = FALSE, "getJointPrecision" = FALSE)
  
  # Return it
  return(vast_fit_out)
}
