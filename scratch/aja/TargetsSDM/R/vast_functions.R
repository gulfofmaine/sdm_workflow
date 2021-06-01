####  Common Resources  ####

pred_template_load<- function(pred_template_dir){
  if(FALSE){
    tar_load(pred_template_dir)
  }
  
  # Load the raster template gird
  pred_template_rast<- raster(paste(pred_template_dir, "mod_pred_template.grd", sep = "/"))
  
  # Convert it to a data frame
  pred_template_df<- as.data.frame(pred_template_rast, xy = TRUE) %>%
    drop_na() %>%
    dplyr::select(., x, y)
  names(pred_template_df)<- c("longitude", "latitude")
  
  # Return it
  return(pred_template_df)
}

high_res_load <- function(high_res_dir) {
  high_res<- raster(paste(high_res_dir, "HighResTemplate.grd", sep = "/"))
  return(high_res)
}

####  Functions  ####
####

#' @title Make VAST prediction dataframe 
#'
#' @description This function creates a dataframe of prediction covariates to combine with the other VAST data
#'
#' @param predict_covariates_processed_dir = The directory holding processed covariate raster stacks
#' @param mask = Shapefile mask
#' @param summarize = Currently, either "annual" or "seasonal" to indicate whether the each dynamic raster stack should be summarized to an annual or seasonal time scale
#' @param ensemble_stat = Either the climate model ensemble statistic to use when working with climate model projections, or NULL. This is only used in naming the output file
#' @param fit_year_min
#' @param fit_year_max
#' @param pred_years
#' @param out_dir = Directory to save the prediction dataframe
#' 
#' @return A dataframe with prediction information. This file is also saved in out_dir. 
#'
#' @export

make_vast_predict_df<- function(predict_covariates_processed_dir, extra_covariates_stack, mask, summarize, ensemble_stat, fit_seasons, fit_year_min, fit_year_max, pred_years, out_dir){
  
  # For debugging
  if(FALSE){
    tar_load(predict_covariates_processed_dir)
    tar_load(static_covariates_stack)
    extra_covariates_stack = static_covariates_stack
    tar_load(shapefile)
    mask = shapefile
    summarize<- "seasonal"
    ensemble_stat<- "mean"
    fit_year_min = fit_year_min
    fit_year_max = fit_year_max
    pred_years = pred_years
    out_dir = here::here("scratch/aja/TargetsSDM/data/predict")
  }
  
  ####
  ## Need to figure out what to do about depth here!!!
  
  # Get raster stack covariate files
  rast_files_load<- list.files(predict_covariates_processed_dir, pattern = paste0(summarize, "_", ensemble_stat, ".grd"), full.names = TRUE)
  
  # Get variable names
  cov_names_full<- list.files(predict_covariates_processed_dir, pattern = paste0(summarize, "_", ensemble_stat, ".grd"), full.names = FALSE)
  predict_covs_names<- gsub(paste("_", ensemble_stat, ".grd", sep = ""), "", gsub("predict_stack_", "", cov_names_full))
  
  # Looping through prediction stack time steps 
  for(i in 1:nlayers(raster::stack(rast_files_load[1]))){
    # Get the time index
    time_ind<- i
    
    # Load corresponding raster layers matching the time index 
    pred_covs_stack_temp<- rotate(raster::stack(raster::stack(rast_files_load[1])[[time_ind]], raster::stack(rast_files_load[2])[[time_ind]], raster::stack(rast_files_load[3])[[time_ind]], raster::stack(rast_files_load[4])[[time_ind]]))
    
    # Mask out values outside area of interest
    pred_covs_stack_temp<- raster::mask(pred_covs_stack_temp, mask = mask)
    
    # Some processing to keep observations within our area of interest and get things in a "tidy-er" prediction dataframe
    time_name<- sub('.[^.]*$', '', names(pred_covs_stack_temp))
    names(pred_covs_stack_temp)<- paste(time_name, predict_covs_names, sep = "_")
    pred_covs_df_temp<- as.data.frame(pred_covs_stack_temp, xy = TRUE) %>%
      drop_na()
    colnames(pred_covs_df_temp)[2:ncol(pred_covs_df_temp)]<- gsub("X", "", gsub("[.]", "_", colnames(pred_covs_df_temp)[2:ncol(pred_covs_df_temp)]))
    colnames(pred_covs_df_temp)[1:2]<- c("DECDEG_BEGLON", "DECDEG_BEGLAT")
    
    pred_covs_df_out_temp<- pred_covs_df_temp %>%
      pivot_longer(., -c(DECDEG_BEGLON, DECDEG_BEGLAT), names_to = c("variable"), values_to = "value") %>%
      separate(., variable, into = c("EST_YEAR", "SEASON", "variable"), sep = "_", extra = "merge") %>%
      pivot_wider(., names_from = variable, values_from = value)
    
    # Adding in some other columns we will want to match up easily with 'vast_data_out'
    pred_covs_df_out_temp<- pred_covs_df_out_temp %>%
      mutate(., EST_YEAR = as.numeric(EST_YEAR),
             DATE = paste(EST_YEAR, case_when(
               SEASON == "Winter" ~ "12-16",
               SEASON == "Spring" ~ "03-16",
               SEASON == "Summer" ~ "07-16",
               SEASON == "Fall" ~ "09-16"), sep = "-"),
             SURVEY = "DUMMY",
             SVVESSEL = "DUMMY",
             NMFS_SVSPP = "DUMMY",
             DFO_SPEC = "DUMMY",
             PRESENCE = 1,
             BIOMASS = 1,
             ABUNDANCE = 1,
             ID = paste("DUMMY", DATE, sep = ""),
             PredTF = TRUE
      )
    
    if(i == 1){
      pred_covs_out<- pred_covs_df_out_temp
    } else {
      pred_covs_out<- bind_rows(pred_covs_out, pred_covs_df_out_temp)
    }
  }
  
  # Only going to keep information from fit_year_min through pred_years...
  pred_covs_out_final<- pred_covs_out %>%
    dplyr::filter(., EST_YEAR >= fit_year_min & EST_YEAR <= max(pred_years))
  
  # Some more work on columns to match up with `vast_data_out`
  pred_covs_out_final<- pred_covs_out_final %>%
    mutate(., VAST_YEAR_COV = ifelse(EST_YEAR > fit_year_max, fit_year_max, EST_YEAR),
           VAST_SEASON = case_when(
             SEASON == "Winter" ~ "DFO",
             SEASON == "Spring" ~ "SPRING",
             SEASON == "Summer" ~ "SUMMER",
             SEASON == "Fall" ~ "FALL"
           ),
           "VAST_YEAR_SEASON" = paste(EST_YEAR, VAST_SEASON, sep = "_"))
  
  # Subset to only seasons of interest...
  pred_covs_out_final<- pred_covs_out_final %>%
    filter(., VAST_SEASON %in% fit_seasons)
  
  # Need to account for new levels in year season...
  all_years<- seq(from = fit_year_min, to = max(pred_years), by = 1)
  all_seasons<- fit_seasons
  year_season_set<- expand.grid("SEASON" = all_seasons, "EST_YEAR" = all_years)
  all_year_season_levels<- apply(year_season_set[,2:1], MARGIN = 1, FUN = paste, collapse = "_")
  
  pred_covs_out_final<- pred_covs_out_final %>%
    mutate(., "VAST_YEAR_SEASON" = factor(VAST_YEAR_SEASON, levels = all_year_season_levels),
           "VAST_SEASON" = factor(VAST_SEASON, levels = all_seasons))
  
  # Name rearrangement!
  # Keep only what we need..
  cov_names<- names(pred_covs_out_final)[-which(names(pred_covs_out_final) %in% c("ID", "DATE", "EST_YEAR", "SEASON", "SURVEY", "SVVESSEL", "DECDEG_BEGLAT", "DECDEG_BEGLON", "NMFS_SVSPP", "DFO_SPEC", "PRESENCE", "BIOMASS", "ABUNDANCE", "PredTF", "VAST_YEAR_COV", "VAST_SEASON", "VAST_YEAR_SEASON"))]
  pred_covs_out_final<- pred_covs_out_final %>%
    dplyr::select(., "ID", "DATE", "EST_YEAR", "SEASON", "SURVEY", "SVVESSEL", "DECDEG_BEGLAT", "DECDEG_BEGLON", "NMFS_SVSPP", "DFO_SPEC", "PRESENCE", "BIOMASS", "ABUNDANCE", "PredTF", "VAST_YEAR_COV", "VAST_SEASON", "VAST_YEAR_SEASON", {{cov_names}})
  
  # Any extra covariates will likely be static...
  if(!is.null(extra_covariates_stack)){
    pred_covs_sf<- points_to_sf(pred_covs_out_final)
    
    pred_covs_out_final<- static_extract_wrapper(static_covariates_list = extra_covariates_stack, sf_points = pred_covs_sf, date_col_name = "DATE", df_sf = "df", out_dir = NULL)
  }
  
  # Drop any NAs, save and return it
  pred_covs_out_final<- pred_covs_out_final %>%
    drop_na(., {{cov_names}}) %>%
    mutate(., "Summarized" = summarize,
           "Ensemble_Stat" = ensemble_stat)
  
  saveRDS(pred_covs_out_final, file = paste(out_dir, "/VAST_pred_df_", summarize, "_", ensemble_stat, ".rds", sep = "" ))
  return(pred_covs_out_final)
}

#' @title Make VAST seasonal dataset
#' 
#' @description This function reads in a tidy model dataset and does some cleaning and processing to generate a new dataset to accommodate fitting a VAST seasonal (or other intra annual) model. These cleaning and processing steps boil down to creating an ordered, continuous, season-year vector, such that the model can then estimate density even in season-years not surveyed.
#'
#' @param tidy_mod_data = A tidy model datafame with all the information (tows, habitat covariates, species occurrences) needed to fit a species distribution model. 
#' @param nmfs_species_code = Numeric NMFS species code
#' @param fit_year_min = Minimum year to keep
#' @param fit_year_max = Maximum year to keep
#' @param pred_df = Either NULL or a dataframe with prediction information as created by `make_vast_predict_df`
#' @param out_dir = Directory to save the tidy model dataframe as an .rds file
#' 
#' @return  A VAST seasonal dataset, ready to be split into a `sample data` dataframe and a `covariate data` dataframe. This file is also saved in out_dir. 
#' 
#' @export

make_vast_seasonal_data<- function(tidy_mod_data, fit_seasons, nmfs_species_code, fit_year_min, fit_year_max, pred_df, out_dir){
  
  # For debugging
  if(FALSE){
    tar_load(tidy_mod_data)
    nmfs_species_code = nmfs_species_code
    fit_year_min = fit_year_min
    fit_year_max = fit_year_max
    fit_seasons = fit_seasons
    pred_years = pred_years
    tar_load(vast_predict_df)
    pred_df = vast_predict_df
    out_dir = here::here("scratch/aja/targets_flow/data/combined/")
  }
  
  # Some work on the time span and seasons
  data_temp<- tidy_mod_data %>%
    filter(., NMFS_SVSPP == nmfs_species_code) %>%
    filter(., EST_YEAR >= fit_year_min & EST_YEAR <= fit_year_max) %>%
    mutate(., "VAST_SEASON" = case_when(
      SURVEY == "DFO" & SEASON == "SPRING" ~ "DFO",
      SURVEY == "NMFS" & SEASON == "SPRING" ~ "SPRING",
      SURVEY == "DFO" & SEASON == "SUMMER" ~ "SUMMER",
      SURVEY == "NMFS" & SEASON == "FALL" ~ "FALL")) %>%
    drop_na(VAST_SEASON)
  
  data_temp<- data_temp %>%
    filter(., VAST_SEASON %in% fit_seasons)

  # Set of years and seasons. The DFO spring survey usually occurs before the NOAA NEFSC spring survey, so ordering accordingly. Pred year max or fit year max??
  all_years<- seq(from = fit_year_min, to = fit_year_max, by = 1)
  all_seasons<- fit_seasons
  yearseason_set<- expand.grid("SEASON" = all_seasons, "EST_YEAR" = all_years)
  all_yearseason_levels<- apply(yearseason_set[,2:1], MARGIN = 1, FUN = paste, collapse = "_")
  
  # year_set<- sort(unique(data_temp$EST_YEAR))
  # season_set<- c("DFO", "SPRING", "FALL")
  # 
  # # Create a grid with all unique combinations of seasons and years and then combine these into one "year_season" variable
  # yearseason_grid<- expand.grid("SEASON" = season_set, "EST_YEAR" = year_set)
  # yearseason_levels<- apply(yearseason_grid[, 2:1], MARGIN = 1, FUN = paste, collapse = "_")
  # yearseason_labels<- round(yearseason_grid$EST_YEAR + (as.numeric(factor(yearseason_grid$VAST_SEASON, levels = season_set))-1)/length(season_set), digits = 1)
  # 
  # Similar process, but for the observations
  yearseason_i<- apply(data_temp[, c("EST_YEAR", "VAST_SEASON")], MARGIN = 1, FUN = paste, collapse = "_")
  yearseason_i<- factor(yearseason_i, levels = all_yearseason_levels)
  
  # Add the year_season factor column to our sampling_data data set
  data_temp$VAST_YEAR_SEASON<- yearseason_i
  data_temp$VAST_SEASON = factor(data_temp$VAST_SEASON, levels = all_seasons)
  
  # VAST year
  data_temp$VAST_YEAR_COV<- ifelse(data_temp$EST_YEAR > fit_year_max, fit_year_max, data_temp$EST_YEAR)
  data_temp$PredTF<- FALSE
  
  # Ordering...
  cov_names<- names(data_temp)[-which(names(data_temp) %in% c("ID", "DATE", "EST_YEAR", "SEASON", "SURVEY", "SVVESSEL", "DECDEG_BEGLAT", "DECDEG_BEGLON", "NMFS_SVSPP", "DFO_SPEC", "PRESENCE", "BIOMASS", "ABUNDANCE", "PredTF", "VAST_YEAR_COV", "VAST_SEASON", "VAST_YEAR_SEASON"))]
  cov_names<- cov_names[-which(cov_names == "Season_Match")]
  data_temp<- data_temp %>%
    dplyr::select("ID", "DATE", "EST_YEAR", "SEASON", "SURVEY", "SVVESSEL", "DECDEG_BEGLAT", "DECDEG_BEGLON", "NMFS_SVSPP", "DFO_SPEC", "PRESENCE", "BIOMASS", "ABUNDANCE", "PredTF", "VAST_YEAR_COV", "VAST_SEASON", "VAST_YEAR_SEASON", {{cov_names}})
  
  # Make dummy data for all year_seasons to estimate gaps in sampling if needed
  dummy_data<- data.frame("ID" = sample(data_temp$ID, size = 1), "DATE" = mean(data_temp$DATE, na.rm = TRUE), "EST_YEAR" = yearseason_set[,'EST_YEAR'], "SEASON" = yearseason_set[,'SEASON'], "SURVEY" = "DUMMY", "SVVESSEL" = "DUMMY", "DECDEG_BEGLAT" = mean(data_temp$DECDEG_BEGLAT, na.rm = TRUE), "DECDEG_BEGLON" = mean(data_temp$DECDEG_BEGLON, na.rm = TRUE), "NMFS_SVSPP" = "DUMMY", "DFO_SPEC" = "DUMMY", "PRESENCE" = 1, "BIOMASS" = 1, "ABUNDANCE" = 1, "PredTF" = TRUE, "VAST_YEAR_COV" = yearseason_set[,'EST_YEAR'], "VAST_SEASON" = yearseason_set[,'SEASON'], "VAST_YEAR_SEASON" = all_yearseason_levels)
  
  # Add in "covariates"
  col_ind<- ncol(dummy_data)
  for(i in seq_along(cov_names)){
    col_ind<- col_ind+1
    cov_vec<- unlist(data_temp[,{{cov_names}}[i]])
    dummy_data[,col_ind]<- mean(cov_vec, na.rm = TRUE)
    names(dummy_data)[col_ind]<- {{cov_names}}[i]
  }
  
  # Combine with original dataset
  vast_data_out<- rbind(data_temp, dummy_data)
  vast_data_out$VAST_YEAR_COV<- factor(vast_data_out$VAST_YEAR_COV, levels = seq(from = fit_year_min, to = fit_year_max, by = 1))
  
  # If we have additional years that we want to predict to and NOT Fit too, we aren't quite done just yet...
  if(!is.null(pred_df)){
    # Name work...
    pred_df<- pred_df %>%
      dplyr::select(., -Summarized, -Ensemble_Stat)
    
    # Add those -- check names first
    check_names<- all(colnames(pred_df) %in% colnames(vast_data_out)) & all(colnames(vast_data_out) %in% colnames(pred_df))
    if(!check_names){
      print("Check data and prediction column names, they don't match")
      stop()
    } else {
      # We only need one observation for each of the times...
      pred_df_bind<- pred_df %>%
        dplyr::select(., colnames(vast_data_out)) %>%
        distinct(., ID, .keep_all = TRUE)
      vast_data_out<- rbind(vast_data_out, pred_df_bind)
    }
  }

  # Save and return it
  saveRDS(vast_data_out, file = paste(out_dir, "vast_data.rds", sep = "/"))
  return(vast_data_out)
}

#' @title Make VAST sample dataset
#' 
#' @description This function creates a VAST sample dataset to pass into calls to `VAST::fit_model`.
#'
#' @param vast_seasonal_data = Description
#' @param out_dir = Description
#' 
#' @return A sample dataframe that includes all of the "sample" or species occurrence information. This file is also saved in out_dir. 
#' 
#' @export

make_vast_sample_data<- function(vast_seasonal_data, fit_seasons, out_dir){
  
  # For debugging
  if(FALSE){
    tar_load(vast_seasonal_data)
    out_dir = here::here("scratch/aja/targets_flow/data/dfo/combined")
  }
  
  # Select columns we want from the "full" vast_seasonal_data dataset
  vast_samp_dat<- data.frame(
    "Year" = as.numeric(vast_seasonal_data$VAST_YEAR_SEASON)-1, 
    "Lat" = vast_seasonal_data$DECDEG_BEGLAT,
    "Lon" = vast_seasonal_data$DECDEG_BEGLON,
    "Biomass" = vast_seasonal_data$BIOMASS,
    "Swept" = rep(0.5, nrow(vast_seasonal_data)),
    "Pred_TF" = vast_seasonal_data$PredTF
  )
  
  # Save and return it
  saveRDS(vast_samp_dat, file = paste(out_dir, "vast_sample_data.rds", sep = "/"))
  return(vast_samp_dat)
}

#' @title Make VAST covariate dataset
#' 
#' @description This function creates a VAST covariate dataset to pass into calls to `VAST::fit_model`.
#'
#' @param vast_seasonal_data = Description 
#' @param out_dir = Description
#' 
#' @return A sample dataframe that includes all of the covariate information at each unique sample. This file is also saved in out_dir. 
#' 
#' @export

make_vast_covariate_data<- function(vast_seasonal_data, out_dir){
  
  # For debugging
  if(FALSE){
    vast_seasonal_data
    out_dir = here::here("scratch/aja/targets_flow/data/dfo/combined")
  }
  
  # Select columns we want from the "full" vast_seasonal_data dataset
  vast_cov_dat<- data.frame(
    "Year" = as.numeric(vast_seasonal_data$VAST_YEAR_SEASON)-1,
    "Year_Cov" = vast_seasonal_data$VAST_YEAR_COV,
    "Season" = vast_seasonal_data$VAST_SEASON,
    "SST_seasonal" = vast_seasonal_data$SST_seasonal,
    "BT_seasonal" = vast_seasonal_data$BT_seasonal,
    "Lat" = vast_seasonal_data$DECDEG_BEGLAT,
    "Lon" = vast_seasonal_data$DECDEG_BEGLON
  )
  
  # Save and return it
  saveRDS(vast_cov_dat, file = paste(out_dir, "vast_covariate_data.rds", sep = "/"))
  return(vast_cov_dat)
}

#' @title Read in shapefile
#' 
#' @description A short function to read in a shapefile given a file path
#'
#' @param file_path = File path to geospatial vector polygon file with .shp extension, specifying the location and shape of the area of interest.
#' @param factor_vars = Names of factor columns that should be checked and converted if necessary
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
#' @param shapefile = A geospatial vector sf polygon file, specifying the location and shape of the area of interest.
#' @param cell_size = The size of grid in meters (since working in UTM). This will control the resolution of the extrapolation grid.
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
#' @param extrap_grid = User created extrapolation grid from vast_make_extrap_grid.
#' @param FieldConfig = A vector defining the number of spatial (Omega) and spatio-temporal (Epsilon) factors to include in the model for each of the linear predictors. For each factor, possible values range from 0 (which effectively turns off a given factor), to the number of categories being modeled. If FieldConfig < number of categories, VAST estimates common factors and then loading matrices. 
#' @param RhoConfig = A vector defining the temporal structure of intercepts (Beta) and spatio-temporal (Epsilon) variation for each of the linear predictors. See `VAST::make_data` for options.
#' @param bias.correct = Logical boolean determining if Epsilon bias-correction should be done. 
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
  settings_out<- make_settings(n_x = n_x_use, Region = "User", purpose = "index2", FieldConfig = FieldConfig, RhoConfig = RhoConfig, ObsModel = c(1, 1), bias.correct = bias.correct, knot_method = "grid", treat_nonencounter_as_zero = FALSE)
  
  # Return it
  return(settings_out)
}

####
#' @title Make VAST covariate effect objects
#' 
#' @description Create covariate effects for both linear predictors
#'
#' @param X1_coveff_vec = A vector specifying the covariate effects for first linear predictor. 
#' @param X2_coveff_vec = A vector specifying the covariate effects for second linear predictor. 
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
#' @param settings = A tagged list with the settings for the model, created with `vast_make_settings`.
#' @param input_grid = An extrapolation grid, created with `vast_make_extrap_grid`.
#' @param samp_dat = A data frame with the biomass sample data for each species at each tow.
#' @param cov_dat = A data frame with the covariate data for each tow.
#' @param X1_formula = A formula for the first linear predictor.
#' @param X2_formula = A formula for the second linear predictor.
#' @param X_contrasts = A tagged list specifying the contrasts to use for factor covariates in the model.
#' @param Xconfig_list = A tagged list specifying the covariate effects for first and second linear predictors.
#'
#' @return A VAST `fit_model` object, with the inputs and built TMB object components.
#'
#' @export

vast_build_sdm <- function(settings, extrap_grid, sample_data, covariate_data, X1_formula, X2_formula, X_contrasts, Xconfig_list){
  
  # For debugging
  if(FALSE){
    library(VAST)
    library(tidyverse)
    library(stringr)
    library(spli)
    
    tar_load(vast_settings)
    settings = vast_settings
    tar_load(vast_extrap_grid)
    extrap_grid = vast_extrap_grid
    tar_load(vast_sample_data)
    sample_data = vast_sample_data
    
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
  samp_dat_names<- c("Lat", "Lon", "Year", "Biomass", "Swept", "Pred_TF")
  if(!(all(samp_dat_names %in% names(sample_data)))){
    stop(paste("Check names in sample data. Must include:", paste0(samp_dat_names, collapse = ","), sep = " "))
  }
  
  cov_dat_names1<- unlist(str_extract_all(X1_formula, boundary("word"))[[2]])
  
  # Remove some stuff associated with the splines...
  spline_words<- c("bs", "degree", "TRUE", "intercept", "2", "FALSE")
  cov_dat_names1<- cov_dat_names1[-which(cov_dat_names1 %in% spline_words)]
  cov_dat_names2<- unlist(str_extract_all(X2_formula, boundary("word"))[[2]])
  cov_dat_names2<- cov_dat_names2[-which(cov_dat_names2 %in% spline_words)]
  cov_dat_names_all<- unique(c(cov_dat_names1, cov_dat_names2))
  if(!(all(cov_dat_names_all %in% names(covariate_data)))){
    stop(paste("Check names in covariate data. Must include", paste0(cov_dat_names_all, collapse = ","), sep = " "))
  }
  
  if(!(all(c("X1config_cp", "X2config_cp") %in% names(Xconfig_list)))){
    stop(paste("Check names of Xconfig_list. Must be", paste0(c("X1config_cp", "X2config_cp"), collapse = ","), sep = ""))
  }
  
  # Run VAST::fit_model with correct info and settings
  vast_build_out<- fit_model("settings" = settings, input_grid = as.matrix(extrap_grid, ncol = 3), "Lat_i" = sample_data[, 'Lat'], "Lon_i" = sample_data[, 'Lon'], "t_i" = sample_data[, 'Year'], "c_i" = rep(0, nrow(sample_data)), "b_i" = sample_data[, 'Biomass'], "a_i" = sample_data[, 'Swept'], "PredTF_i" = sample_data[, 'Pred_TF'], "X1config_cp" = Xconfig_list[['X1config_cp']], "X2config_cp" = Xconfig_list[['X2config_cp']], "covariate_data" = covariate_data, "X1_formula" = X1_formula, "X2_formula" = X2_formula, X_contrasts = X_contrasts, "newtonsteps" = 1, "getsd" = TRUE, "getReportCovariance" = TRUE, "run_model" = FALSE, "test_fit" = FALSE,  "Use_REML" = FALSE, "getJointPrecision" = FALSE)
  
  # Return it
  return(vast_build_out)
}

####
#' @title Adjust VAST SDM
#' 
#' @description Make adjustments to VAST SDM and the model returned in `vast_build_sdm`. This can either be the exact same as the one built using `vast_build_sdm`, or it can update that model with adjustments provided in a tagged list. 
#'
#' @param vast_build = A VAST `fit_model` object.
#' @param adjustments = Either NULL (default) or a tagged list identifying adjustments that should be made to the vast_build `fit_model` object. If NULL, the identical model defined by the `vast_build` is run and fitted.
#'
#' @return A VAST fit_model object, with the inputs and built TMB object components.
#' 
#' @export

vast_make_adjustments <- function(vast_build, adjustments = NULL){
  
  # For debugging
  if(FALSE){
    tar_load(vast_build0)
    vast_build = vast_build0
    adjustments = list("log_sigmaXi1_cp" = factor(c(rep(1, 3), rep(4, nlevels(vast_covariate_data$Year_Cov)), rep(NA, 2))), "log_sigmaXi2_cp" = factor(c(rep(1, 3), rep(4, nlevels(vast_covariate_data$Year_Cov)), rep(NA, 2))))
  }
  
  # If no adjustments are needed, just need to pull information from vast_build and then set "run_model" to TRUE
  if(is.null(adjustments)){
    vast_build_adjust_out<- fit_model("settings" = vast_build$settings, input_grid = vast_build$input_args$data_args_input$input_grid, "Lat_i" = vast_build$data_frame[, 'Lat_i'], "Lon_i" = vast_build$data_frame[, 'Lon_i'], "t_i" = vast_build$data_frame[, 't_i'], "c_iz" = vast_build$data_frame[, 'c_iz'], "b_i" = vast_build$data_frame[, 'b_i'], "a_i" = vast_build$data_frame[, 'a_i'], "PredTF_i" = vast_build$data_list[['PredTF_i']], "X1config_cp" = vast_build$input_args$data_args_input[['X1config_cp']], "X2config_cp" = vast_build$input_args$data_args_input[['X2config_cp']], "covariate_data" = vast_build$input_args$data_args_input$covariate_data, "X1_formula" = vast_build$input_args$data_args_input$X1_formula, "X2_formula" = vast_build$input_args$data_args_input$X2_formula, X_contrasts = vast_build$input_args$data_args_input$X_contrasts, "newtonsteps" = 1, "getsd" = TRUE, "getReportCovariance" = TRUE, "run_model" = FALSE, "test_fit" = FALSE,  "Use_REML" = FALSE, "getJointPrecision" = FALSE)
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
      vast_build_adjust_out<- fit_model("settings" = vast_build$settings, input_grid = vast_build$input_args$data_args_input$input_grid, "Lat_i" = vast_build$data_frame[, 'Lat_i'], "Lon_i" = vast_build$data_frame[, 'Lon_i'], "t_i" = vast_build$data_frame[, 't_i'], "c_iz" = vast_build$data_frame[, 'c_iz'], "b_i" = vast_build$data_frame[, 'b_i'], "a_i" = vast_build$data_frame[, 'a_i'], "PredTF_i" = vast_build$data_list[['PredTF_i']], "X1config_cp" = vast_build$input_args$data_args_input[['X1config_cp']], "X2config_cp" = vast_build$input_args$data_args_input[['X2config_cp']], "covariate_data" = vast_build$input_args$data_args_input$covariate_data, "X1_formula" = vast_build$input_args$data_args_input$X1_formula, "X2_formula" = vast_build$input_args$data_args_input$X2_formula, X_contrasts = vast_build$input_args$data_args_input$X_contrasts, "Map" = map_adjust_out, "newtonsteps" = 1, "getsd" = TRUE, "getReportCovariance" = TRUE, "run_model" = FALSE, "test_fit" = FALSE,  "Use_REML" = FALSE, "getJointPrecision" = FALSE)
    } else {
      # No need for Map argument, just build and fit
      vast_build_adjust_out<- fit_model("settings" = vast_build$settings, input_grid = vast_build$input_args$data_args_input$input_grid, "Lat_i" = vast_build$data_frame[, 'Lat_i'], "Lon_i" = vast_build$data_frame[, 'Lon_i'], "t_i" = vast_build$data_frame[, 't_i'], "c_iz" = vast_build$data_frame[, 'c_iz'], "b_i" = vast_build$data_frame[, 'b_i'], "a_i" = vast_build$data_frame[, 'a_i'], "PredTF_i" = vast_build$data_list[['PredTF_i']], "X1config_cp" = vast_build$input_args$data_args_input[['X1config_cp']], "X2config_cp" = vast_build$input_args$data_args_input[['X2config_cp']], "covariate_data" = vast_build$input_args$data_args_input$covariate_data, "X1_formula" = vast_build$input_args$data_args_input$X1_formula, "X2_formula" = vast_build$input_args$data_args_input$X2_formula, X_contrasts = vast_build$input_args$data_args_input$X_contrasts, "newtonsteps" = 1, "getsd" = TRUE, "getReportCovariance" = TRUE, "run_model" = FALSE, "test_fit" = FALSE,  "Use_REML" = FALSE, "getJointPrecision" = FALSE)
    }
  }
  # Return it
  return(vast_build_adjust_out)
}

#' @title Fit VAST SDM
#' 
#' @description Fit VAST species distribution model
#'
#' @param vast_build = A VAST `fit_model` object.
#'
#' @return A VAST fit_model object, with the inputs and and outputs, including parameter estimates, extrapolation gid info, spatial list info, data info, and TMB info.
#'
#' @export

vast_fit_sdm <- function(vast_build_adjust, nmfs_species_code, out_dir){
  
  # For debugging
  if(FALSE){
    tar_load(vast_adjust)
    vast_build_adjust = vast_adjust
  }
  
  # Build and fit model
  vast_fit_out<- fit_model("settings" = vast_build_adjust$settings, input_grid = vast_build_adjust$input_args$data_args_input$input_grid, "Lat_i" = vast_build_adjust$data_frame[, 'Lat_i'], "Lon_i" = vast_build_adjust$data_frame[, 'Lon_i'], "t_i" = vast_build_adjust$data_frame[, 't_i'], "c_iz" = vast_build_adjust$data_frame[, 'c_iz'], "b_i" = vast_build_adjust$data_frame[, 'b_i'], "a_i" = vast_build_adjust$data_frame[, 'a_i'], "PredTF_i" = vast_build_adjust$data_list[['PredTF_i']], "X1config_cp" = vast_build_adjust$input_args$data_args_input[['X1config_cp']], "X2config_cp" = vast_build_adjust$input_args$data_args_input[['X2config_cp']], "covariate_data" = vast_build_adjust$input_args$data_args_input$covariate_data, "X1_formula" = vast_build_adjust$input_args$data_args_input$X1_formula, "X2_formula" = vast_build_adjust$input_args$data_args_input$X2_formula, "X_contrasts" = vast_build_adjust$input_args$data_args_input$X_contrasts, "Map" = vast_build_adjust$tmb_list$Map, "newtonsteps" = 1, "getsd" = TRUE, "getReportCovariance" = TRUE, "run_model" = TRUE, "test_fit" = FALSE,  "Use_REML" = FALSE, "getJointPrecision" = FALSE)
  
  # Save and return it
  saveRDS(vast_fit_out, file = paste(out_dir, "/", nmfs_species_code, "_", "fitted_vast.rds", sep = "" ))
  return(vast_fit_out)
}

#' @title Plot VAST model predicted density surfaces
#' 
#' @description Creates either a panel plot or a gif of VAST model predicted density surfaces
#'
#' @param vast_fit = A VAST `fit_model` object.
#' @param all_times = A vector of all of the unique time steps available from the VAST fitted model
#' @param plot_times = Either NULL to make a plot for each time in `all_times` or a vector of all of the times to plot, which must be a subset of `all_times`
#' @param land_sf = Land sf object
#' @param xlim = A two element vector with the min and max longitudes 
#' @param ylim = A two element vector with the min and max latitudes 
#' @param panel_or_gif = A character string of either "panel" or "gif" indicating how the multiple plots across time steps should be displayed
#' @param out_dir = Output directory to save the panel plot or gif
#' 
#' @return A VAST fit_model object, with the inputs and and outputs, including parameter estimates, extrapolation gid info, spatial list info, data info, and TMB info.
#'
#' @export
vast_fit_plot_density<- function(vast_fit, nmfs_species_code, mask, all_times = all_times, plot_times = NULL, land_sf, xlim, ylim, panel_or_gif = "gif", out_dir, land_color = "#d9d9d9", panel_cols = NULL, panel_rows = NULL, ...){
  if(FALSE){
    tar_load(vast_fit)
    template = raster("~/GitHub/sdm_workflow/scratch/aja/TargetsSDM/data/supporting/HighResTemplate.grd")
    tar_load(vast_seasonal_data)
    all_times = as.character(levels(vast_seasonal_data$YEAR_SEASON))
    plot_times = NULL
    tar_load(land_sf)
    tar_load(shapefile)
    mask = shapefile
    land_color = "#d9d9d9"
    res_data_path = "~/Box/RES_Data/"
    xlim = c(-85, -55)
    ylim = c(30, 50)
    panel_or_gif = "gif"
    panel_cols = NULL
    panel_rows = NULL
  }
  
  # Plotting at spatial knots...
  # Getting prediction array
  pred_array<- log(vast_fit$Report$D_gct+1)
  
  # Getting time info
  if(!is.null(plot_times)){
    plot_times<- all_times[which(all_times) %in% plot_times]
  } else {
    plot_times<- all_times
  }
  
  # Getting spatial information
  spat_data<- vast_fit$extrapolation_list
  loc_g<- spat_data$Data_Extrap[which(spat_data$Data_Extrap[, "Include"] > 0), c("Lon", "Lat")]
  CRS_orig<- sp::CRS("+proj=longlat")
  CRS_proj<- sp::CRS(spat_data$projargs)
  land_sf<- st_crop(land_sf, xmin = xlim[1], ymin = ylim[1], xmax = xlim[2], ymax = ylim[2])
  
  # Looping through...
  rasts_out<- vector("list", dim(pred_array)[3])
  rasts_range<- pred_array
  rast_lims<- c(round(min(rasts_range)-0.000001, 2), round(max(rasts_range) + 0.0000001, 2))
  
  if(dim(pred_array)[3] == 1){
    df<- data.frame(loc_g, z = pred_array[,1,])
    points_ll = st_as_sf(data_df, coords = c("Lon", "Lat"), crs = CRS_orig)
    points_proj = points_ll %>%
      st_transform(., crs = CRS_proj)
    points_bbox<- st_bbox(points_proj)
    raster_proj<- st_rasterize(points_proj)
    raster_proj<- resample(raster_proj, raster(template))
    
    plot_out<- ggplot() +
      geom_stars(data = raster_proj, aes(x = x, y = y, fill = z)) +
      scale_fill_viridis_c(name = "Density", option = "viridis", na.value = "transparent", limits = rast_lims) +
      geom_sf(data = land_sf_proj, fill = land_color, lwd = 0.2) +
      coord_sf(xlim = points_bbox[c(1,3)], ylim = points_bbox[c(2,4)], expand = FALSE, datum = sf::st_crs(CRS_proj)) 
    theme(panel.background = element_rect(fill = "white"), panel.border = element_rect(fill = NA), axis.text.x=element_blank(), axis.text.y=element_blank(), axis.ticks=element_blank(), axis.title = element_blank(), plot.margin = margin(t = 0.05, r = 0.05, b = 0.05, l = 0.05))
    
    ggsave(filename = paste(out_dir, file_name, ".png", sep = ""), plot_out, width = 11, height = 8, units = "in")
  } else {
    
    for (tI in 1:dim(pred_array)[3]) {
      data_df<- data.frame(loc_g, z = pred_array[,1,tI])
      
      # Interpolation
      pred_df<- na.omit(data.frame("x" = data_df$Lon, "y" = data_df$Lat, "layer" = data_df$z))
      pred_df_interp<- interp(pred_df[,1], pred_df[,2], pred_df[,3], duplicate = "mean", extrap = TRUE,
                              xo=seq(-87.99457, -57.4307, length = 115),
                              yo=seq(22.27352, 48.11657, length = 133))
      pred_df_interp_final<- data.frame(expand.grid(x = pred_df_interp$x, y = pred_df_interp$y), z = c(round(pred_df_interp$z, 2)))
      pred_sp<- st_as_sf(pred_df_interp_final, coords = c("x", "y"), crs = CRS_orig)
      
      pred_df_temp<- pred_sp[which(st_intersects(pred_sp, mask, sparse = FALSE) == TRUE),]
      coords_keep<- as.data.frame(st_coordinates(pred_df_temp))
      row.names(coords_keep)<- NULL
      pred_df_use<- data.frame(cbind(coords_keep, "z" = as.numeric(pred_df_temp$z)))
      names(pred_df_use)<- c("x", "y", "z")
      
      # raster_proj<- raster::rasterize(as_Spatial(points_ll), template, field = "z", fun = mean)
      # raster_proj<- as.data.frame(raster_proj, xy = TRUE)
      # 
      time_plot_use<- plot_times[tI]
      
      rasts_out[[tI]]<- ggplot() +
        geom_tile(data = pred_df_use, aes(x = x, y = y, fill = z)) +
        scale_fill_viridis_c(name = "Log (density+1)", option = "viridis", na.value = "transparent", limits = rast_lims) +
        annotate("text", x = -65, y = 37.5, label = time_plot_use) +
        geom_sf(data = land_sf, fill = land_color, lwd = 0.2, na.rm = TRUE) +
        coord_sf(xlim = xlim, ylim = ylim, expand = FALSE) +
      theme(panel.background = element_rect(fill = "white"), panel.border = element_rect(fill = NA), axis.text.x=element_blank(), axis.text.y=element_blank(), axis.ticks=element_blank(), axis.title = element_blank(), plot.margin = margin(t = 0.05, r = 0.05, b = 0.05, l = 0.05, unit = "pt"))
    }
    if(panel_or_gif == "panel"){
      # Panel plot
      all_plot<- wrap_plots(rasts_out, ncol = panel_cols, nrow = panel_rows, guides = "collect", theme(plot.margin = margin(t = 0.05, r = 0.05, b = 0.05, l = 0.05, unit = "pt")))
      ggsave(filename = paste(working_dir, file_name, ".png", sep = ""), all.plot, width = 11, height = 8, units = "in")
    } else {
      # Make a gif
      plot_loop_func<- function(plot_list){
        for (i in seq_along(plot_list)) {
          plot_use<- plot_list[[i]]
          print(plot_use)
        }
      }
      invisible(save_gif(plot_loop_func(rasts_out), paste0(out_dir, nmfs_species_code, "_LogDensity.gif"), delay = 0.75, progress = FALSE))
    }
  }
}

#' @title Predict fitted VAST model
#'
#' @description This function makes predictions from a fitted VAST SDM to new locations using VAST::predict.fit_model. Importantly, to use this feature for new times, at least one location for each time of interest needs to be included during the model fitting process. This dummy observation should have a PredTF value of 1 so that the observation is only used in the predicted probability and NOT eastimating the likelihood.
#'
#' @param vast_fitted_sdm = A fitted VAST SDM object, as returned with `vast_fit_sdm`
#' @param predict_variable = Which variable should be predicted, default is density (D_i)
#' @param predict_category = Which category (species/age/size) should be predicted, default is 0
#' @param predict_vessel = Which sampling category should be predicted, default is 0
#' @param predict_covariates_df = A stack of the prediction covariates
#' @param out_dir = Output directory to save...
#' 
#' @return
#'
#' @export
predict_vast<- function(vast_fitted_sdm, nmfs_species_code, predict_variable = "D_i", predict_category = 0, predict_vessel = 0, predict_covariates_df_all, out_dir){

  # For debugging
  if(FALSE){
    tar_load(vast_fit)
    predict_variable = "D_i"
    predict_covariates_df_all = pred_covs_out_final
    predict_category = 0
    predict_vessel = 0
  }
  
  # Collecting necessary bits from the prediction covariates -- lat, lon, time
  pred_lats<- predict_covariates_df_all$DECDEG_BEGLAT
  pred_lons<- predict_covariates_df_all$DECDEG_BEGLON
  pred_times<- as.numeric(predict_covariates_df_all$VAST_YEAR_SEASON)-1
  
  # Catch stuff...
  pred_sampled_areas<- rep(1, length(pred_lats))
  pred_category<- rep(predict_category, length(pred_lats))
  pred_vessel<- rep(predict_vessel, length(pred_lats))
  
  # Renaming predict_covariates_df_all to match vast_fit_covariate_data
  pred_cov_dat_use<- data.frame(
    "Year" = as.numeric(predict_covariates_df_all$VAST_YEAR_SEASON)-1,
    "Year_Cov" = predict_covariates_df_all$VAST_YEAR_COV,
    "Season" = predict_covariates_df_all$VAST_SEASON,
    "SST_seasonal" = predict_covariates_df_all$SST_seasonal,
    "BT_seasonal" = predict_covariates_df_all$BT_seasonal,
    "Lat" = predict_covariates_df_all$DECDEG_BEGLAT,
    "Lon" = predict_covariates_df_all$DECDEG_BEGLON
  )
  
  # Make the predictions
  preds_out<- predict(x = vast_fitted_sdm, what = predict_variable, Lat_i = pred_lats, Lon_i = pred_lons, t_i = pred_times, a_i = pred_sampled_areas, c_iz = pred_category, v_i = pred_vessel, new_covariate_data = pred_cov_dat_use, do_checks = FALSE)
  
  # Get everything as a dataframe to make plotting easier...
  pred_df_out<- data.frame("Lat" = pred_lats, "Lon" = pred_lons, "Time" = predict_covariates_df_all$VAST_YEAR_SEASON, "Pred" = preds_out)
  
  # Save and return it
  saveRDS(pred_df_out, file = paste(out_dir, "/pred_Di_", nmfs_species_code, "_", unique(predict_covariates_df_all$Summarized), "_", unique(predict_covariates_df_all$Ensemble_Stat), ".rds", sep = "" ))
  return(pred_df_out)
}


#' @title Plot predicted density surfaces from data frame
#' 
#' @description Creates either a panel plot or a gif of predicted density surfaces from a data frame that has location and time information
#'
#' @param pred_df = A dataframe with Lat, Lon, Time and Pred columns
#' @param mask = Land mask
#' @param plot_times = Either NULL to make a plot for each time in `pred_df$Time` or a vector of all of the times to plot, which must be a subset of `pred_df$Time`
#' @param land_sf = Land sf object
#' @param xlim = A two element vector with the min and max longitudes 
#' @param ylim = A two element vector with the min and max latitudes 
#' @param panel_or_gif = A character string of either "panel" or "gif" indicating how the multiple plots across time steps should be displayed
#' @param out_dir = Output directory to save the panel plot or gif
#' 
#' @return NULL. Panel or gif plot is saved in out_dir.
#'
#' @export
df_plot_density<- function(pred_df, nmfs_species_code, mask, plot_times = NULL, land_sf, xlim, ylim, panel_or_gif = "gif", out_dir, land_color = "#d9d9d9", panel_cols = NULL, panel_rows = NULL, ...){
  if(FALSE){
    tar_load(vast_predictions)
    pred_df = vast_predictions
    plot_times = NULL
    tar_load(land_sf)
    tar_load(shapefile)
    mask = shapefile
    land_color = "#d9d9d9"
    res_data_path = "~/Box/RES_Data/"
    xlim = c(-80, -55)
    ylim = c(35, 50)
    panel_or_gif = "gif"
    panel_cols = NULL
    panel_rows = NULL
  }
  
  # Time ID column for filtering
  pred_df<- pred_df %>%
    mutate(., "Time_Filter" = as.numeric(Time))
  
  # Log transform pred_df$Pred 
  pred_df$Pred<- log(pred_df$Pred+1)
  
  # Getting all unique times
  all_times<- unique(pred_df$Time)
  
  # Getting time info
  if(!is.null(plot_times)){
    plot_times<- all_times[which(all_times) %in% plot_times]
  } else {
    plot_times<- all_times
  }
  
  # Getting spatial information
  land_sf<- st_crop(land_sf, xmin = xlim[1], ymin = ylim[1], xmax = xlim[2], ymax = ylim[2])
  
  # Looping through...
  rasts_out<- vector("list", length(plot_times))
  rasts_range<- pred_df$Pred
  rast_lims<- c(round(min(rasts_range)-0.000001, 2), round(max(rasts_range) + 0.0000001, 2))
  
  for (tI in 1:length(plot_times)) {
    pred_df_temp<- pred_df %>%
      dplyr::filter(., Time_Filter == tI)
    
    # Interpolation
    pred_df_temp<- na.omit(data.frame("x" = pred_df_temp$Lon, "y" = pred_df_temp$Lat, "layer" = pred_df_temp$Pred))
    pred_df_interp<- interp(pred_df_temp[,1], pred_df_temp[,2], pred_df_temp[,3], duplicate = "mean", extrap = TRUE,
                            xo=seq(-87.99457, -57.4307, length = 115),
                            yo=seq(22.27352, 48.11657, length = 133))
    pred_df_interp_final<- data.frame(expand.grid(x = pred_df_interp$x, y = pred_df_interp$y), z = c(round(pred_df_interp$z, 2)))
    pred_sp<- st_as_sf(pred_df_interp_final, coords = c("x", "y"), crs = 4326)
    
    pred_df_temp2<- pred_sp[which(st_intersects(pred_sp, mask, sparse = FALSE) == TRUE),]
    coords_keep<- as.data.frame(st_coordinates(pred_df_temp2))
    row.names(coords_keep)<- NULL
    pred_df_use<- data.frame(cbind(coords_keep, "z" = as.numeric(pred_df_temp2$z)))
    names(pred_df_use)<- c("x", "y", "z")
    
    # raster_proj<- raster::rasterize(as_Spatial(points_ll), template, field = "z", fun = mean)
    # raster_proj<- as.data.frame(raster_proj, xy = TRUE)
    # 
    time_plot_use<- plot_times[tI]
    
    rasts_out[[tI]]<- ggplot() +
      geom_tile(data = pred_df_use, aes(x = x, y = y, fill = z)) +
      scale_fill_viridis_c(name = "Log (density+1)", option = "viridis", na.value = "transparent", limits = rast_lims) +
      annotate("text", x = -65, y = 37.5, label = time_plot_use) +
      geom_sf(data = land_sf, fill = land_color, lwd = 0.2, na.rm = TRUE) +
      coord_sf(xlim = xlim, ylim = ylim, expand = FALSE) +
      theme(panel.background = element_rect(fill = "white"), panel.border = element_rect(fill = NA), axis.text.x=element_blank(), axis.text.y=element_blank(), axis.ticks=element_blank(), axis.title = element_blank(), plot.margin = margin(t = 0.05, r = 0.05, b = 0.05, l = 0.05, unit = "pt"))
  }
  
  if(panel_or_gif == "panel"){
    # Panel plot
    all_plot<- wrap_plots(rasts_out, ncol = panel_cols, nrow = panel_rows, guides = "collect", theme(plot.margin = margin(t = 0.05, r = 0.05, b = 0.05, l = 0.05, unit = "pt")))
    ggsave(filename = paste(working_dir, file_name, ".png", sep = ""), all.plot, width = 11, height = 8, units = "in")
    } else {
      # Make a gif
      plot_loop_func<- function(plot_list){
        for (i in seq_along(plot_list)) {
          plot_use<- plot_list[[i]]
          print(plot_use)
        }
      }
      invisible(save_gif(plot_loop_func(rasts_out), paste0(out_dir, nmfs_species_code, "_LogDensity.gif"), delay = 0.75, progress = FALSE))
    }
  }

#' Predict density for new samples (\emph{Beta version; may change without notice})
#'
#' \code{predict.fit_model} calculates predictions given new data
#'
#' \code{predict.fit_model} is designed with two purposes in mind:
#' \enumerate{
#' \item If \code{new_covariate_data=NULL} as by default, then the model uses the covariate values supplied during original model fits,
#'       and interpolates as needed from those supplied values to new predicted locations.  This then uses *exactly* the same information
#'       as was available during model fitting.
#' \item If \code{new_covariate_data} is supplied with new values (e.g., at locations for predictions), then these values are used in
#'       combination with original covariate values when interpolating to new values.  However, supplying \code{new_oovariate_data}
#'       at the same Lat-Lon-Year combination as any original covariate value will delete those matches in the latter, such that originally fitted data
#'       can be predicted using alternative values for covariates (e.g., when calculating partial dependence plots)
#' }
#'
#' @inheritParams make_covariates
#' @inheritParams VAST::make_data
#' @param x Output from \code{\link{fit_model}}
#' @param what Which output from \code{fit$Report} should be extracted; default is predicted density
#' @param keep_old_covariates Whether to add new_covariate_data to existing data.
#'        This is useful when predicting values at new locations, but does not work
#'        when predicting data are locations with existing data (because the interpolation of
#'        covariate values will conflict for existing and new covariate values), e.g.,
#'        when calculating partial dependence plots for existing data.
#'
#' @return NULL
#'
#' @examples
#' \dontrun{
#'
#' # Showing use of package pdp for partial dependence plots
#' pred.fun = function( object, newdata ){
#'   predict( x=object,
#'     Lat_i = object$data_frame$Lat_i,
#'     Lon_i = object$data_frame$Lon_i,
#'     t_i = object$data_frame$t_i,
#'     a_i = object$data_frame$a_i,
#'     what = "P1_iz",
#'     new_covariate_data = newdata,
#'     do_checks = FALSE )
#' }
#'
#' library(ggplot2)
#' library(pdp)
#' Partial = partial( object = fit,
#'                    pred.var = "BOT_DEPTH",
#'                    pred.fun = pred.fun,
#'                    train = fit$covariate_data )
#' autoplot(Partial)
#'
#' }
#'
#' @method predict fit_model
#' @export

predict.fit_model_aja<- function(x, what = "D_i", Lat_i, Lon_i, t_i, a_i, c_iz = rep(0,length(t_i)), v_i = rep(0,length(t_i)), new_covariate_data = NULL, Xcontrasts_pred, X1config_cp_pred, X2config_cp_pred, new_catchability_data = NULL, do_checks = TRUE, working_dir = paste0(getwd(),"/")){
  
  if(FALSE){
    x = vast_fit
    what = predict_variable
    Lat_i = pred_lats
    Lon_i = pred_lons
    t_i = pred_times
    a_i = pred_sampled_areas
    c_iz = pred_category
    v_i = pred_vessel
    new_covariate_data = pred_cov_dat_use
    new_catchability_data = NULL
    working_dir = here::here("scratch/aja/TargetsSDM")
  }
  
  message("`predict.fit_model(.)` is in beta-testing, and please explore results carefully prior to using")
  
  # Check issues
  if( !(what%in%names(x$Report)) || (length(x$Report[[what]])!=x$data_list$n_i) ){
    stop("`what` can only take a few options")
  }
  if( !is.null(new_covariate_data) ){
    # Confirm all columns are available
    if( !all(colnames(x$covariate_data) %in% colnames(new_covariate_data)) ){
      stop("Please ensure that all columns of `x$covariate_data` are present in `new_covariate_data`")
    }
    # Eliminate unnecessary columns
    new_covariate_data = new_covariate_data[,match(colnames(x$covariate_data),colnames(new_covariate_data))]
    # Eliminate old-covariates that are also present in new_covariate_data
    NN = RANN::nn2( query=x$covariate_data[,c('Lat','Lon','Year')], data=new_covariate_data[,c('Lat','Lon','Year')], k=1 )
    if( any(NN$nn.dist==0) ){
      x$covariate_data = x$covariate_data[-which(NN$nn.dist==0),,drop=FALSE]
    }
  }
  if( !is.null(new_catchability_data) ){
    stop("Option not implemented")
  }
  
  # Process covariates
  catchability_data = rbind( x$catchability_data, new_catchability_data )
  covariate_data = rbind( x$covariate_data, new_covariate_data )
  
  # Process inputs
  PredTF_i = c( x$data_list$PredTF_i, rep(1,length(t_i)) )
  b_i = c( x$data_frame[,"b_i"], rep(1,length(t_i)) )
  c_iz = rbind( matrix(x$data_frame[,grep("c_iz",names(x$data_frame))]), matrix(c_iz) )
  Lat_i = c( x$data_frame[,"Lat_i"], Lat_i )
  Lon_i = c( x$data_frame[,"Lon_i"], Lon_i )
  a_i = c( x$data_frame[,"a_i"], a_i )
  v_i = c( x$data_frame[,"v_i"], v_i )
  t_i = c( x$data_frame[,"t_i"], t_i )
  #assign("b_i", b_i, envir=.GlobalEnv)
  
  # Build information regarding spatial location and correlation
  message("\n### Re-making spatial information")
  spatial_args_new = list("anisotropic_mesh"=x$spatial_list$MeshList$anisotropic_mesh, "Kmeans"=x$spatial_list$Kmeans,
                          "Lon_i"=Lon_i, "Lat_i"=Lat_i )
  spatial_args_input = combine_lists( input=spatial_args_new, default=x$input_args$spatial_args_input )
  spatial_list = do.call( what=make_spatial_info, args=spatial_args_input )
  
  # Check spatial_list
  if( !all.equal(spatial_list$MeshList,x$spatial_list$MeshList) ){
    stop("`MeshList` generated during `predict.fit_model` doesn't match that of original fit; please email package author to report issue")
  }
  
  # Build data
  # Do *not* restrict inputs to formalArgs(make_data) because other potential inputs are still parsed by make_data for backwards compatibility
  message("\n### Re-making data object")
  data_args_new = list( "c_iz"=c_iz, "b_i"=b_i, "a_i"=a_i, "v_i"=v_i, "PredTF_i"=PredTF_i,
                        "t_i"=t_i, "spatial_list"=spatial_list,
                        "covariate_data"=covariate_data, "catchability_data"=catchability_data )
  data_args_input = combine_lists( input=data_args_new, default=x$input_args$data_args_input )  # Do *not* use args_to_use
  data_list = do.call( what=make_data, args=data_args_input )
  data_list$n_g = 0
  
  # Build object
  message("\n### Re-making TMB object")
  model_args_default = list("TmbData"=data_list, "RunDir"=working_dir, "Version"=x$settings$Version,
                            "RhoConfig"=x$settings$RhoConfig, "loc_x"=spatial_list$loc_x, "Method"=spatial_list$Method)
  model_args_input = combine_lists( input=list("Parameters"=x$ParHat),
                                    default=model_args_default, args_to_use=formalArgs(make_model) )
  tmb_list = do.call( what=make_model, args=model_args_input )
  
  # Extract output
  Report = tmb_list$Obj$report()
  Y_i = Report[[what]][(1+nrow(x$data_frame)):length(Report$D_i)]
  
  # sanity check
  #if( all.equal(covariate_data,x$covariate_data) & Report$jnll!=x$Report$jnll){
  if( do_checks==TRUE && (Report$jnll!=x$Report$jnll) ){
    message("Problem detected in `predict.fit_model`; returning outputs for diagnostic purposes")
    Return = list("Report"=Report, "data_list"=data_list)
    return(Return)
  }
  
  # return prediction
  return(Y_i)
}