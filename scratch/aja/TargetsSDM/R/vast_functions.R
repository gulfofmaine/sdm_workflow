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
#' @param predict_covariates_stack_agg = The directory holding processed covariate raster stacks
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

make_vast_predict_df<- function(predict_covariates_stack_agg, extra_covariates_stack, covs_rescale = c("Depth", "BS_seasonal", "BT_seasonal", "SS_seasonal", "SST_seasonal"), rescale_params, mask, summarize, ensemble_stat, fit_seasons, fit_year_min, fit_year_max, pred_years, out_dir){
  
  # For debugging
  if(FALSE){
    tar_load(predict_covariates_stack_agg_out)
    predict_covariates_stack_agg<- predict_covariates_stack_agg_out
    tar_load(static_covariates_stack)
    extra_covariates_stack = static_covariates_stack
    tar_load(rescale_params)
    tar_load(region_shapefile)
    mask = region_shapefile
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
  rast_files_load<- list.files(predict_covariates_stack_agg, pattern = paste0(summarize, "_", ensemble_stat, ".grd"), full.names = TRUE)
  
  # Get variable names
  cov_names_full<- list.files(predict_covariates_stack_agg, pattern = paste0(summarize, "_", ensemble_stat, ".grd"), full.names = FALSE)
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
  
  # New implementation...
  pred_covs_out_final<- pred_covs_out_final %>%
    mutate(., VAST_YEAR_COV = ifelse(EST_YEAR > fit_year_max, fit_year_max, EST_YEAR),
           VAST_SEASON = case_when(
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
  
  # Drop any NAs
  pred_covs_out_final<- pred_covs_out_final %>%
    drop_na(., {{cov_names}}) %>%
    mutate(., "Summarized" = summarize,
           "Ensemble_Stat" = ensemble_stat)
  
  # Rescale
  if(!is.null(rescale_params)){
    for(i in seq_along(covs_rescale)){
      match_mean<- rescale_params[which(names(rescale_params) == paste(covs_rescale[i], "Mean", sep = "_"))]
      match_sd<- rescale_params[which(names(rescale_params) == paste(covs_rescale[i], "SD", sep = "_"))]
      pred_covs_out_final<- pred_covs_out_final %>%
        mutate_at(., {{covs_rescale[i]}}, .funs = covariate_rescale_func, type = "AJA", center = match_mean, scale = match_sd)
    }
  }
  
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
  # Previous implementation before trying to include both surveys within a given season
  # data_temp<- tidy_mod_data %>%
  #   filter(., NMFS_SVSPP == nmfs_species_code) %>%
  #   filter(., EST_YEAR >= fit_year_min & EST_YEAR <= fit_year_max) %>%
  #   mutate(., "VAST_SEASON" = case_when(
  #     SURVEY == "DFO" & SEASON == "SPRING" ~ "DFO",
  #     SURVEY == "NMFS" & SEASON == "SPRING" ~ "SPRING",
  #     SURVEY == "DFO" & SEASON == "SUMMER" ~ "SUMMER",
  #     SURVEY == "NMFS" & SEASON == "FALL" ~ "FALL")) %>%
  #   drop_na(VAST_SEASON)
  
  # New implementatiom...
  data_temp<- tidy_mod_data %>%
    filter(., NMFS_SVSPP == nmfs_species_code) %>%
    filter(., EST_YEAR >= fit_year_min & EST_YEAR <= fit_year_max) %>%
    mutate(., "VAST_SEASON" = case_when(
      SURVEY == "DFO" & SEASON == "SPRING" ~ "SPRING",
      SURVEY == "NMFS" & SEASON == "SPRING" ~ "SPRING",
      SURVEY == "DFO" & SEASON == "SUMMER" ~ "SUMMER",
      SURVEY == "NMFS" & SEASON == "FALL" ~ "FALL",
      SURVEY == "DFO" & SEASON == "FALL" ~ "FALL")) %>%
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
  
  # Select columns we want from the "full" vast_seasonal_data dataset. Area swept Marine fish diversity on the Scotian Shelf, Canada
  vast_samp_dat<- data.frame(
    "Year" = as.numeric(vast_seasonal_data$VAST_YEAR_SEASON)-1, 
    "Lat" = vast_seasonal_data$DECDEG_BEGLAT,
    "Lon" = vast_seasonal_data$DECDEG_BEGLON,
    "Biomass" = vast_seasonal_data$BIOMASS,
    "Swept" = ifelse(vast_seasonal_data$SURVEY == "NMFS", 0.0384, 0.0404),
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
#' @param rescale = Logical indicating whether or not the covariates should be rescaled. 
#' @param out_dir = Description
#' 
#' @return A sample dataframe that includes all of the covariate information at each unique sample. This file is also saved in out_dir. 
#' 
#' @export

make_vast_covariate_data<- function(vast_seasonal_data, out_dir){
  
  # For debugging
  if(FALSE){
    tar_load(vast_seasonal_data)
    rescale = 
    out_dir = here::here("scratch/aja/targets_flow/data/dfo/combined")
  }
  
  # Select columns we want from the "full" vast_seasonal_data dataset
  vast_cov_dat<- data.frame(
    "Year" = as.numeric(vast_seasonal_data$VAST_YEAR_SEASON)-1,
    "Year_Cov" = vast_seasonal_data$VAST_YEAR_COV,
    "Season" = vast_seasonal_data$VAST_SEASON,
    "Depth" = vast_seasonal_data$Depth,
    "SST_seasonal" = vast_seasonal_data$SST_seasonal,
    "BT_seasonal" = vast_seasonal_data$BT_seasonal,
    "Lat" = vast_seasonal_data$DECDEG_BEGLAT,
    "Lon" = vast_seasonal_data$DECDEG_BEGLON
  )
  
  # Save and return 
  saveRDS(vast_cov_dat, file = paste(out_dir, "vast_covariate_data.rds", sep = "/"))
  return(vast_cov_dat)
} 

#' @title Make VAST catachability
#' 
#' @description This function creates a VAST catachability dataset to pass into calls to `VAST::fit_model`.
#'
#' @param vast_seasonal_data = Description 
#' @param out_dir = Description
#' 
#' @return A sample dataframe that includes all of the covariate information at each unique sample. This file is also saved in out_dir. 
#' 
#' @export

make_vast_catchability_data<- function(vast_seasonal_data, out_dir){
  
  # For debugging
  if(FALSE){
    vast_seasonal_data
    out_dir = here::here("scratch/aja/targets_flow/data/dfo/combined")
  }
  
  # Select columns we want from the "full" vast_seasonal_data dataset
  vast_catch_dat<- data.frame(
    "Year" = as.numeric(vast_seasonal_data$VAST_YEAR_SEASON)-1,
    "Year_Cov" = vast_seasonal_data$VAST_YEAR_COV,
    "Season" = vast_seasonal_data$VAST_SEASON,
    "Lat" = vast_seasonal_data$DECDEG_BEGLAT,
    "Lon" = vast_seasonal_data$DECDEG_BEGLON,
    "Survey" = factor(vast_seasonal_data$SURVEY, levels = c("NMFS", "DFO", "DUMMY"))
  )
  
  # Save and return it
  saveRDS(vast_catch_dat, file = paste(out_dir, "vast_catchability_data.rds", sep = "/"))
  return(vast_catch_dat)
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

vast_make_extrap_grid<- function(region_shapefile, index_shapes, strata.limits, cell_size){
  
  # For debugging
  if(FALSE){
    tar_load(region_shapefile)
    tar_load(index_shapefiles)
    index_shapes = index_shapefiles
    strata.limits = strata_use
    cell_size = 25000
  }
  
  # Transform crs of shapefile to common WGS84 lon/lat format.
  region_wgs84<- st_transform(region_shapefile, crs = "+proj=longlat +lat_0=90 +lon_0=180 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0 ")
  
  # Get UTM zone
  lon<- sum(st_bbox(region_wgs84)[c(1,3)])/2
  utm_zone<- floor((lon + 180)/6)+1
  
  # Transform to the UTM zone
  crs_utm<- st_crs(paste0("+proj=utm +zone=",utm_zone," +ellps=WGS84 +datum=WGS84 +units=m +no_defs "))
  region_utm<- st_transform(region_wgs84, crs = crs_utm)
  
  # Make extrapolation grid with sf
  region_grid<- st_as_sf(st_make_grid(region_utm, cellsize = cell_size, what = "centers"), crs = crs_utm) 
  
  # Now get only the points that fall within the shape polygon
  points_keep<- data.frame("pt_row" = seq(from = 1, to = nrow(region_grid), by = 1), "in_out" = st_intersects(region_grid, region_utm, sparse = FALSE))            
  region_grid<- region_grid %>%
    mutate(., "in_poly" = st_intersects(region_grid, region_utm, sparse = FALSE)) %>%
    filter(., in_poly == TRUE)
  
  # Convert back to WGS84 lon/lat, as that is what VAST expects.
  extrap_grid<- region_grid %>%
    st_transform(., crs = "+proj=longlat +lat_0=90 +lon_0=180 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0 ") %>%
    st_join(., index_shapes, join = st_within) %>%
    mutate(., "Lon" = st_coordinates(.)[,1],
           "Lat" = st_coordinates(.)[,2]) %>%
    st_drop_geometry() %>%
    dplyr::select(., Lon, Lat, Region) %>%
    mutate(., Area_km2=((cell_size/1000)^2),
           STRATA = factor(Region, levels = index_shapes$Region, labels = index_shapes$Region))
    
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
#' @param Options = Tagged vector to turn on or off specific options (e.g., SD_site_logdensity, Effective area, etc)
#' @param strata.limits
#'
#' @return Tagged list containing settings needed to fit a VAST model of species occurrence.
#' 
#' @export

vast_make_settings <- function(extrap_grid, FieldConfig, RhoConfig, OverdispersionConfig, bias.correct, Options, strata.limits){
  
  # For debugging
  if(FALSE){
    tar_load(vast_extrap_grid)
    extrap_grid = vast_extrap_grid
    FieldConfig = c("Omega1" = 1, "Epsilon1" = 1, "Omega2" = 1, "Epsilon2" = 1)
    RhoConfig = c("Beta1" = 3, "Beta2" = 3, "Epsilon1" = 2, "Epsilon2" = 2)
    OverdispersionConfig = c(0, 0)
    bias.correct = FALSE
    Options = NULL
    strata.limits = strata_use
  }
  
  # Get number of vertices in the mesh, which is based on our input_grid
  n_x_use<- nrow(extrap_grid)
  
  # Run FishStatsUtils::make_settings
  settings_out<- make_settings(n_x = n_x_use, Region = "User", purpose = "index2", FieldConfig = FieldConfig, RhoConfig = RhoConfig, ObsModel = c(1, 1), OverdispersionConfig = OverdispersionConfig, bias.correct = bias.correct, knot_method = "grid", treat_nonencounter_as_zero = FALSE, strata.limits = strata.limits)
  
  # Adjust options?
  options_new<- settings_out$Options
  if(!is.null(Options)){
    for(i in seq_along(Options)){
      options_adjust_i<- Options[[i]]
      options_new[[{{options_adjust_i}}]]<- Options[[{{options_adjust_i}}]]
    }
    settings_out<- make_settings(n_x = n_x_use, Region = "User", purpose = "index2", FieldConfig = FieldConfig, RhoConfig = RhoConfig, ObsModel = c(1, 1), OverdispersionConfig = OverdispersionConfig, bias.correct = bias.correct, knot_method = "grid", treat_nonencounter_as_zero = FALSE, Options = options_new)
  }
  
  # Return it
  return(settings_out)
}

####
#' @title Make VAST covariate effect objects
#' 
#' @description Create covariate effects for both linear predictors
#'
#' @param X1_coveff_vec = A vector specifying the habitat covariate effects for first linear predictor. 
#' @param X2_coveff_vec = A vector specifying the habitat covariate effects for second linear predictor.
#' @param Q1_coveff_vec = A vector specifying the catchability covariate effects for first linear predictor. 
#' @param Q2_coveff_vec = A vector specifying the catchability covariate effects for second linear predictor. 
#'
#' @return A list with covariate effects for the habitat covariates and first linear predictor (first list slot), habitat covariates and second linear predictor (second list slot), catchability covariates and first linear predictor (third slot) and catchability covariates and second linear predictor (fourth slot).  
#'
#' @export

vast_make_coveff<- function(X1_coveff_vec, X2_coveff_vec, Q1_coveff_vec, Q2_coveff_vec){
  
  # For debugging
  if(FALSE){
   X1_coveff_vec = c(2, 3, 3, 2, rep(3, 32))
   X2_coveff_vec = c(2, 3, 3, 2, rep(3, 32))
   Q1_coveff_vec = NULL 
   Q2_coveff_vec = NULL
  }
  
  # Combine into a list and name it 
  if(is.null(Q1_coveff_vec) | is.null(Q2_coveff_vec)){
    coveff_out<- list("X1config_cp" = matrix(X1_coveff_vec, nrow = 1), "X2config_cp" = matrix(X2_coveff_vec, nrow = 1), "Q1config_k" = NULL, "Q2config_k" = NULL)
  } else {
    coveff_out<- list("X1config_cp" = matrix(X1_coveff_vec, nrow = 1), "X2config_cp" = matrix(X2_coveff_vec, nrow = 1), "Q1config_k" = matrix(Q1_coveff_vec, nrow = 1), "Q2config_k" = matrix(Q2_coveff_vec, nrow = 1))
  }
  
  
  # Return it
  return(coveff_out)
}

####
#' @title Build VAST SDM
#' 
#' @description Build VAST species distribution model, without running it. This can be helpful to check settings before running `vast_fit_sdm`. Additionally, it can be helpful for making subsequent modifications, particularly to mapping.
#'
#' @param settings = A tagged list with the settings for the model, created with `vast_make_settings`.
#' @param extrap_grid = An extrapolation grid, created with `vast_make_extrap_grid`.
#' @param sample_dat = A data frame with the biomass sample data for each species at each tow.
#' @param covariate_dat = A data frame with the habitat covariate data for each tow.
#' @param X1_formula = A formula for the habitat covariates and first linear predictor.
#' @param X2_formula = A formula for the habitat covariates and second linear predictor.
#' @param X_contrasts = A tagged list specifying the contrasts to use for factor covariates in the model.
#' @param Xconfig_list = A tagged list specifying the habitat and catchability covariate effects for first and second linear predictors.
#' @param catchability_data = A data frame with the catchability data for every sample
#' @param Q1_formula = A formula for the catchability covariates and first linear predictor.
#' @param Q2_formula = A formula for the catchability covariates and second linear predictor.
#' @param index_shapefiles = A sf object with rows for each of the regions of interest
#'
#' @return A VAST `fit_model` object, with the inputs and built TMB object components.
#'
#' @export

vast_build_sdm <- function(settings, extrap_grid, sample_data, covariate_data, X1_formula, X2_formula, X_contrasts, Xconfig_list, catchability_data, Q1_formula, Q2_formula, index_shapes){
  
  # For debugging
  if(FALSE){
    library(VAST)
    library(tidyverse)
    library(stringr)
    
    # Seasonal
    tar_load(vast_settings)
    settings = vast_settings
    tar_load(vast_extrap_grid)
    extrap_grid = vast_extrap_grid
    tar_load(vast_sample_data)
    sample_data = vast_sample_data
    tar_load(vast_covariate_data)
    covariate_data = vast_covariate_data
    X1_formula = hab_formula
    X2_formula = hab_formula
    hab_env_coeffs_n = hab_env_coeffs_n
    tar_load(vast_catchability_data)
    catchability_data = vast_catchability_data
    catch_formula<- ~ Survey
    Q1_formula = catch_formula
    Q2_formula = catch_formula
    X_contrasts = list(Season = contrasts(vast_covariate_data$Season, contrasts = FALSE), Year_Cov = contrasts(vast_covariate_data$Year_Cov, contrasts = FALSE))
    # X_contrasts = list(Year_Cov = contrasts(vast_covariate_data$Year_Cov, contrasts = FALSE))
    tar_load(vast_coveff)
    Xconfig_list = vast_coveff
    tar_load(index_shapefiles)
    
    # Annual
    tar_load(vast_settings)
    settings = vast_settings
    tar_load(vast_extrap_grid)
    extrap_grid = vast_extrap_grid
    tar_load(vast_sample_data)
    sample_data = vast_sample_data
    tar_load(vast_covariate_data)
    covariate_data = vast_covariate_data
    X1_formula = hab_formula
    X2_formula = hab_formula
    hab_env_coeffs_n = hab_env_coeffs_n
    tar_load(vast_catchability_data)
    catchability_data = vast_catchability_data
    catch_formula<- ~ Survey
    Q1_formula = catch_formula
    Q2_formula = catch_formula
    X_contrasts = list(Year_Cov = contrasts(vast_covariate_data$Year_Cov, contrasts = FALSE))
    tar_load(vast_coveff)
    Xconfig_list = vast_coveff
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
  
  if(!(all(c("X1config_cp", "X2config_cp", "Q1config_k", "Q2config_k") %in% names(Xconfig_list)))){
    stop(paste("Check names of Xconfig_list. Must be", paste0(c("X1config_cp", "X2config_cp", "Q1config_k", "Q2config_k"), collapse = ","), sep = ""))
  }
  
  # Run VAST::fit_model with correct info and settings
  vast_build_out<- fit_model_aja("settings" = settings, "input_grid" = extrap_grid, "Lat_i" = sample_data[, 'Lat'], "Lon_i" = sample_data[, 'Lon'], "t_i" = sample_data[, 'Year'], "c_i" = rep(0, nrow(sample_data)), "b_i" = sample_data[, 'Biomass'], "a_i" = sample_data[, 'Swept'], "PredTF_i" = sample_data[, 'Pred_TF'], "X1config_cp" = Xconfig_list[['X1config_cp']], "X2config_cp" = Xconfig_list[['X2config_cp']], "covariate_data" = covariate_data, "X1_formula" = X1_formula, "X2_formula" = X2_formula, "X_contrasts" = X_contrasts, "catchability_data" = catchability_data, "Q1_formula" = Q1_formula, "Q2_formula" = Q2_formula, "Q1config_k" = Xconfig_list[['Q1config_k']], "Q2config_k" = Xconfig_list[['Q2config_k']], "newtonsteps" = 1, "getsd" = TRUE, "getReportCovariance" = TRUE, "run_model" = FALSE, "test_fit" = FALSE,  "Use_REML" = FALSE, "getJointPrecision" = FALSE, "index_shapes" = index_shapes)
  
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
#' @param index_shapefiles = A sf object with rows for each of the regions of interest
#'
#' @return A VAST fit_model object, with the inputs and built TMB object components.
#' 
#' @export

vast_make_adjustments <- function(vast_build, index_shapes, adjustments = NULL){
  
  # For debugging
  if(FALSE){
    tar_load(vast_build0)
    vast_build = vast_build0
    tar_load(vast_covariate_data)
    adjustments = list("log_sigmaXi1_cp" = factor(c(rep(1, 3), rep(4, nlevels(vast_covariate_data$Year_Cov)), rep(NA, 2))), "log_sigmaXi2_cp" = factor(c(rep(1, 3), rep(4, nlevels(vast_covariate_data$Year_Cov)), rep(NA, 2))))
    adjustments = list("log_sigmaXi1_cp" = factor(c(rep(4, nlevels(vast_covariate_data$Year_Cov)), rep(NA, 2*hab_env_coeffs_n))), "log_sigmaXi2_cp" = factor(c(rep(4, nlevels(vast_covariate_data$Year_Cov)), rep(NA, 2*hab_env_coeffs_n))), "lambda1_k" = factor(c(NA, 1)), "lambda2_k" = factor(c(NA, 1)), "RhoConfig" = c("Beta1" = 3, "Beta2" = 3, "Epsilon1" = 3, "Epsilon2" = 3))
  }
  
  # If no adjustments are needed, just need to pull information from vast_build and then set "run_model" to TRUE
  if(is.null(adjustments)){
    vast_build_adjust_out<- fit_model_aja("settings" = vast_build$settings, "input_grid" = vast_build$input_args$data_args_input$input_grid, "Lat_i" = vast_build$data_frame[, 'Lat_i'], "Lon_i" = vast_build$data_frame[, 'Lon_i'], "t_i" = vast_build$data_frame[, 't_i'], "c_iz" = vast_build$data_frame[, 'c_iz'], "b_i" = vast_build$data_frame[, 'b_i'], "a_i" = vast_build$data_frame[, 'a_i'], "PredTF_i" = vast_build$data_list[['PredTF_i']], "X1config_cp" = vast_build$input_args$data_args_input[['X1config_cp']], "X2config_cp" = vast_build$input_args$data_args_input[['X2config_cp']], "covariate_data" = vast_build$input_args$data_args_input$covariate_data, "X1_formula" = vast_build$input_args$data_args_input$X1_formula, "X2_formula" = vast_build$input_args$data_args_input$X2_formula, "X_contrasts" = vast_build$input_args$data_args_input$X_contrasts, "catchability_data" = vast_build$input_args$data_args_input$catchability_data, "Q1_formula" = vast_build$input_args$data_args_input$Q1_formula, "Q2_formula" = vast_build$input_args$data_args_input$Q2_formula, "Q1config_k" = vast_build$input_args$data_args_input[['Q1config_cp']], "Q2config_k" = vast_build$input_args$data_args_input[['Q2config_k']], "newtonsteps" = 1, "getsd" = TRUE, "getReportCovariance" = TRUE, "run_model" = FALSE, "test_fit" = FALSE,  "Use_REML" = FALSE, "getJointPrecision" = FALSE, "index_shapes" = index_shapes)
  }
  
  # If there are adjustments, need to make those and then re run model. 
  if(!is.null(adjustments)){
    # Check names -- trying to think of what the possible adjustment flags would be in the named list
    adjust_names<- c("FieldConfig", "RhoConfig", "X1_formula", "X2_formula", "X1config_cp", "X2config_cp", "X_contrasts", "log_sigmaXi1_cp", "log_sigmaXi2_cp", "lambda1_k", "lambda2_k", "Q1_formula", "Q2_formula", "Q1config_k", "Q2config_k")
    
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
    if(any(names(adjustments) %in% c("log_sigmaXi1_cp", "log_sigmaXi2_cp", "lambda1_k", "lambda2_k"))){
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
    if(any(names(adjustments) %in% c("log_sigmaXi1_cp", "log_sigmaXi2_cp", "lambda1_k", "lambda2_k"))){
      # Adding Map argument
      vast_build_adjust_out<- fit_model_aja("settings" = vast_build$settings, "input_grid" = vast_build$input_args$data_args_input$input_grid, "Lat_i" = vast_build$data_frame[, 'Lat_i'], "Lon_i" = vast_build$data_frame[, 'Lon_i'], "t_i" = vast_build$data_frame[, 't_i'], "c_iz" = vast_build$data_frame[, 'c_iz'], "b_i" = vast_build$data_frame[, 'b_i'], "a_i" = vast_build$data_frame[, 'a_i'], "PredTF_i" = vast_build$data_list[['PredTF_i']], "X1config_cp" = vast_build$input_args$data_args_input[['X1config_cp']], "X2config_cp" = vast_build$input_args$data_args_input[['X2config_cp']], "covariate_data" = vast_build$input_args$data_args_input$covariate_data, "X1_formula" = vast_build$input_args$data_args_input$X1_formula, "X2_formula" = vast_build$input_args$data_args_input$X2_formula, "X_contrasts" = vast_build$input_args$data_args_input$X_contrasts, "catchability_data" = vast_build$input_args$data_args_input$catchability_data, "Q1_formula" = vast_build$input_args$data_args_input$Q1_formula, "Q2_formula" = vast_build$input_args$data_args_input$Q2_formula, "Q1config_k" = vast_build$input_args$data_args_input[['Q1config_k']], "Q2config_k" = vast_build$input_args$data_args_input[['Q2config_k']], "Map" = map_adjust_out, "newtonsteps" = 1, "getsd" = TRUE, "getReportCovariance" = TRUE, "run_model" = FALSE, "test_fit" = FALSE,  "Use_REML" = FALSE, "getJointPrecision" = FALSE, "index_shapes" = index_shapes)
    } else {
      # No need for Map argument, just build and fit
      vast_build_adjust_out<- fit_model_aja("settings" = vast_build$settings, "input_grid" = vast_build$input_args$data_args_input$input_grid, "Lat_i" = vast_build$data_frame[, 'Lat_i'], "Lon_i" = vast_build$data_frame[, 'Lon_i'], "t_i" = vast_build$data_frame[, 't_i'], "c_iz" = vast_build$data_frame[, 'c_iz'], "b_i" = vast_build$data_frame[, 'b_i'], "a_i" = vast_build$data_frame[, 'a_i'], "PredTF_i" = vast_build$data_list[['PredTF_i']], "X1config_cp" = vast_build$input_args$data_args_input[['X1config_cp']], "X2config_cp" = vast_build$input_args$data_args_input[['X2config_cp']], "covariate_data" = vast_build$input_args$data_args_input$covariate_data, "X1_formula" = vast_build$input_args$data_args_input$X1_formula, "X2_formula" = vast_build$input_args$data_args_input$X2_formula, "X_contrasts" = vast_build$input_args$data_args_input$X_contrasts, "catchability_data" = vast_build$input_args$data_args_input$catchability_data, "Q1_formula" = vast_build$input_args$data_args_input$Q1_formula, "Q2_formula" = vast_build$input_args$data_args_input$Q2_formula, "Q1config_cp" = vast_build$input_args$data_args_input[['Q1config_cp']], "Q2config_cp" = vast_build$input_args$data_args_input[['Q2config_cp']], "newtonsteps" = 1, "getsd" = TRUE, "getReportCovariance" = TRUE, "run_model" = FALSE, "test_fit" = FALSE,  "Use_REML" = FALSE, "getJointPrecision" = FALSE, "index_shapes" = index_shapes)
    }
  }
  # Return it
  return(vast_build_adjust_out)
}

#' @title Fit VAST SDM
#' 
#' @description Fit VAST species distribution model
#'
#' @param vast_build_adjust = A VAST `fit_model` object.
#' @param nmfs_species_code
#' @param index_shapefiles = A sf object with rows for each of the regions of interest
#' @param out_dir
#'
#' @return A VAST fit_model object, with the inputs and and outputs, including parameter estimates, extrapolation gid info, spatial list info, data info, and TMB info.
#'
#' @export

vast_fit_sdm <- function(vast_build_adjust, nmfs_species_code, index_shapes, out_dir){
  
  # For debugging
  if(FALSE){
    tar_load(vast_adjust)
    vast_build_adjust = vast_adjust
  }
  
  # Build and fit model
  vast_fit_out<- fit_model_aja("settings" = vast_build_adjust$settings, "input_grid" = vast_build_adjust$input_args$data_args_input$input_grid, "Lat_i" = vast_build_adjust$data_frame[, 'Lat_i'], "Lon_i" = vast_build_adjust$data_frame[, 'Lon_i'], "t_i" = vast_build_adjust$data_frame[, 't_i'], "c_iz" = vast_build_adjust$data_frame[, 'c_iz'], "b_i" = vast_build_adjust$data_frame[, 'b_i'], "a_i" = vast_build_adjust$data_frame[, 'a_i'], "PredTF_i" = vast_build_adjust$data_list[['PredTF_i']], "X1config_cp" = vast_build_adjust$input_args$data_args_input[['X1config_cp']], "X2config_cp" = vast_build_adjust$input_args$data_args_input[['X2config_cp']], "covariate_data" = vast_build_adjust$input_args$data_args_input$covariate_data, "X1_formula" = vast_build_adjust$input_args$data_args_input$X1_formula, "X2_formula" = vast_build_adjust$input_args$data_args_input$X2_formula, "X_contrasts" = vast_build_adjust$input_args$data_args_input$X_contrasts, "catchability_data" = vast_build_adjust$input_args$data_args_input$catchability_data, "Q1_formula" = vast_build_adjust$input_args$data_args_input$Q1_formula, "Q2_formula" = vast_build_adjust$input_args$data_args_input$Q2_formula, "Q1config_cp" = vast_build_adjust$input_args$data_args_input[['Q1config_cp']], "Q2config_cp" = vast_build_adjust$input_args$data_args_input[['Q2config_cp']], "Map" = vast_build_adjust$tmb_list$Map, "newtonsteps" = 1, "getsd" = TRUE, "getReportCovariance" = TRUE, "run_model" = TRUE, "test_fit" = FALSE,  "Use_REML" = FALSE, "getJointPrecision" = FALSE, "index_shapes" = index_shapes)
  
  # Save and return it
  saveRDS(vast_fit_out, file = paste(out_dir, "/", nmfs_species_code, "_", "fitted_vast.rds", sep = "" ))
  return(vast_fit_out)
}

#' @title Predict fitted VAST model
#'
#' @description This function makes predictions from a fitted VAST SDM to new locations using VAST::predict.fit_model. Importantly, to use this feature for new times, at least one location for each time of interest needs to be included during the model fitting process. This dummy observation should have a PredTF value of 1 so that the observation is only used in the predicted probability and NOT estimating the likelihood.
#'
#' @param vast_fitted_sdm = A fitted VAST SDM object, as returned with `vast_fit_sdm`
#' @param nmfs_species_code = A numeric species code
#' @param predict_variable = Which variable should be predicted, default is density (D_i)
#' @param predict_category = Which category (species/age/size) should be predicted, default is 0
#' @param predict_vessel = Which sampling category should be predicted, default is 0
#' @param predict_covariates_df_all = A long data frame with all of the prediction covariates
#' @param memory_save = Logical. If TRUE, then predictions are only made to knots as defined within the vast_fitted_sdm object. This is done by finding the prediction locations that are nearest neighbors to each knot. If FALSE, then predictions are made to each of the locations in the predict_covariates_df_all.  
#' @param out_dir = Output directory to save...
#' 
#' @return
#'
#' @export

predict_vast<- function(vast_fitted_sdm, nmfs_species_code, predict_variable = "D_i", predict_category = 0, predict_vessel = 0, predict_covariates_df_all, out_dir){

  # For debugging
  if(FALSE){
    tar_load(vast_fit)
    vast_fitted_sdm = vast_fit
    nmfs_species_code = nmfs_species_code
    predict_variable = "D_i"
    predict_category = 0
    predict_vessel = 0
    tar_load(vast_predict_df)
    predict_covariates_df_all = vast_predict_df
  }
  
  #### Not the biggest fan of this, but for now, building in a work around to resolve some of the memory issues that we were running into by supplying a 0.25 degree grid and trying to predict/project for each season-year from 1980-2100. To overcome this issue, going to try to just make the projections to knots and do the smoothing later.
  # First, need to get the knot locations
  knot_locs<- data.frame(vast_fitted_sdm$spatial_list$latlon_g) %>%
    st_as_sf(., coords = c("Lat", "Lon"), remove = FALSE) %>%
    mutate(., "Pt_Id" = 1:nrow(.))
  
  # Nearest knot to each point?
  pred_sf<- predict_covariates_df_all %>%
    st_as_sf(., coords = c("DECDEG_BEGLAT", "DECDEG_BEGLON"), remove = FALSE)
  
  pred_sf<- pred_sf %>%
    mutate(., "Nearest_Knot" = st_nearest_feature(., knot_locs))
  
  # Average the points...
  pred_sf_knots<- pred_sf %>%
    group_by(., ID, DATE, EST_YEAR, SEASON, SURVEY, SVVESSEL, NMFS_SVSPP, DFO_SPEC, PRESENCE, BIOMASS, ABUNDANCE, PredTF, VAST_YEAR_COV, VAST_SEASON, VAST_YEAR_SEASON, Summarized, Ensemble_Stat, Nearest_Knot) %>%
    summarize(., "BS_seasonal" = mean(BS_seasonal, na.rm = TRUE),
              "BT_seasonal" = mean(BT_seasonal, na.rm = TRUE),
              "SS_seasonal" = mean(SS_seasonal, na.rm = TRUE),
              "SST_seasonal" = mean(SST_seasonal, na.rm = TRUE),
              "Depth" = mean(Depth, na.rm = TRUE)) %>%
    st_drop_geometry() %>%
    left_join(., st_drop_geometry(knot_locs), by = c("Nearest_Knot" = "Pt_Id")) %>%
    ungroup()
  
  # Collecting necessary bits from the prediction covariates -- lat, lon, time
  pred_lats<- pred_sf_knots$Lat
  pred_lons<- pred_sf_knots$Lon
  pred_times<- as.numeric(pred_sf_knots$VAST_YEAR_SEASON)-1
  
  # Catch stuff...
  pred_sampled_areas<- rep(1, length(pred_lats))
  pred_category<- rep(predict_category, length(pred_lats))
  pred_vessel<- rep(predict_vessel, length(pred_lats))
  
  # Renaming predict_covariates_df_all to match vast_fit_covariate_data
  pred_cov_dat_use<- data.frame(
    "Year" = pred_times,
    "Year_Cov" = pred_sf_knots$VAST_YEAR_COV,
    "Season" = pred_sf_knots$VAST_SEASON,
    "Depth" = pred_sf_knots$Depth,
    "SST_seasonal" = pred_sf_knots$SST_seasonal,
    "BT_seasonal" = pred_sf_knots$BT_seasonal,
    "Survey" = pred_vessel,
    "Lat" = pred_lats,
    "Lon" = pred_lons
  )
 
  pred_catch_dat_use<- pred_cov_dat_use %>%
    dplyr::select(., c(Year, Year_Cov, Season, Lat, Lon, Survey)
  )
  pred_catch_dat_use$Survey<- rep("NMFS", nrow(pred_catch_dat_use))
  pred_catch_dat_use$Survey<- factor(pred_catch_dat_use$Survey, levels = c("NMFS", "DFO", "DUMMY"))

  # Make the predictions
  preds_out<- predict.fit_model_aja(x = vast_fitted_sdm, what = predict_variable, Lat_i = pred_lats, Lon_i = pred_lons, t_i = pred_times, a_i = pred_sampled_areas, c_iz = pred_category, NULL, new_covariate_data = pred_cov_dat_use, new_catchability_data = pred_catch_dat_use, do_checks = FALSE)
  
  # Get everything as a dataframe to make plotting easier...
  pred_df_out<- data.frame("Lat" = pred_lats, "Lon" = pred_lons, "Time" = pred_sf_knots$VAST_YEAR_SEASON, "Pred" = preds_out)
  
  # Save and return it
  saveRDS(pred_df_out, file = paste(out_dir, "/pred_Di_", nmfs_species_code, "_", unique(pred_sf_knots$Summarized), "_", unique(pred_sf_knots$Ensemble_Stat), ".rds", sep = "" ))
  return(pred_df_out)
}

#' @title Prediction spatial summary 
#' 
#' @description Calculates average "availability" of fish biomass from SDM predictions within spatial area of interest
#'
#' @param pred_df = A dataframe with Lat, Lon, Time and Pred columns
#' @param spatial_areas = 
#' @return What does this function return?
#'
#' @export

pred_spatial_summary<- function(pred_df, spatial_areas){
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

vast_fit_plot_density<- function(vast_fit, nice_category_names, mask, all_times = all_times, plot_times = NULL, land_sf, xlim, ylim, panel_or_gif = "gif", out_dir, land_color = "#d9d9d9", panel_cols = NULL, panel_rows = NULL, ...){
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
    
    ggsave(filename = paste(out_dir, file_name, ".png", sep = "/"), plot_out, width = 11, height = 8, units = "in")
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
      ggsave(filename = paste(working_dir, file_name, ".png", sep = "/"), all.plot, width = 11, height = 8, units = "in")
    } else {
      # Make a gif
      plot_loop_func<- function(plot_list){
        for (i in seq_along(plot_list)) {
          plot_use<- plot_list[[i]]
          print(plot_use)
        }
      }
      invisible(save_gif(plot_loop_func(rasts_out), paste0(out_dir, "/", nice_category_names, "_LogDensity.gif"), delay = 0.75, progress = FALSE))
    }
  }
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

predict.fit_model_aja<- function(x, what = "D_i", Lat_i, Lon_i, t_i, a_i, c_iz = rep(0,length(t_i)), v_i = rep(0,length(t_i)), new_covariate_data = NULL, new_catchability_data = NULL, do_checks = TRUE, working_dir = paste0(getwd(),"/")){
  
  if(FALSE){
    tar_load(vast_fit)
    x = vast_fit
    what = "D_i"
    Lat_i = x$data_frame$Lat_i
    #Lat_i = pred_cov_dat_use$Lat
    Lon_i = x$data_frame$Lon_i
    #Lon_i = pred_cov_dat_use$Lon
    t_i = x$data_frame$t_i
    #t_i = pred_cov_dat_use$Year
    a_i<- x$data_frame$a_i
    #a_i<- rep(unique(pred_sampled_areas), length(Lat_i))
    c_iz = rep(0,length(t_i))
    #c_iz<- rep(unique(predict_category), length(Lat_i))
    v_i = rep(0,length(t_i))
    #v_i<- rep(unique(predict_vessel), length(t_i))
    new_covariate_data = NULL
    #new_covariate_data = pred_cov_dat_use
    new_catchability_data = NULL
    #new_catchability_data = pred_catch_dat_use
    do_checks = FALSE
    
    x = vast_fitted_sdm
    what = predict_variable
    Lat_i = pred_lats
    Lon_i = pred_lons
    t_i = pred_times
    a_i = pred_sampled_areas
    c_iz = pred_category
    v_i = pred_vessel
    new_covariate_data = pred_cov_dat_use
    new_catchability_data = pred_catch_dat_use
    do_checks = FALSE
    
    object = vast_fit
    x = object
    Lat_i = object$data_frame$Lat_i
    Lon_i = object$data_frame$Lon_i
    t_i = object$data_frame$t_i
    a_i = object$data_frame$a_i
    c_iz = rep(0,length(t_i))
    v_i = rep(0,length(t_i))
    what = "P1_iz"
    new_covariate_data = object$covariate_data
    new_catchability_data = object$catchability_data
    do_checks = FALSE
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
    # Confirm all columns are available
    if( !all(colnames(x$catchability_data) %in% colnames(new_catchability_data)) ){
      stop("Please ensure that all columns of `x$catchability_data` are present in `new_covariate_data`")
    }
    # Eliminate unnecessary columns
    new_catchability_data = new_catchability_data[,match(colnames(x$catchability_data),colnames(new_catchability_data))]
    # Eliminate old-covariates that are also present in new_covariate_data
    NN = RANN::nn2( query=x$catchability_data[,c('Lat','Lon','Year')], data=new_catchability_data[,c('Lat','Lon','Year')], k=1 )
    if( any(NN$nn.dist==0) ){
      x$catchability_data = x$catchability_data[-which(NN$nn.dist==0),,drop=FALSE]
    }
  }

  # Process covariates
  covariate_data = rbind( x$covariate_data, new_covariate_data )
  catchability_data = rbind( x$catchability_data, new_catchability_data )
  
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
  spatial_args_new = list("anisotropic_mesh"=x$spatial_list$MeshList$anisotropic_mesh, "Kmeans"=x$spatial_list$Kmeans, "Lon_i"=Lon_i, "Lat_i"=Lat_i )
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
  model_args_default = list("TmbData"=data_list, "RunDir"=working_dir, "Version"=x$settings$Version, "RhoConfig"=x$settings$RhoConfig, "loc_x"=spatial_list$loc_x, "Method"=spatial_list$Method, "Map" = x$tmb_list$Map)
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

match_strata_fn_aja <- function(points, strata_dataframe, index_shapes) {
  if(FALSE){
    points = Tmp
    strata_dataframe = strata.limits[l, , drop = FALSE]
    index_shapes = index_shapes
  }
  if(is.null(index_shapes)){
    # Default all strata
    match_latitude_TF = match_longitude_TF = match_depth_TF = rep( TRUE, nrow(strata_dataframe))
    if( all(c("south_border","north_border") %in% names(strata_dataframe)) ){
      match_latitude_TF = as.numeric(x["BEST_LAT_DD"])>strata_dataframe[,'south_border'] & as.numeric(x["BEST_LAT_DD"])<=strata_dataframe[,'north_border']
    }
    if( all(c("west_border","east_border") %in% names(strata_dataframe)) ){
      match_longitude_TF = as.numeric(x["BEST_LON_DD"])>strata_dataframe[,'west_border'] & as.numeric(x["BEST_LON_DD"])<=strata_dataframe[,'east_border']
    }
    if( all(c("shallow_border","deep_border") %in% names(strata_dataframe)) ){
      match_depth_TF = as.numeric(x["BEST_DEPTH_M"])>strata_dataframe[,'shallow_border'] & as.numeric(x["BEST_DEPTH_M"])<=strata_dataframe[,'deep_border']
    }
    # Return stuff
    Char = as.character(strata_dataframe[match_latitude_TF & match_longitude_TF & match_depth_TF,"STRATA"]) 
    return(ifelse(length(Char)==0,NA,Char))
  }
  
  # Andrew edit...
  if(!is.null(index_shapes)){
    Tmp_sf<- data.frame(points) %>%
      st_as_sf(., coords = c("BEST_LON_DD", "BEST_LAT_DD"), crs = st_crs(index_shapes), remove = FALSE)
    match_shape<- Tmp_sf %>%
      st_join(., index_shapes, join = st_within) %>%
      mutate(., "Row_ID" = seq(from = 1, to = nrow(.))) %>%
      st_drop_geometry() %>%
      dplyr::select(., Region) %>%
      as.vector()
    return(match_shape)
  }
}

Prepare_User_Extrapolation_Data_Fn_aja<- function (input_grid, strata.limits = NULL, projargs = NA, zone = NA, flip_around_dateline = TRUE, index_shapes, ...) {
  if (is.null(strata.limits)) {
    strata.limits = data.frame(STRATA = "All_areas")
  }
  message("Using strata ", strata.limits)
  Data_Extrap <- input_grid
  Area_km2_x = Data_Extrap[, "Area_km2"]
  Tmp = cbind(BEST_LAT_DD = Data_Extrap[, "Lat"], BEST_LON_DD = Data_Extrap[, "Lon"])
  if ("Depth" %in% colnames(Data_Extrap)) {
    Tmp = cbind(Tmp, BEST_DEPTH_M = Data_Extrap[, "Depth"])
  }
  a_el = as.data.frame(matrix(NA, nrow = nrow(Data_Extrap), ncol = nrow(strata.limits), dimnames = list(NULL, strata.limits[, "STRATA"])))
  for (l in 1:ncol(a_el)) {
    a_el[, l] = match_strata_fn_aja(points = Tmp, strata_dataframe = strata.limits[l, , drop = FALSE], index_shapes = index_shapes[l,])
    a_el[, l] = ifelse(is.na(a_el[, l]), 0, Area_km2_x)
  }
  tmpUTM = project_coordinates(X = Data_Extrap[, "Lon"], Y = Data_Extrap[, "Lat"], projargs = projargs, zone = zone, flip_around_dateline = flip_around_dateline)
  Data_Extrap = cbind(Data_Extrap, Include = 1)
  if (all(c("E_km", "N_km") %in% colnames(Data_Extrap))) {
    Data_Extrap[, c("E_km", "N_km")] = tmpUTM[, c("X", "Y")]
  } else {
    Data_Extrap = cbind(Data_Extrap, E_km = tmpUTM[, "X"], N_km = tmpUTM[, "Y"])
  }
  Return = list(a_el = a_el, Data_Extrap = Data_Extrap, zone = attr(tmpUTM, "zone"), projargs = attr(tmpUTM, "projargs"), flip_around_dateline = flip_around_dateline, Area_km2_x = Area_km2_x)
  return(Return)
}

make_extrapolation_info_aja<- function (Region, projargs = NA, zone = NA, strata.limits = data.frame(STRATA = "All_areas"), create_strata_per_region = FALSE, max_cells = NULL, input_grid = NULL, observations_LL = NULL, grid_dim_km = c(2, 2), maximum_distance_from_sample = NULL, grid_in_UTM = TRUE, grid_dim_LL = c(0.1, 0.1), region = c("south_coast", "west_coast"), strata_to_use = c("SOG", "WCVI", "QCS", "HS", "WCHG"), epu_to_use = c("All", "Georges_Bank", "Mid_Atlantic_Bight", "Scotian_Shelf", "Gulf_of_Maine", "Other")[1], survey = "Chatham_rise", surveyname = "propInWCGBTS", flip_around_dateline, nstart = 100, area_tolerance = 0.05, backwards_compatible_kmeans = FALSE, DirPath = paste0(getwd(), "/"), index_shapes, ...) {
  if (is.null(max_cells)) 
    max_cells = Inf
  for (rI in seq_along(Region)) {
    Extrapolation_List = NULL
    if (tolower(Region[rI]) == "california_current") {
      if (missing(flip_around_dateline)) 
        flip_around_dateline = FALSE
      Extrapolation_List = Prepare_WCGBTS_Extrapolation_Data_Fn(strata.limits = strata.limits, surveyname = surveyname, projargs = projargs, zone = zone, flip_around_dateline = flip_around_dateline, ...)
    }
    if (tolower(Region[rI]) %in% c("wcghl", "wcghl_domain", "west_coast_hook_and_line")) {
      if (missing(flip_around_dateline)) 
        flip_around_dateline = FALSE
      Extrapolation_List = Prepare_WCGHL_Extrapolation_Data_Fn(strata.limits = strata.limits, projargs = projargs, zone = zone, flip_around_dateline = flip_around_dateline, ...)
    }
    if (tolower(Region[rI]) == "british_columbia") {
      if (missing(flip_around_dateline)) 
        flip_around_dateline = FALSE
      Extrapolation_List = Prepare_BC_Coast_Extrapolation_Data_Fn(strata.limits = strata.limits, strata_to_use = strata_to_use, projargs = projargs, zone = zone, flip_around_dateline = flip_around_dateline, ...)
    }
    if (tolower(Region[rI]) == "eastern_bering_sea") {
      if (missing(flip_around_dateline)) 
        flip_around_dateline = TRUE
      Extrapolation_List = Prepare_EBS_Extrapolation_Data_Fn(strata.limits = strata.limits, projargs = projargs, zone = zone, flip_around_dateline = flip_around_dateline, ...)
    }
    if (tolower(Region[rI]) == "northern_bering_sea") {
      if (missing(flip_around_dateline)) 
        flip_around_dateline = FALSE
      Extrapolation_List = Prepare_NBS_Extrapolation_Data_Fn(strata.limits = strata.limits, projargs = projargs, zone = zone, flip_around_dateline = flip_around_dateline, ...)
    }
    if (tolower(Region[rI]) == "bering_sea_slope") {
      if (missing(flip_around_dateline)) 
        flip_around_dateline = FALSE
      Extrapolation_List = Prepare_BSslope_Extrapolation_Data_Fn(strata.limits = strata.limits, projargs = projargs, zone = zone, flip_around_dateline = flip_around_dateline, ...)
    }
    if (tolower(Region[rI]) %in% c("st_matthews_island", "smi")) {
      if (missing(flip_around_dateline)) 
        flip_around_dateline = TRUE
      Extrapolation_List = Prepare_SMI_Extrapolation_Data_Fn(strata.limits = strata.limits, projargs = projargs, zone = zone, flip_around_dateline = flip_around_dateline, ...)
    }
    if (tolower(Region[rI]) == "aleutian_islands") {
      if (missing(flip_around_dateline)) 
        flip_around_dateline = TRUE
      Extrapolation_List = Prepare_AI_Extrapolation_Data_Fn(strata.limits = strata.limits, projargs = projargs, zone = zone, flip_around_dateline = flip_around_dateline, ...)
    }
    if (tolower(Region[rI]) == "gulf_of_alaska") {
      if (missing(flip_around_dateline)) 
        flip_around_dateline = FALSE
      Extrapolation_List = Prepare_GOA_Extrapolation_Data_Fn(strata.limits = strata.limits, projargs = projargs, zone = zone, flip_around_dateline = flip_around_dateline, ...)
    }
    if (tolower(Region[rI]) == tolower("BFISH_MHI")) {
      if (missing(flip_around_dateline)) 
        flip_around_dateline = FALSE
      Extrapolation_List = Prepare_BFISH_MHI_Extrapolation_Data_Fn(strata.limits = strata.limits, projargs = projargs, zone = zone, flip_around_dateline = flip_around_dateline, ...)
    }
    if (tolower(Region[rI]) == "northwest_atlantic") {
      if (missing(flip_around_dateline)) 
        flip_around_dateline = FALSE
      Extrapolation_List = Prepare_NWA_Extrapolation_Data_Fn(strata.limits = strata.limits, epu_to_use = epu_to_use, projargs = projargs, zone = zone, flip_around_dateline = flip_around_dateline, ...)
    }
    if (tolower(Region[rI]) == "south_africa") {
      if (missing(flip_around_dateline)) 
        flip_around_dateline = FALSE
      Extrapolation_List = Prepare_SA_Extrapolation_Data_Fn(strata.limits = strata.limits, region = region, projargs = projargs, zone = zone, flip_around_dateline = flip_around_dateline, ...)
    }
    if (tolower(Region[rI]) == "gulf_of_st_lawrence") {
      if (missing(flip_around_dateline)) 
        flip_around_dateline = FALSE
      Extrapolation_List = Prepare_GSL_Extrapolation_Data_Fn(strata.limits = strata.limits, projargs = projargs, zone = zone, flip_around_dateline = flip_around_dateline, ...)
    }
    if (tolower(Region[rI]) == "new_zealand") {
      if (missing(flip_around_dateline)) 
        flip_around_dateline = FALSE
      Extrapolation_List = Prepare_NZ_Extrapolation_Data_Fn(strata.limits = strata.limits, survey = survey, projargs = projargs, zone = zone, flip_around_dateline = flip_around_dateline, ...)
    }
    if (tolower(Region[rI]) == "habcam") {
      if (missing(flip_around_dateline)) 
        flip_around_dateline = FALSE
      Extrapolation_List = Prepare_HabCam_Extrapolation_Data_Fn(strata.limits = strata.limits, projargs = projargs, zone = zone, flip_around_dateline = flip_around_dateline, ...)
    }
    if (tolower(Region[rI]) == "gulf_of_mexico") {
      if (missing(flip_around_dateline)) 
        flip_around_dateline = FALSE
      Extrapolation_List = Prepare_GOM_Extrapolation_Data_Fn(strata.limits = strata.limits, projargs = projargs, zone = zone, flip_around_dateline = flip_around_dateline, ...)
    }
    Shapefile_set = c("ATL-IBTS-Q1", "ATL-IBTS-Q4", "BITS", 
                      "BTS", "BTS-VIIA", "EVHOE", "IE-IGFS", "NIGFS", "NS_IBTS", 
                      "PT-IBTS", "SP-ARSA", "SP-NORTH", "SP-PORC", "CalCOFI_Winter-Spring", 
                      "CalCOFI-IMECOCAL_Winter-Spring", "IMECOCAL_Winter-Spring", 
                      "CalCOFI-IMECOCAL_Summer", "rockfish_recruitment_coastwide", 
                      "rockfish_recruitment_core")
    if (toupper(Region[rI]) %in% toupper(Shapefile_set)) {
      if (Region[rI] == "SP-ARSA") 
        stop("There's some problem with `SP-ARSA` which precludes it's use")
      Conversion = convert_shapefile(file_path = paste0(system.file("region_shapefiles", package = "FishStatsUtils"), "/", toupper(Region[rI]), "/Shapefile.shp"), projargs_for_shapefile = "+proj=longlat +ellps=WGS84 +no_defs", projargs = projargs, grid_dim_km = grid_dim_km, area_tolerance = area_tolerance, ...)
      Extrapolation_List = list(a_el = matrix(Conversion$extrapolation_grid[, "Area_km2"], ncol = 1), Data_Extrap = Conversion$extrapolation_grid, zone = NA, projargs = Conversion$projargs, flip_around_dateline = FALSE, Area_km2_x = Conversion$extrapolation_grid[, "Area_km2"])
    }
    if (file.exists(Region[rI])) {
      Conversion = convert_shapefile(file_path = Region[rI], projargs = projargs, grid_dim_km = grid_dim_km, area_tolerance = area_tolerance, ...)
      Extrapolation_List = list(a_el = matrix(Conversion$extrapolation_grid[, "Area_km2"], ncol = 1), Data_Extrap = Conversion$extrapolation_grid, zone = NA, projargs = Conversion$projargs, flip_around_dateline = FALSE, Area_km2_x = Conversion$extrapolation_grid[, "Area_km2"])
    }
    if (tolower(Region[rI]) == "stream_network") {
      if (is.null(input_grid)) {
        stop("Because you're using a stream network, please provide 'input_grid' input")
      }
      if (!(all(c("Lat", "Lon", "Area_km2", "child_i") %in% colnames(input_grid)))) {
        stop("'input_grid' must contain columns named 'Lat', 'Lon', 'Area_km2', and 'child_i'")
      }
      if (missing(flip_around_dateline)) 
        flip_around_dateline = FALSE
      Extrapolation_List = Prepare_User_Extrapolation_Data_Fn(strata.limits = strata.limits, input_grid = input_grid, projargs = projargs, zone = zone, flip_around_dateline = flip_around_dateline, ...)
    }
    if (tolower(Region[rI]) == "user") {
      if (is.null(input_grid)) {
        stop("Because you're using a user-supplied region, please provide 'input_grid' input")
      }
      if (!(all(c("Lat", "Lon", "Area_km2") %in% colnames(input_grid)))) {
        stop("'input_grid' must contain columns named 'Lat', 'Lon', and 'Area_km2'")
      }
      if (missing(flip_around_dateline)) 
        flip_around_dateline = FALSE
      Extrapolation_List = Prepare_User_Extrapolation_Data_Fn_aja(strata.limits = strata.limits, input_grid = input_grid, projargs = projargs, zone = zone, flip_around_dateline = flip_around_dateline, index_shapes = index_shapes, ...)
    }
    if (is.null(Extrapolation_List)) {
      if (is.null(observations_LL)) {
        stop("Because you're using a new Region[rI], please provide 'observations_LL' input with columns named `Lat` and `Lon`")
      }
      if (missing(flip_around_dateline)) 
        flip_around_dateline = FALSE
      Extrapolation_List = Prepare_Other_Extrapolation_Data_Fn(strata.limits = strata.limits, observations_LL = observations_LL, grid_dim_km = grid_dim_km, maximum_distance_from_sample = maximum_distance_from_sample, grid_in_UTM = grid_in_UTM, grid_dim_LL = grid_dim_LL, projargs = projargs, zone = zone, flip_around_dateline = flip_around_dateline, ...)
    }
    if (rI == 1) {
      Return = Extrapolation_List
    } else {
      Return = combine_extrapolation_info(Return, Extrapolation_List, create_strata_per_region = create_strata_per_region)
    }
  }
  if (max_cells < nrow(Return$Data_Extrap)) {
    message("# Reducing extrapolation-grid from ", nrow(Return$Data_Extrap), " to ", max_cells, " cells for Region(s): ", paste(Region, collapse = ", "))
    loc_orig = Return$Data_Extrap[, c("E_km", "N_km")]
    loc_orig = loc_orig[which(Return$Area_km2_x > 0), ]
    Kmeans = make_kmeans(n_x = max_cells, loc_orig = loc_orig, nstart = nstart, randomseed = 1, iter.max = 1000, DirPath = DirPath, Save_Results = TRUE, kmeans_purpose = "extrapolation", backwards_compatible_kmeans = backwards_compatible_kmeans)
    Kmeans[["cluster"]] = RANN::nn2(data = Kmeans[["centers"]], query = Return$Data_Extrap[, c("E_km", "N_km")], k = 1)$nn.idx[, 1]
    aggregate_vector = function(values_x, index_x, max_index, FUN = sum) { 
      tapply(values_x, INDEX = factor(index_x, levels = 1:max_index), FUN = FUN)
    }
    a_el = matrix(NA, nrow = max_cells, ncol = ncol(Return$a_el))
    for (lI in 1:ncol(Return$a_el)) {
      a_el[, lI] = aggregate_vector(values_x = Return$a_el[, lI], index_x = Kmeans$cluster, max_index = max_cells)
    }
    Area_km2_x = aggregate_vector(values_x = Return$Area_km2_x, index_x = Kmeans$cluster, max_index = max_cells)
    Include = aggregate_vector(values_x = Return$Data_Extrap[, "Include"], index_x = Kmeans$cluster, max_index = max_cells, FUN = function(vec) {
      any(vec > 0)
      })
    lonlat_g = project_coordinates(X = Kmeans$centers[, "E_km"], Y = Kmeans$centers[, "N_km"], projargs = "+proj=longlat +ellps=WGS84", origargs = Return$projargs)
    Data_Extrap = cbind(Lon = lonlat_g[, 1], Lat = lonlat_g[, 2], Include = Include, Kmeans$centers)
    Return = list(a_el = a_el, Data_Extrap = Data_Extrap, zone = Return$zone, projargs = Return$projargs, flip_around_dateline = Return$flip_around_dateline, Area_km2_x = Area_km2_x)
  }
  if (length(Region) > 1 & create_strata_per_region == TRUE) {
    Return$a_el = cbind(Total = rowSums(Return$a_el), Return$a_el)
  }
  class(Return) = "make_extrapolation_info"
  return(Return)
}

fit_model_aja<- function (settings, Lat_i, Lon_i, t_i, b_i, a_i, c_iz = rep(0, length(b_i)), v_i = rep(0, length(b_i)), working_dir = paste0(getwd(), "/"), X1config_cp = NULL, X2config_cp = NULL, covariate_data, X1_formula = ~0, X2_formula = ~0, Q1config_k = NULL, Q2config_k = NULL, catchability_data, Q1_formula = ~0, Q2_formula = ~0, newtonsteps = 1, silent = TRUE, build_model = TRUE, run_model = TRUE, test_fit = TRUE, ...) {
  if(FALSE){
    vast_build_adjust = vast_adjust
    nmfs_species_code = nmfs_species_code
    out_dir = here::here("results/mod_fits")
    index_shapes = index_shapefiles
    
    tar_load(vast_settings)
    settings = vast_settings
    tar_load(vast_extrap_grid)
    extrap_grid = input_grid
    tar_load(vast_sample_data)
    sample_data = vast_sample_data
    
    
    "settings" = vast_settings
    "input_grid" = extrap_grid
    "Lat_i" = sample_data[, 'Lat']
    "Lon_i" = sample_data[, 'Lon']
    "t_i" = sample_data[, 'Year']
    "c_i" = rep(0, nrow(sample_data))
    "b_i" = sample_data[, 'Biomass']
    "c_iz" = rep(0, length(b_i))
    "v_i" = rep(0, length(b_i))
    "a_i" = sample_data[, 'Swept']
    "PredTF_i" = sample_data[, 'Pred_TF']
    "X1config_cp" = Xconfig_list[['X1config_cp']]
    "X2config_cp" = Xconfig_list[['X2config_cp']]
    "covariate_data" = covariate_data
    "X1_formula" = X1_formula
    "X2_formula" = X2_formula
    "X_contrasts" = X_contrasts
    "catchability_data" = catchability_data
    "Q1_formula" = Q1_formula
    "Q2_formula" = Q2_formula
    "Q1config_k" = Xconfig_list[['Q1config_k']]
    "Q2config_k" = Xconfig_list[['Q2config_k']]
    "newtonsteps" = 1
    "getsd" = TRUE
    "getReportCovariance" = TRUE
    "run_model" = FALSE
    "test_fit" = FALSE
    "Use_REML" = FALSE
    "getJointPrecision" = FALSE
    "index_shapes" = index_shapefiles
  }
  extra_args = list(...)
  extra_args = c(extra_args, extra_args$extrapolation_args, extra_args$spatial_args, extra_args$optimize_args, extra_args$model_args)
  data_frame = data.frame(Lat_i = Lat_i, Lon_i = Lon_i, a_i = a_i, v_i = v_i, b_i = b_i, t_i = t_i, c_iz = c_iz)
  year_labels = seq(min(t_i), max(t_i))
  years_to_plot = which(year_labels %in% t_i)
  message("\n### Writing output from `fit_model` in directory: ", working_dir)
  dir.create(working_dir, showWarnings = FALSE, recursive = TRUE)
  capture.output(settings, file = file.path(working_dir, "settings.txt"))
  message("\n### Making extrapolation-grid")
  extrapolation_args_default = list(Region = settings$Region, strata.limits = settings$strata.limits, zone = settings$zone, max_cells = settings$max_cells, DirPath = working_dir)
  extrapolation_args_input = combine_lists(input = extra_args, default = extrapolation_args_default, args_to_use = formalArgs(make_extrapolation_info_aja))
  extrapolation_list = do.call(what = make_extrapolation_info_aja, args = extrapolation_args_input)
  message("\n### Making spatial information")
  spatial_args_default = list(grid_size_km = settings$grid_size_km, n_x = settings$n_x, Method = settings$Method, Lon_i = Lon_i, Lat_i = Lat_i, Extrapolation_List = extrapolation_list, DirPath = working_dir, Save_Results = TRUE, fine_scale = settings$fine_scale, knot_method = settings$knot_method)
  spatial_args_input = combine_lists(input = extra_args, default = spatial_args_default, args_to_use = c(formalArgs(make_spatial_info), formalArgs(INLA::inla.mesh.create)))
  spatial_list = do.call(what = make_spatial_info, args = spatial_args_input)
  message("\n### Making data object")
  if (missing(covariate_data)) 
    covariate_data = NULL
  if (missing(catchability_data)) 
    catchability_data = NULL
  data_args_default = list(Version = settings$Version, FieldConfig = settings$FieldConfig, OverdispersionConfig = settings$OverdispersionConfig, RhoConfig = settings$RhoConfig, VamConfig = settings$VamConfig, ObsModel = settings$ObsModel, c_iz = c_iz, b_i = b_i, a_i = a_i, v_i = v_i, s_i = spatial_list$knot_i - 1, t_i = t_i, spatial_list = spatial_list, Options = settings$Options, Aniso = settings$use_anisotropy, X1config_cp = X1config_cp, X2config_cp = X2config_cp, covariate_data = covariate_data, X1_formula = X1_formula, X2_formula = X2_formula, Q1config_k = Q1config_k, Q2config_k = Q2config_k, catchability_data = catchability_data, Q1_formula = Q1_formula, Q2_formula = Q2_formula)
  data_args_input = combine_lists(input = extra_args, default = data_args_default)
  data_list = do.call(what = make_data, args = data_args_input)
  message("\n### Making TMB object")
  model_args_default = list(TmbData = data_list, RunDir = working_dir, Version = settings$Version, RhoConfig = settings$RhoConfig, loc_x = spatial_list$loc_x, Method = spatial_list$Method, build_model = build_model)
  model_args_input = combine_lists(input = extra_args, default = model_args_default, args_to_use = formalArgs(make_model))
  tmb_list = do.call(what = make_model, args = model_args_input)
  if (run_model == FALSE | build_model == FALSE) {
    input_args = list(extra_args = extra_args, extrapolation_args_input = extrapolation_args_input, model_args_input = model_args_input, spatial_args_input = spatial_args_input, data_args_input = data_args_input)
    Return = list(data_frame = data_frame, extrapolation_list = extrapolation_list, spatial_list = spatial_list, data_list = data_list, tmb_list = tmb_list, year_labels = year_labels, years_to_plot = years_to_plot, settings = settings, input_args = input_args)
    class(Return) = "fit_model"
    return(Return)
  }
  if (silent == TRUE) 
    tmb_list$Obj$env$beSilent()
  if (test_fit == TRUE) {
    message("\n### Testing model at initial values")
    LogLike0 = tmb_list$Obj$fn(tmb_list$Obj$par)
    Gradient0 = tmb_list$Obj$gr(tmb_list$Obj$par)
    if (any(Gradient0 == 0)) {
      message("\n")
      stop("Please check model structure; some parameter has a gradient of zero at starting values\n", 
           call. = FALSE)
    } else {
      message("Looks good: All fixed effects have a nonzero gradient")
    }
  }
  message("\n### Estimating parameters")
  optimize_args_default1 = list(lower = tmb_list$Lower, upper = tmb_list$Upper, loopnum = 2)
  optimize_args_default1 = combine_lists(default = optimize_args_default1, input = extra_args, args_to_use = formalArgs(TMBhelper::fit_tmb))
  optimize_args_input1 = list(obj = tmb_list$Obj, savedir = NULL, newtonsteps = 0, bias.correct = FALSE, control = list(eval.max = 10000, iter.max = 10000, trace = 1), quiet = TRUE, getsd = FALSE)
  optimize_args_input1 = combine_lists(default = optimize_args_default1, input = optimize_args_input1, args_to_use = formalArgs(TMBhelper::fit_tmb))
  parameter_estimates = do.call(what = TMBhelper::fit_tmb, args = optimize_args_input1)
  if (exists("check_fit") & test_fit == TRUE) {
    problem_found = VAST::check_fit(parameter_estimates)
    if (problem_found == TRUE) {
      message("\n")
      stop("Please change model structure to avoid problems with parameter estimates and then re-try; see details in `?check_fit`\n", call. = FALSE)
    }
  }
  optimize_args_default2 = list(obj = tmb_list$Obj, lower = tmb_list$Lower, upper = tmb_list$Upper, savedir = working_dir, bias.correct = settings$bias.correct, newtonsteps = newtonsteps, bias.correct.control = list(sd = FALSE, split = NULL, nsplit = 1, vars_to_correct = settings$vars_to_correct), control = list(eval.max = 10000, iter.max = 10000, trace = 1), loopnum = 1, getJointPrecision = TRUE)
  optimize_args_input2 = combine_lists(input = extra_args, default = optimize_args_default2, args_to_use = formalArgs(TMBhelper::fit_tmb))
  optimize_args_input2 = combine_lists(input = list(startpar = parameter_estimates$par), default = optimize_args_input2)
  parameter_estimates = do.call(what = TMBhelper::fit_tmb, args = optimize_args_input2)
  if ("par" %in% names(parameter_estimates)) {
    Report = tmb_list$Obj$report()
    ParHat = tmb_list$Obj$env$parList(parameter_estimates$par)
  } else {
    Report = ParHat = "Model is not converged"
  }
  input_args = list(extra_args = extra_args, extrapolation_args_input = extrapolation_args_input, model_args_input = model_args_input, spatial_args_input = spatial_args_input, optimize_args_input1 = optimize_args_input1, optimize_args_input2 = optimize_args_input2, data_args_input = data_args_input)
  Return = list(data_frame = data_frame, extrapolation_list = extrapolation_list, spatial_list = spatial_list, data_list = data_list, tmb_list = tmb_list, parameter_estimates = parameter_estimates, Report = Report, ParHat = ParHat, year_labels = year_labels, years_to_plot = years_to_plot, settings = settings, input_args = input_args, X1config_cp = X1config_cp, X2config_cp = X2config_cp, covariate_data = covariate_data, X1_formula = X1_formula, X2_formula = X2_formula, Q1config_k = Q1config_k, Q2config_k = Q1config_k, catchability_data = catchability_data, Q1_formula = Q1_formula, Q2_formula = Q2_formula)
  Return$effects = list()
  if (!is.null(catchability_data)) {
    catchability_data_full = data.frame(catchability_data, linear_predictor = 0)
    Q1_formula_full = update.formula(Q1_formula, linear_predictor ~ . + 0)
    call_Q1 = lm(Q1_formula_full, data = catchability_data_full)$call
    Q2_formula_full = update.formula(Q2_formula, linear_predictor ~ . + 0)
    call_Q2 = lm(Q2_formula_full, data = catchability_data_full)$call
    Return$effects = c(Return$effects, list(call_Q1 = call_Q1, call_Q2 = call_Q2, catchability_data_full = catchability_data_full))
  }
  if (!is.null(covariate_data)) {
    covariate_data_full = data.frame(covariate_data, linear_predictor = 0)
    X1_formula_full = update.formula(X1_formula, linear_predictor ~ . + 0)
    call_X1 = lm(X1_formula_full, data = covariate_data_full)$call
    X2_formula_full = update.formula(X2_formula, linear_predictor ~ . + 0)
    call_X2 = lm(X2_formula_full, data = covariate_data_full)$call
    Return$effects = c(Return$effects, list(call_X1 = call_X1, call_X2 = call_X2, covariate_data_full = covariate_data_full))
  }
  class(Return) = "fit_model"
  return(Return)
}

vast_read_region_shape<- function(region_shapefile_dir){
  region_file<- list.files(region_shapefile_dir, pattern = ".shp", full.names = TRUE)
  region_sf<- st_read(region_file)
  return(region_sf)
} 

vast_read_index_shapes<- function(index_shapefiles_dir){
  
  if(FALSE){
    index_shapefiles_dir<- "~/GitHub/sdm_workflow/scratch/aja/TargetsSDM/data/supporting/index_shapefiles/"
  }
  
  index_files<- list.files(index_shapefiles_dir, pattern = ".shp", full.names = TRUE)
  
  for(i in seq_along(index_files)){
    index_shapes_temp<- st_read(index_files[i])
    if(i == 1){
      index_shapes_out<- index_shapes_temp
    } else {
      index_shapes_out<- bind_rows(index_shapes_out, index_shapes_temp)
    }
  }
  return(index_shapes_out)
}

######
## Getting abundance index time series
######
get_vast_index_timeseries<- function(vast_fit, nice_category_names, index_scale = c("raw", "log"), out_dir){
  
  if(FALSE){
    tar_load(vast_fit)
    nice_category_names = "American lobster"
    index_scale = "log"
    out_dir = here::here("scratch/aja/TargetsSDM/results/tables")
  }
  
  TmbData<- vast_fit$data_list
  Sdreport<- vast_fit$parameter_estimates$SD
  
  # Time series steps
  time_ind<- 1:TmbData$n_t
  time_labels<- sort(unique(vast_fit$data_frame$t_i)[time_ind])
  
  # Index regions
  index_regions_ind<- 1:TmbData$n_l
  index_regions<- vast_fit$settings$strata.limits$STRATA[index_regions_ind]
  
  # Categories
  categories_ind<- 1:TmbData$n_c
  
  # Get the index information
  SD<- TMB::summary.sdreport(Sdreport)
  SD_stderr<- TMB:::as.list.sdreport(Sdreport, what = "Std. Error", report = TRUE)
  SD_estimate<- TMB:::as.list.sdreport(Sdreport, what = "Estimate", report = TRUE)
  if(vast_fit$settings$bias.correct == TRUE && "unbiased" %in% names(Sdreport)){
    SD_estimate_biascorrect<- TMB:::as.list.sdreport(Sdreport, what = "Std. (bias.correct)", report = TRUE)
  }
  
  # Now, populate array with values
  Index_ctl = log_Index_ctl = array(NA, dim = c(unlist(TmbData[c('n_c','n_t','n_l')]), 2), dimnames = list(categories_ind, time_labels, index_regions, c('Estimate','Std. Error')))
  
  if(index_scale == "raw"){
    if(vast_fit$settings$bias.correct == TRUE && "unbiased" %in% names(Sdreport)){
      Index_ctl[] = SD[which(rownames(SD) == "Index_ctl"),c('Est. (bias.correct)','Std. Error')]
    } else {
      Index_ctl[]<- SD[which(rownames(SD) == "Index_ctl"), c('Estimate','Std. Error')]
    }
    index_res_array<- Index_ctl
  } else {
    if(vast_fit$settings$bias.correct == TRUE && "unbiased" %in% names(Sdreport)){
      log_Index_ctl[] = SD[which(rownames(SD) == "ln_Index_ctl"),c('Est. (bias.correct)','Std. Error')]
    } else {
      log_Index_ctl[]<- SD[which(rownames(SD) == "ln_Index_ctl"), c('Estimate','Std. Error')]
    }
    index_res_array<- log_Index_ctl
  }
  
  # Data manipulation to get out out the array and to something more "plottable"
  for(i in seq_along(categories_ind)){
    index_array_temp<- index_res_array[i, , , ]
    index_res_temp_est<- data.frame("Time" = as.numeric(rownames(index_array_temp[,,1])), "Category" = categories_ind[i], index_array_temp[,,1]) %>%
      gather(., "Index_Region", "Index_Estimate", -Time, -Category)
    index_res_temp_sd<- data.frame("Time" = as.numeric(rownames(index_array_temp[,,1])), "Category" = categories_ind[i], index_array_temp[,,2]) %>%
      gather(., "Index_Region", "Index_SD", -Time, -Category)
    index_res_temp_out<- index_res_temp_est %>%
      left_join(., index_res_temp_sd)
    
    if(i == 1){
      index_res_out<- index_res_temp_out
    } else {
      index_res_out<- bind_rows(index_res_out, index_res_temp_out)
    }
  }
  
  # Get date info instead of time..
  year_start<- min(as.numeric(as.character(vast_fit$covariate_data$Year_Cov)))
  seasons<- nlevels(unique(vast_fit$covariate_data$Season))
  
  if(seasons == 3 & max(time_labels) == 347){
    time_labels_use<- paste(rep(seq(from = year_start, to = 2100), each = 3), rep(c("SPRING", "SUMMER", "FALL")), sep = "-")
  }
  
  index_res_out$Date<- factor(c(time_labels_use, time_labels_use), levels = time_labels_use)
  
  
  # Save and return it
  write.csv(index_res_out, file = paste(out_dir, "/Biomass_Index_", index_scale, "_", nice_category_names, ".csv", sep = ""))
  return(index_res_out)
}

plot_vast_index_timeseries<- function(index_res_df, index_scale, nice_category_names, nice_xlab, nice_ylab, paneling = c("Category", "Index_Region", "None"), color_pal = c('#66c2a5','#fc8d62','#8da0cb'), out_dir){
  
  if(FALSE){
    tar_load(biomass_indices)
    index_res_df<- index_res_out
    index_res_df<- biomass_indices
    nice_category_names<- "American lobster"
    nice_xlab = "Year"
    nice_ylab = "Biomass index (metric tons)"
    color_pal = NULL
    paneling<- "none"
    date_breaks<- "5 year"
    out_dir = here::here("scratch/aja/TargetsSDM/results/")
  }
  
  if(paneling == "none"){
    if(!is.null(color_pal)){
      colors_use<- color_pal
    } else {
      color_pal<- c('#66c2a5','#fc8d62','#8da0cb','#e78ac3','#a6d854')
      colors_use<- color_pal[1:length(unique(index_res_df$Index_Region))]
    }
    
    index_res_df<- index_res_df %>%
      mutate(., Index_Region = factor(Index_Region, levels = unique(Index_Region), labels = unique(Index_Region)),
             Year = as.numeric(gsub("([0-9]+).*$", "\\1", Date)))
    
    # Date info
    index_res_df$Date<- as.Date(paste(index_res_df$Year, ifelse(grepl("SPRING", index_res_df$Date), "-04-15",
                                                        ifelse(grepl("SUMMER", index_res_df$Date), "-07-15", "-10-15")), sep = ""))
    
    
    plot_out<- ggplot() +
      geom_errorbar(data = index_res_df, aes(x = Date, ymin = (Index_Estimate - Index_SD), ymax = (Index_Estimate + Index_SD), color = Index_Region, group = Index_Region)) + 
      geom_point(data = index_res_df, aes(x = Date, y = Index_Estimate, color = Index_Region)) +
      scale_color_manual(values = colors_use) +
      scale_x_date(date_breaks = "5 year", date_labels =  "%Y", limits = c(min(index_res_df$Date), max(index_res_df$Date)), expand = c(0, 0)) +
      xlab({{nice_xlab}}) +
      ylab({{nice_ylab}}) +
      ggtitle({{nice_category_names}}) + 
      theme_bw() +
      theme(legend.title = element_blank(), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) 
  }
  
  # Save and return the plot
  ggsave(plot_out, file = paste(out_dir, "/Biomass_Index_", index_scale, "_", nice_category_names, ".jpg", sep = ""))
}


######
## Plot parameter effects...
######
plot_vast_parameter_curves<- function(vast_fit, nice_category_names){
  if(FALSE){
    tar_load(vast_fit)
    nice_category_names = "American lobster"
  }
  
  # Must add data-frames to global environment (hope to fix in future)
  if(is.null(vast_fit$effects)){
    vast_fit$effects<- list()
    if(!is.null(vast_fit$catchability_data)) {
      catchability_data_full = data.frame(vast_fit$catchability_data, linear_predictor = 0)
      Q1_formula_full = update.formula(vast_fit$Q1_formula, linear_predictor ~ . + 0)
      call_Q1 = lm(Q1_formula_full, data = catchability_data_full)$call
      Q2_formula_full = update.formula(vast_fit$Q2_formula, linear_predictor ~ . + 0)
      call_Q2 = lm(Q2_formula_full, data = catchability_data_full)$call
      vast_fit$effects = c(Return$effects, list(call_Q1 = call_Q1, call_Q2 = call_Q2, catchability_data_full = catchability_data_full))
    }
    if (!is.null(vast_fit$covariate_data)) {
      covariate_data_full = data.frame(vast_fit$covariate_data, linear_predictor = 0)
      X1_formula_full = update.formula(vast_fit$X1_formula, linear_predictor ~  . + 0)
      call_X1 = lm(X1_formula_full, data = covariate_data_full)$call
      X2_formula_full = update.formula(vast_fit$X2_formula, linear_predictor ~  . + 0)
      call_X2 = lm(X2_formula_full, data = covariate_data_full)$call
      vast_fit$effects = c(vast_fit$effects, list(call_X1 = call_X1, call_X2 = call_X2, covariate_data_full = covariate_data_full))
    }
  }
  
  # Must add data-frames to global environment
  covariate_data_full = vast_fit$effects$covariate_data_full
  catchability_data_full = vast_fit$effects$catchability_data_full
  
  # Define formula.
  X1_formula = vast_fit$X1_formula
  X2_formula = vast_fit$X2_formula
  
  # Get effects...
  depth_plot<- data.frame(Effect.fit_model(focal.predictors = c("Depth"), mod = vast_fit, which_formula = "X2", xlevels = 100))
  
  ## What's happening here...breaking into Effect.fit_model
  focal.predictors = c("Depth")
  mod = vast_fit
  which_formula = "X2"
  xlevels = 100
  
  if(mod$data_list$n_c > 1) {
    stop("`Effect.fit_model` is not currently designed for multivariate models")
  }
  if(!all(c("covariate_data_full", "catchability_data_full") %in% ls(.GlobalEnv))) {
    stop("Please load `covariate_data_full` and `catchability_data_full` into global memory")
  }
  if(!requireNamespace("effects")) {
    stop("please install the effects package")
  }
  if(!("effects" %in% names(mod))) {
    stop("`effects` slot not detected in input to `Effects.fit_model`. Please update model using later package version.")
  }
  if(which_formula == "X1") {
    formula_orig = mod$X1_formula
    parname = "gamma1_cp"
    mod$call = mod$effects$call_X1
  } else if(which_formula == "X2") {
    formula_orig = mod$X2_formula
    parname = "gamma2_cp"
    mod$call = mod$effects$call_X2
  } else if(which_formula == "Q1") {
    formula_orig = mod$Q1_formula
    parname = "lambda1_k"
    mod$call = mod$effects$call_Q1
  } else if(which_formula == "Q2") {
    formula_orig = mod$Q2_formula
    parname = "lambda2_k"
    mod$call = mod$effects$call_Q2
  } else {
    stop("Check `which_formula` input")
  }
  whichnum = which(names(mod$parameter_estimates$par) == parname)
  mod$parhat = mod$parameter_estimates$par[whichnum]
  mod$covhat = mod$parameter_estimates$SD$cov.fixed[whichnum, whichnum, drop = FALSE]
  if(parname %in% names(mod$tmb_list$Obj$env$map)) {
    mod$parhat = mod$parhat[mod$tmb_list$Obj$env$map[[parname]]]
    mod$covhat = mod$covhat[mod$tmb_list$Obj$env$map[[parname]], mod$tmb_list$Obj$env$map[[parname]], drop = FALSE]
    mod$parhat = ifelse(is.na(mod$parhat), 0, mod$parhat)
    mod$covhat = ifelse(is.na(mod$covhat), 0, mod$covhat)
  }
  names(mod$parhat)[] = parname
  rownames(mod$covhat) = colnames(mod$covhat) = names(mod$parhat)
  formula_full = stats::update.formula(formula_orig, linear_predictor ~ . + 0)
  mod$coefficients = mod$parhat
  mod$vcov = mod$covhat
  mod$formula = formula_full
  mod$family = stats::gaussian(link = "identity")
  family.fit_model = function(x, ...) x$family
  vcov.fit_model = function(x, ...) x$vcov
  dummyfuns = list(variance = function(mu) mu, initialize = expression(mustart = y +  0.1), dev.resids = function(...) stats::poisson()$dev.res(...))
  fam = mod$family
  for (i in names(dummyfuns)) {
    if (is.null(fam[[i]])) 
      fam[[i]] = dummyfuns[[i]]
  }
  if (length(formals(fam$variance)) > 1) {
    warning("overriding variance function for effects: computed variances may be incorrect")
    fam$variance = dummyfuns$variance
  }
  args = list(call = mod$call, coefficients = mod$coefficients, vcov = mod$vcov, family = fam, formula = formula_full)
  
  ## Effects::Effect.default
  focal.predictors<- focal.predictors
  mod<- mod
  sources<- args
  
  sources<- if(missing(sources)) 
    effSources(mod) else sources
  
  formula<- if(is.null(sources$formula)) 
    insight::find_formula(mod)$conditional else sources$formula
  
  if(is.null(focal.predictors)) 
    return(formula)
  
  cl<- if(is.null(sources$call)) {
    if(isS4(mod)) 
      mod@call else mod$call
  } else sources$call
  cl$formula <- formula
  
  type<- if(is.null(sources$type)) 
    "glm" else sources$type
  fam<- try(family(mod), silent = TRUE)
  if(inherits(fam, "try-error")) 
    fam <- NULL
  if(!is.null(sources$family)) {
    fam <- sources$family
  }
  if(!is.null(fam)) {
    fam$aic <- function(...) NULL
    if (!is.null(fam$variance)) {
      if (length(formals(fam$variance)) > 1) 
        stop("Effect plots are not implemented for families with more than\n             one parameter in the variance function (e.g., negitave binomials).")
    }
  }
  cl$family<- fam
  coefficients<- if(is.null(sources$coefficients)) 
    effCoef(mod) else sources$coefficients
  vcov<- if(is.null(sources$vcov)) 
    as.matrix(vcov(mod, complete = TRUE)) else sources$vcov
  zeta<- if(is.null(sources$zeta)) 
    NULL else sources$zeta
  cl$control<- switch(type, glm = glm.control(epsilon = Inf, maxit = 1), polr = list(maxit = 1), multinom = c(maxit = 1))
  cl$method<- sources$method
  .m<- switch(type, glm = match(c("formula", "data", "family", "contrasts", "subset", "control", "offset"), names(cl), 0L), polr = match(c("formula", "data", "family", "contrasts", "subset", "control", "method"), names(cl), 0L), multinom = match(c("formula", "data", "family", "contrasts", "subset", "family", "maxit", "offset"), names(cl), 0L))
  cl<- cl[c(1L, .m)]
  cl[[1L]]<- as.name(type)
  mod2<- eval(cl)
  mod2$coefficients<- coefficients
  mod2$vcov<- vcov
  if (!is.null(zeta)) 
    mod2$zeta <- zeta
  if (type == "glm") {
    mod2$weights <- as.vector(with(mod2, prior.weights * (family$mu.eta(linear.predictors)^2/family$variance(fitted.values))))
  }
  class(mod2) <- c("fakeeffmod", class(mod2))
  Effect(focal.predictors, mod2, ...)
  
  ## What's happening there...Effect.lm
  # Levels for the different covariates -- this is interesting, there's 100 for Season and Year_Cov. 
  if (is.numeric(xlevels)){
    if (length(xlevels) > 1 || round(xlevels != xlevels)) stop("xlevels must be a single whole number or a list")
    form <- Effect.default(NULL, mod) #returns the fixed-effects formula
    terms <- attr(terms(form), "term.labels")
    predictors <- all.vars(parse(text=terms))
    xlevs <- list()
    for (pred in predictors){
      xlevs[[pred]] <- xlevels
    }
    xlevels <- xlevs
  }
  
  if (!missing(partial.residuals)) residuals <- partial.residuals
  partial.residuals <- residuals
  if (missing(transformation)) 
    transformation <- list(link = family(mod)$linkfun, inverse = family(mod)$linkinv)
  if (missing(fixed.predictors)) fixed.predictors <- NULL
  fixed.predictors <- applyDefaults(fixed.predictors, list(given.values=NULL, typical=mean, apply.typical.to.factors=FALSE, offset=mean), arg="fixed.predictors")
  if (missing(given.values)) given.values <- fixed.predictors$given.values
  # new 1/22/18 to allow for automatical equal weighting of factor levels
  if(!is.null(given.values)){
    if (given.values[1] == "default") given.values <- NULL
    if (given.values[1] == "equal") given.values <- .set.given.equal(mod)}
  # end new code
  if (missing(typical)) typical <- fixed.predictors$typical
  if (missing(offset)) offset <- fixed.predictors$offset
  apply.typical.to.factors <- fixed.predictors$apply.typical.to.factors
  if (!missing(confint)) se <- confint
  confint <- applyDefaults(se, list(compute=TRUE, level=.95, type="pointwise"),
                           onFALSE=list(compute=FALSE, level=.95, type="pointwise"),
                           arg="se")
  se <- confint$compute
  if (missing(confidence.level)) confidence.level <- confint$level
  confidence.type <- match.arg(confint$type, c("pointwise", "Scheffe", "scheffe"))
  default.levels <- NULL # just for backwards compatibility
  data <- if (partial.residuals){
    all.vars <- all.vars(formula(mod))
    expand.model.frame(mod, all.vars)[, all.vars]
  }
  else NULL
  if (!is.null(given.values) && !all(which <- names(given.values) %in% names(coef(mod))))
    stop("given.values (", names(given.values[!which]), ") not in the model")
  off <- if (is.numeric(offset) && length(offset) == 1) offset
  else if (is.function(offset)) {
    mod.off <- model.offset(model.frame(mod))
    if (is.null(mod.off)) 0 else offset(mod.off)
  }
  else stop("offset must be a function or a number")
  formula.rhs <- formula(mod)[[3]]
  if (!missing(x.var)){
    if (!is.numeric(x.var)) {
      x.var.name <- x.var
      x.var <- which(x.var == focal.predictors)
    }
    if (length(x.var) == 0) stop("'", x.var.name, "' is not among the focal predictors")
    if (length(x.var) > 1) stop("x.var argument must be of length 1")
  }
  model.components <- Analyze.model(focal.predictors, mod, xlevels, default.levels, formula.rhs,
                                    partial.residuals=partial.residuals, quantiles=quantiles, x.var=x.var, data=data, typical=typical)
  excluded.predictors <- model.components$excluded.predictors
  predict.data <- model.components$predict.data
  predict.data.all.rounded <- predict.data.all <- if (partial.residuals) na.omit(data[, all.vars(formula(mod))]) else NULL
  factor.levels <- model.components$factor.levels
  factor.cols <- model.components$factor.cols
  n.focal <- model.components$n.focal
  x <- model.components$x
  X.mod <- model.components$X.mod
  cnames <- model.components$cnames
  X <- model.components$X
  x.var <- model.components$x.var
  formula.rhs <- formula(mod)[c(1, 3)]
  Terms <- delete.response(terms(mod))
  mf <- model.frame(Terms, predict.data, xlev = factor.levels, na.action=NULL)
  mod.matrix <- model.matrix(formula.rhs, data = mf, contrasts.arg = mod$contrasts)
  if (is.null(x.var)) partial.residuals <- FALSE
  factors <- sapply(predict.data, is.factor)
  if (partial.residuals){
    for (predictor in focal.predictors[-x.var]){
      if (!factors[predictor]){
        values <- unique(predict.data[, predictor])
        predict.data.all.rounded[, predictor] <- values[apply(outer(predict.data.all[, predictor], values, function(x, y) (x - y)^2), 1, which.min)]
      }
    }
  }
  mod.matrix.all <- model.matrix(mod)
  wts <- weights(mod)
  if (is.null(wts))
    wts <- rep(1, length(residuals(mod)))
  mod.matrix <- Fixup.model.matrix(mod, mod.matrix, mod.matrix.all,
                                   X.mod, factor.cols, cnames, focal.predictors, 
                                   excluded.predictors, typical, given.values, 
                                   apply.typical.to.factors) 
  # 11/3/2017.  Check to see if the model is full rank
  # Compute a basis for the null space, using estimibility package
  null.basis <- estimability::nonest.basis(mod)  # returns basis for null space
  # check to see if each row of mod.matrix is estimable
  is.estimable <- estimability::is.estble(mod.matrix, null.basis) # TRUE if effect is estimable else FALSE
  # substitute 0 for NA in coef vector and compute effects
  scoef <- ifelse(is.na(mod$coefficients), 0L, mod$coefficients)
  effect <- off + mod.matrix %*% scoef
  effect[!is.estimable] <- NA  # set all non-estimable effects to NA
  # end estimability check
  if (partial.residuals){
    res <- na.omit(residuals(mod, type="working"))
    fitted <- na.omit(if (inherits(mod, "glm")) predict(mod, type="link") else predict(mod))
    partial.residuals.range <- range(fitted + res)
  }
  else {
    res <- partial.residuals.range <- NULL
  }
  result <- list(term = paste(focal.predictors, collapse="*"),
                 formula = formula(mod), response = response.name(mod),
                 variables = x, fit = effect, x = predict.data[, 1:n.focal, drop=FALSE],
                 x.all=predict.data.all.rounded[, focal.predictors, drop=FALSE],
                 model.matrix = mod.matrix,
                 data = X,
                 discrepancy = 0, offset=off,
                 residuals=res, partial.residuals.range=partial.residuals.range,
                 x.var=x.var)
  if (se) {
    if (any(family(mod)$family == c("binomial", "poisson"))) {
      z <- if (confidence.type == "pointwise") {
        qnorm(1 - (1 - confidence.level)/2)
      } else {
        p <- length(na.omit(coef(mod)))
        scheffe(confidence.level, p)
      }
    }
    else {
      z <- if (confidence.type == "pointwise") {
        qt(1 - (1 - confidence.level)/2, df = mod$df.residual)
      } else {
        p <- length(na.omit(coef(mod)))
        scheffe(confidence.level, p, mod$df.residual)
      }
    }
    V <- vcov.(mod, complete=FALSE)
    mmat <- mod.matrix[, !is.na(mod$coefficients)] # remove non-cols with NA coeffs
    eff.vcov <- mmat %*% V %*% t(mmat)
    rownames(eff.vcov) <- colnames(eff.vcov) <- NULL
    var <- diag(eff.vcov)
    result$vcov <- eff.vcov
    result$se <- sqrt(var)
    result$se[!is.estimable] <- NA
    result$lower <- effect - z * result$se
    result$upper <- effect + z * result$se
    result$confidence.level <- confidence.level
  }
  if (is.null(transformation$link) && is.null(transformation$inverse)) {
    transformation$link <- I
    transformation$inverse <- I
  }
  result$transformation <- transformation
  result$family <- family(mod)$family
  # 2018-10-08 result$family kept to work with legacy code  
  result$link <- family(mod) 
  class(result) <- "eff"
  result
  
  
  
  
  
  
  effects::Effect.default(focal.predictors, mod, ..., sources = args)
  
  
  depth_plot$Variable<- rep("Depth", nrow(depth_plot))
  names(depth_plot)[1]<- "VarValue"
  sst_plot<- data.frame(Effect.fit_model(focal.predictors = c("Seasonal_SST"), mod = vast_fit, which_formula = "X2", xlevels = 100))
  sst_plot$Variable<- rep("Seasonal SST", nrow(sst_plot))
  names(sst_plot)[1]<- "VarValue"
  bt_plot<-  data.frame(Effect.fit_model(focal.predictors = c("Seasonal_BT"), mod = vast_fit, which_formula = "X2", xlevels = 100))
  bt_plot$Variable<- rep("Seasonal BT", nrow(bt_plot))
  names(bt_plot)[1]<- "VarValue"
  
  plot_all<- bind_rows(depth_plot, sst_plot, bt_plot)
  plot_all$Variable<- factor(plot_all$Variable, levels = c("Depth", "Seasonal SST", "Season BT"))
  
  
 
  plot(pred)
  
}

library(pdp)

# Make function to interface with pdp
pred.fun = function( object, newdata ){
  predict( x=object,
           Lat_i = object$data_frame$Lat_i,
           Lon_i = object$data_frame$Lon_i,
           t_i = object$data_frame$t_i,
           a_i = object$data_frame$a_i,
           what = "P1_iz",
           new_covariate_data = NULL,
           new_catchability_data = NULL,
           do_checks = FALSE)
}

# Run partial
Partial = partial( object = vast_fit,
                   pred.var = "Depth",
                   pred.fun = pred.fun,
                   new_covariate_data = NULL,
                   train = vast_fit$covariate_data,
                   new_catchability_data = NULL)

# Error in catchability_data when NULL and when `vast_fit$catchability_data`






for(i in seq_along(res_folders)){
  
  for(j in seq_along(fore_challenges)){
    # Get the extension...
    fore_challenge_use<- paste(fore_challenges[j], "to2019", sep = "")
    
    # Find model fit file...
    res_use<- res_files[which(grepl(fore_challenge_use, res_files))]
    
    # Load it     
    t1<- readRDS(res_use)
    
    
  }
}

Effect.fit_model<- function(focal.predictors, mod, which_formula = "X1", ...) 
  
  if(FALSE){
    mod = vast_fit
    focal.predictors = c("Depth", "SST_seasonal", "BT_seasonal")
    which_formula = "X1"
    xlevels = 100
  }
{
  if (mod$data_list$n_c > 1) {
    stop("`Effect.fit_model` is not currently designed for multivariate models")
  }
  if (!all(c("covariate_data_full", "catchability_data_full") %in% 
           ls(.GlobalEnv))) {
    stop("Please load `covariate_data_full` and `catchability_data_full` into global memory")
  }
  if (!requireNamespace("effects")) {
    stop("please install the effects package")
  }
  if (!("effects" %in% names(mod))) {
    stop("`effects` slot not detected in input to `Effects.fit_model`. Please update model using later package version.")
  }
  if (which_formula == "X1") {
    formula_orig = mod$X1_formula
    parname = "gamma1_cp"
    mod$call = mod$effects$call_X1
  } else if (which_formula == "X2") {
    formula_orig = mod$X2_formula
    parname = "gamma2_cp"
    mod$call = mod$effects$call_X2
  } else if (which_formula == "Q1") {
    formula_orig = mod$Q1_formula
    parname = "lambda1_k"
    mod$call = mod$effects$call_Q1
  } else if (which_formula == "Q2") {
    formula_orig = mod$Q2_formula
    parname = "lambda2_k"
    mod$call = mod$effects$call_Q2
  } else {
    stop("Check `which_formula` input")
  }
  whichnum = which(names(mod$parameter_estimates$par) == parname)
  mod$parhat = mod$parameter_estimates$par[whichnum]
  mod$covhat = mod$parameter_estimates$SD$cov.fixed[whichnum, whichnum, drop = FALSE]
  if (parname %in% names(mod$tmb_list$Obj$env$map)) {
    mod$parhat = mod$parhat[mod$tmb_list$Obj$env$map[[parname]]]
    mod$covhat = mod$covhat[mod$tmb_list$Obj$env$map[[parname]], mod$tmb_list$Obj$env$map[[parname]], drop = FALSE]
    mod$parhat = ifelse(is.na(mod$parhat), 0, mod$parhat)
    mod$covhat = ifelse(is.na(mod$covhat), 0, mod$covhat)
  }
  names(mod$parhat)[] = parname
  rownames(mod$covhat) = colnames(mod$covhat) = names(mod$parhat)
  formula_full = stats::update.formula(formula_orig, linear_predictor ~  . + 0)
  mod$coefficients = mod$parhat
  mod$vcov = mod$covhat
  mod$formula = formula_full
  mod$family = stats::gaussian(link = "identity")
  family.fit_model = function(x, ...) x$family
  vcov.fit_model = function(x, ...) x$vcov
  dummyfuns = list(variance = function(mu) mu, initialize = expression(mustart = y + 0.1), dev.resids = function(...) stats::poisson()$dev.res(...))
  fam = mod$family
  for (i in names(dummyfuns)) {
    if (is.null(fam[[i]])) 
      fam[[i]] = dummyfuns[[i]]
  }
  if (length(formals(fam$variance)) > 1) {
    warning("overriding variance function for effects: computed variances may be incorrect")
    fam$variance = dummyfuns$variance
  }
  args = list(call = mod$call, coefficients = mod$coefficients, 
              vcov = mod$vcov, family = fam, formula = formula_full)
  effects::Effect.default(focal.predictors, mod, ..., sources = args)
}