####  Common Resources  ####


####  Functions  ####
####
## Dynamic loading function -- each takes in a directory path and reads in files from it. 
static_covariates_read<- function(static_covariates_dir){
  
  # Get the files...
  rast_read_layers<- list.files(static_covariates_dir)
  
  # Output storage
  static_covariate_rasters_out<- vector(mode = "list", length = length(rast_read_layers))
  
  # Loop through each, save em in the list
  for(i in seq_along(rast_read_layers)){
    static_covariate_rasters_out[[i]]<- raster::stack(paste(static_covariates_dir, rast_read_layers[i], sep = "/"))
    names(static_covariate_rasters_out)[i]<- file_path_sans_ext(rast_read_layers[i])
  }
  
  # Return it
  return(static_covariate_rasters_out)
}

dynamic_covariates_read<- function(dynamic_covariates_dir){
  
  # Get the files...
  rast_read_stacks<- list.files(dynamic_covariates_dir, pattern = ".grd")
  
  # Output storage
  dynamic_covariate_rasters_out<- vector(mode = "list", length = length(rast_read_stacks))
  
  # Loop through each, save em in the list
  for(i in seq_along(rast_read_stacks)){
    dynamic_covariate_rasters_out[[i]]<- raster::stack(paste(dynamic_covariates_dir, rast_read_stacks[i], sep = "/"))
    names(dynamic_covariate_rasters_out)[i]<- file_path_sans_ext(rast_read_stacks[i])
  }
  
  # Return it
  return(dynamic_covariate_rasters_out)
}

## Processing functions
#' @title Extract static covariates at point locations
#' 
#' @description This function reads a raster layer and a sf spatial points object. The function then does a simple extraction by overlaying the spatial points onto the raster layer and gathering the value for each raster layer in the stack at all of the point locations. The function then returns the original sf_points data set as either an sf object or dataframe with columns for the raster layer as a tibble and also saves this file as "model_covs.rds".
#'
#' @param rast = A raster layer of one static covariate variable 
#' @param cov_name = The column name to use for the extracted value for the static covariate 
#' @param sf_points = SF spatial points object specifying the locations where we want to extract raster stack values
#' @param date_col_name = The column name of sf_points with the date information, formatted as YYYY-MM-DD
#' @param df_sf = Character string of "df" or "sf" signaling whether the returned object should be a dataframe or sf object 
#' 
#' @return Either an SF object or dataframe, with the covariate value appended as a new column to the point locations file. This file is also saved in out_dir. 
#' 
#' @export
static_extract<- function(rast, cov_name, sf_points, date_col_name, df_sf){
  # For debugging
  if(FALSE){
    rast = raster(paste(shared.path(os.use = os_use, group = "root", folder = "RES Data/Shapefiles/"), "NEShelf_Etopo1_bathy.tiff", sep = ""))
    cov_name = "DEPTH"
    sf_points = trawl_sf
    out_path = here::here("/scratch/aja/data/")
    df_sf = "sf"
  }
  
  # A few checks...
  # Do the projections match?
  if(!st_crs(rast) == st_crs(sf_points)){
    print("Check that projection in raster stack and spatial points match")
    stop()
  }
  
  # # Are there duplicate records?
  # if(any(duplicated(sf_points, by = c(EST_DATE, geometry)))){
  #   print("Check `sf_points` and remove duplicated observations to reduce extraction time")
  #   stop()
  # }
  
  # Extraction -----------------------------------------------------------
  sf_extract<- data.frame(raster::extract(rast, sf_points))
  names(sf_extract)<- cov_name
  
  # Bind to sf_unique 
  sf_points<- bind_cols(sf_points, sf_extract)
  
  # Return processed file -----------------------------------------------------------
  if(df_sf == "sf"){
    # Keep it as sf object
    out<- sf_points
    return(out)
  } else {
    # Drop the geometry and save the data frame
    out<- st_drop_geometry(sf_points)
    return(out)
  }
  
}

#' @title Static extract at point locations wrapper
#' 
#' @description This short wrapper function is used to process all potential covariates in the static_covariates_dir.
#'
#' @param static_covariates_list = List with static raster layers
#' @param all_tows = All tows dataframe, created by `bind_nmfs_dfo_tows()`
#' @param out_dir = Directory to save the combined dataframe as an .rds file
#' 
#' @return A datafame with information of all unique tows and appended columns for each of the static habitat covariates. This file is also saved in out_dir. 
#' 
#' @export
static_extract_wrapper<- function(static_covariates_list, sf_points, date_col_name, df_sf, out_dir){
  # For debugging
  if(FALSE){
    static_covariates_list
    sf_points
    date_col_name
    df_sf
    out_dir
  }
  
  # For each of the covariate layers in `static_covariates_list` we are going to run our `static_extract` function. I think this might require some creative manipulating of the `sf_points` so we don't lose covariates as we move through them...
  sf_points_run<- sf_points
  for(i in seq_along(static_covariates_list)){
    rast_use<- static_covariates_list[[i]]
    temp_tows_with_covs<- static_extract(rast = rast_use, cov_name = names(static_covariates_list)[[i]], sf_points = sf_points_run, date_col_name = date_col_name, df_sf = "sf")
    
    # If there are more files to go, update sf_points_run
    if(i < length(static_covariates_list)){
      sf_points_run<- temp_tows_with_covs
    }
    
    # If "last" one, good to return it
    if(i == length(static_covariates_list)){
      tows_with_covs_out<- temp_tows_with_covs
    }
  }
  
  if(df_sf == "sf"){
    # Keep it as sf object
    out<- tows_with_covs_out
    saveRDS(out, file = paste(out_dir, "all_tows_with_static_covs_sf.rds", sep = "/"))
    return(out)
  } else {
    # Drop the geometry and save the data frame
    out<- st_drop_geometry(tows_with_covs_out)
    saveRDS(out, file = paste(out_dir, "all_tows_with_static_covs.rds", sep = "/"))
    return(out)
  }
}

#' @title Extract dynamic covariates at point locations
#' 
#' @description This function summarizes spatio-temporal 2d gridded data at sample point locations. The function takes as the primary inputs a stack of one dynamic variable and a sf spatial points object. The function then does the extraction by overlaying the spatial points onto the raster stack and gathering the value for each raster layer in the stack at all of the point locations. Then, depending on the t_summ and t_position arguments, the function summarizes the fully extracted data (all points all locations) at the appropriate temporal scale and position (leading up to, saddling, or in the future) relative to the observation time. The function returns either the original sf_points data set with a new column for the dynamic raster layer values or a dataframe of the same structure and also saves this file either as "model_covs.rds" if new_file_name is NULL or as the name provided in new_file_name.
#'
#' @param rast_ts_stack = A raster stack of one covariate variable, where each layer represents a different time step
#' @param stack_name = The column name to use for the extracted value for the dynamic covariate 
#' @param t_summ = A numeric value or character string that indicates what temporal resolution should be used in summarizing the covariate values. If numeric, the function will simply summarize the values in the raster stack from the matching period back `t_summ` numeric time steps. For example, if the `rast_ts_stack` provided daily values and t_summ = 90, the function would calculate a 90 day average, where the 90-day window would either be leading up to and including the day of observation, saddled around the observation, or include the day of the observation and 89 days into the future. If a character string, should be one of "daily", monthly", "seasonal", or "annual". These options are built into the function to provide a bit easier specification to quickly calculate the monthly/seasonal/annual summaries of the raster stack values. When used, this automatically defines t_position = saddle.
#' @param t_position = A character vector of either NULL, "past", "saddle", or "future". If NULL, then values are extracted based on matching up the observation point with the dynamic raster stack at the level specified in the t_summ character vector (e.g., "daily", "monthly", "seasonal", "annual". If not, then summaries are calculated leading up to and including the observation time ("past"), saddled around around the observation time ("saddle"), or including and in the future of the observation time ("future").
#' @param sf_points = SF spatial points object specifying the locations where we want to extract raster stack values
#' @param date_col_name = The column name of sf_points with the date information, formatted as YYYY-MM-DD
#' @param df_sf = Character string of "df" or "sf" signaling whether the returned object should be a dataframe or sf object 
#' 
#' @return Either an SF object or dataframe, with the covariate values appended as new columns to the point locations file. This file is also saved in out_dir. 
#' 
#' @export

dynamic_2d_extract<- function(rast_ts_stack, stack_name, t_summ, t_position, sf_points, date_col_name, df_sf){
 
  # Custom rowMeans like function
  rowMeans_cust<- function(row_id, start_col, end_col, full_data){
    out<- rowMeans(full_data[row_id,start_col:end_col])
    return(out)
  }
  
  # For debugging
  if(FALSE){
    rast_ts_stack = raster::stack(paste(shared.path(os.use = "unix", group = "root", folder = "RES_Data/OISST/"), "ThroughFeb2020.grd", sep = ""))
    # Annoyingly, dates not preserved...
    rast_ts_dates<- seq(from = ymd('1981-09-01'), to = ymd('1981-09-01') + nlayers(rast_ts_stack)-1, by = 'day')
    names(rast_ts_stack)<- rast_ts_dates
    writeRaster(projectRaster(rast_ts_stack, crs = 4326), filename = here::here("scratch/aja/targets_flow/data/covariates/dynamic/SST.grd"), overwrite = TRUE, format = "raster")
    t<- raster::stack(here::here("scratch/aja/targets_flow/data/covariates/dynamic/SST.grd"))
    stack_name = "SST"
    t_summ = "seasonal"
    t_position = NULL
    sf_points = trawl_covs
    out_path = here::here("/scratch/aja/data/")
    new_file_name = NULL
    
    rast_ts_stack = oisst_stack
    stack_name = "sst"
    t_summ = "seasonal"
    t_position = NULL
    sf_points = cov_sf
    df_sf = df
    out_path = here::here("scratch/aja/data")
  }
  
  # A few checks...
  # If t_summ is a character, does it match one of "daily", monthly", "seasonal", or "annual" AND is t_position set to NULL?
  if(is.character(t_summ)){
    t_summ_check<- t_summ %in% c("monthly", "seasonal", "annual") & is.null(t_position)
    if(!t_summ_check){
      print("Check 't_summ' argument and 't_position'. 't_summ' must be one of 'daily' monthly', 'seasonal' or 'annual' and 't_position' = NULL")
      stop()
    }
  }
  
  # Check t_position argument
  t_position_check<- is.null(t_position) || t_position %in% c("past", "saddle", "future") 
  if(!t_position_check){
    print("Check 't_position' argument. Must be one of 'past', 'saddle' or 'future'")
    stop()
  }
  
  # # Are there duplicate records?
  # if(any(duplicated(sf_points, by = c(EST_DATE, geometry)))){
  #   print("Check `sf_points` and remove duplicated observations to reduce extraction time")
  #   stop()
  # }
  
  # Finally, do the projections match?
  if(!st_crs(rast_ts_stack) == st_crs(sf_points)){
    print("Check that projection in raster stack and spatial points match")
    stop()
  }
  
  # To do the matching for seasons, need some type of look up table.
  month_season_table<- data.frame("Month" = str_pad(seq(from = 1, to = 12, by = 1), 2, "left", 0), "Season" = c("Winter", "Winter", "Spring", "Spring", "Spring", "Summer", "Summer", "Summer", "Fall", "Fall", "Fall", "Winter"))
  
  # Rename the date column to be "EST_DATE" to match the function set up
  sf_points<- sf_points %>%
    rename(., "EST_DATE" := {{date_col_name}})

  # Full extraction, all points and layers -----------------------------------------------------------
  sf_extract<- data.frame(raster::extract(rast_ts_stack, sf_points))
  
  # Getting the values we want from full extraction -------------------------
  # This leaves us with a data frame where the rows are the sf_points and then the columns are the covariate value for a given location for each of the layers (i.e., time steps) in the stack. To start with calculating any summary, we first need to know the which column provides an exact match between the time of the observation and the time in the raster stack. Then, depending on if it is a past, saddle, or future, and how many columns are specified by t_pos, we can average the necessary columns.
  colnames(sf_extract)<- gsub("[.]", "-", gsub("X", "", colnames(sf_extract)))
  
  # Now, accounting for the different windows and positions. As mentioned earlier, this is going to be completely different depending on if t_summ is a character sting ("daily", "monthly", "seasonal", "annual") or if it is a numeric integer. So, here comes that split...
  summ_df_list<- vector("list", length(t_summ))
  
  #####
  ## t_summm as a character string
  #####
  if(is.character(t_summ)){
    
    for(i in seq_along(t_summ)){
      # Get t_summ for this iteration
      t_summ_use<- t_summ[i]
      
      # Create a data frame to store the start and end...
      summ_df<- data.frame("Point" = 1:nrow(sf_extract), "Date_Match" = rep(NA, nrow(sf_extract)), "Start_Summ" = rep(NA, nrow(sf_extract)), "End_Summ" = rep(NA, nrow(sf_extract)))
      
      # Get the exact date match
      summ_df$Date_Match<- match(sf_points$EST_DATE, as.Date(colnames(sf_extract)))
      
      # Column index start based on t_summ_use. Seasonal needs more finagling.
      if(t_summ_use != "seasonal"){
        summ_df$Start_Summ<- switch(t_summ_use,
                                    "daily" = summ_df$Date_Match,
                                    "monthly" = match(format(sf_points$EST_DATE, "%Y-%m"), format(rast_ts_match, "%Y-%m")),
                                    "annual" = match(format(sf_points$EST_DATE, "%Y"), format(rast_ts_match, "%Y")))
        
        # Now index end
        summ_df$End_Summ<- switch(t_summ_use,
                                  "daily" = summ_df$Date_Match,
                                  "monthly" = sapply(format(sf_points$EST_DATE, "%Y-%m"), FUN = function(x) max(which(x == format(rast_ts_match, "%Y-%m")))),
                                  "annual" = sapply(format(sf_points$EST_DATE, "%Y"), FUN = function(x) max(which(x == format(rast_ts_match, "%Y")))))
      } else {
        # Season bits, need a year - season combo for both the points AND columns of sf_extract.
        colnames_orig<- as.Date(colnames(sf_extract))
        colnames_season<- month_season_table$Season[match(format(colnames_orig, "%m"), month_season_table$Month)]
        colnames_season_match<- paste(format(colnames_orig, "%Y"), colnames_season, sep = "-")
        
        sf_points$Season_Match<- paste(format(sf_points$EST_DATE, "%Y"), month_season_table$Season[match(format(sf_points$EST_DATE, "%m"), month_season_table$Month)], sep = "-")
        
        # Start index
        summ_df$Start_Summ<- match(sf_points$Season_Match, colnames_season_match)
        
        # Now index end
        summ_df$End_Summ<- sapply(sf_points$Season_Match, FUN = function(x) max(which(x == colnames_season_match)))
      }
      
      # Deal with the "-Inf" indices, arising when there are no matches for "max"
      summ_df[summ_df == -Inf]<- NA
      
      # Store it
      summ_df_list[[i]]<- summ_df
      names(summ_df_list)[i]<- paste(stack_name, t_summ_use, sep = "_")
    }
  } 
  
  #####
  ## t_summ as a numeric vector
  #####
  if(is.numeric(t_summ)){
    for(i in seq_along(t_summ)){
      # Create a data frame to store the start and end...
      summ_df<- data.frame("Point" = 1:nrow(sf_extract), "Date_Match" = rep(NA, nrow(sf_extract)), "Start_Summ" = rep(NA, nrow(sf_extract)), "End_Summ" = rep(NA, nrow(sf_extract)))
      
      # Get the exact date match
      summ_df$Date_Match<- match(sf_points$EST_DATE, as.Date(colnames(sf_extract)))
      
      # Get summary window time range
      t_summ_use<- t_summ[i]
      
      # Column index start based on the t_position and t_summ_use
      summ_df$Start_Summ<- switch(t_position,
                                  "past" = summ_df$Date_Match - (t_summ_use),
                                  "saddle" = summ_df$Date_Match - (t_summ_use/2),
                                  "future" = summ_df$Date_Match)
      # Now index end
      summ_df$End_Summ<- switch(t_position,
                                "past" = summ_df$Date_Match,
                                "saddle" = summ_df$Date_Match + (t_summ_use/2),
                                "future" = summ_df$Date_Match + (t_summ_use))
      
      # Store it
      summ_df_list[[i]]<- summ_df
      names(summ_df_list)[i]<- paste(stack_name, t_summ_use, sep = "_")
    }
  }
  
  # Subsetting full extraction, taking row means across the columns signaled by Start_Summ and End_Summ. First, need to drop any NAs as those are going to through an error in the indexing AND any negative index numbers (this would indicate a situation where the date of observation falls during the beginning of the dynamic variable time series, such that there isn't data for the prior or saddle option t_summ_use time steps prior to the observation)
  # Has to be a better option for how to do this, but a loop for now..
  # Adding Point ID column, used for merging
  sf_points$Point<- seq(from = 1, to = nrow(sf_points))
  
  for(j in seq_along(summ_df_list)){
    
    summ_df_comp<- summ_df_list[[j]] %>%
      drop_na() %>%
      dplyr::filter(., Start_Summ > 0)
    
    # Calculate mean given start and end column index of extraction file and row based on observation point ID
    sf_extract_mean<- summ_df_comp %>% 
      mutate(., "Summ_Val" = pmap_dbl(list(row_id = Point, start_col = Start_Summ, end_col = End_Summ, full_data = list(sf_extract)), rowMeans_cust))
    
    # Don't need to keep all the columns, just point and summ_val. Rename to match stack and t_summ
    point_mean<- sf_extract_mean %>%
      dplyr::select(., Point, Summ_Val)
    names(point_mean)[2]<- names(summ_df_list)[j]
    
    # Join back to original trawl sf_points
    sf_points<- sf_points %>%
      left_join(., point_mean, by = c("Point" = "Point")) %>%
      dplyr::select(., -Point)
  }
  
  # Save and return processed file -----------------------------------------------------------
  if(df_sf == "sf"){
    # Keep it as sf object
    out<- sf_points
    return(out)
  } else {
    # Drop the geometry and save the data frame
    out<- st_drop_geometry(sf_points)
    return(out)
  }
  #########
  ## End
  #########
}

#' @title Extract dynamic covariates at point locations wrapper
#' 
#' @description This short wrapper function is used to process all potential covariates in the dynamic_covariates_dir.
#'
#' @param dynamic_covariates_list = 
#' @param t_summ
#' @param t_position 
#' @param sf_points
#' @param date_col_name
#' @param df_sf
#' @param out_dir = Directory to save the combined dataframe as an .rds file
#' 
#' @return A datafame with information of all unique tows with dynamic covariate values as new columns. This file is also saved in out_dir. 
#' 
#' @export
dynamic_2d_extract_wrapper<- function(dynamic_covariates_list, t_summ, t_position, sf_points, date_col_name, df_sf, out_dir){
  # For debugging
  if(FALSE){
    dynamic_covariates_dir
    stack_names
    t_summ
    t_position
    sf_points
    date_col_name
    df_sf
    out_dir
  }
  
  # For each of the covariate layers in `dynamic_covariates_list` we are going to run our `dynamic_2d_extract` function. I think this might require some creative manipulating of the `sf_points` so we don't lose covariates as we move through them...
  sf_points_run<- sf_points
  for(i in seq_along(dynamic_covariates_list)){
    stack_use<- dynamic_covariates_list[[i]]
    temp_tows_with_covs<- dynamic_2d_extract(rast_ts_stack = stack_use, stack_name = names(dynamic_covariates_list)[[i]], t_summ = "seasonal", t_position = NULL, sf_points = sf_points_run, date_col_name = date_col_name, df_sf = "sf")

    # If there are more files to go, update sf_points_run
    if(i < length(dynamic_covariates_list)){
      sf_points_run<- temp_tows_with_covs
    }
    
    # If "last" one, good to return it
    if(i == length(dynamic_covariates_list)){
      tows_with_covs_out<- temp_tows_with_covs
    }
  }
  
  if(df_sf == "sf"){
    # Keep it as sf object
    out<- tows_with_covs_out
    saveRDS(out, file = paste(out_dir, "all_tows_all_covs_sf.rds", sep = "/"))
    return(out)
  } else {
    # Drop the geometry and save the data frame
    out<- st_drop_geometry(tows_with_covs_out)
    saveRDS(out, file = paste(out_dir, "all_tows_all_covs.rds", sep = "/"))
    return(out)
  }
}