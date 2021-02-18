#### Functions for cmip6 processing  ####
# This file should be sourced to provide the stepping stones for the sdm_workflows pipeline
# Used to either facilitate a {targets} workflow or to build a makefile workflow.


####  Common Resources  ####

# Common cropping area for all the netcdf files used in this project
study_area <- extent(c(260, 320, 20, 70))

# Common month abbreviations
month_abbrevs <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")

# Months as 0-padded numerics
months_numeric <- str_pad(seq(from = 1, to = 12, by = 1), width = 2, side = "left", pad = 0)


####  OISST Processing  ####


#' @title Import OISST Climatology
#' 
#' @description Load OISST climatology from box as a raster stack
#'
#' @param climatology_period Choices: "1991-2020" or "1982-2011"
#' @param os.use Specification for gmRi::shared.path for Mac/Windows paths
#'
#' @return Rster stack of OISST Daily Climatology cropped to study area

import_oisst_clim <- function(climatology_period = "1991-2020", os.use = "unix"){
  
  # Path to OISST on Box
  oisst_path <- shared.path(os.use = os.use, group = "RES_Data", folder = "OISST/oisst_mainstays")
  
  # re-direct based on desired climatology
  if(climatology_period == "1991-2020"){
    clim_stack <- stack(paste0(oisst_path, "daily_climatologies/daily_clims_1991to2020.nc"))
  
    }else if(climatology_period == "1982-2011"){
      clim_stack <- stack(paste0(oisst_path, "daily_climatologies/daily_clims_1982to2011.nc"))}
  
  # Crop it to study area
  clim_cropped <- crop(clim_stack, study_area)
  return(clim_cropped)
}




#' @title Load CMIP6 NetCDF data.
#' 
#' @description Load CMIP6 data off BOX by path to specific file(s) in RES_Data/CMIP6/.
#' Default loads the CMIP6 SST subset used for testing. 
#'
#' @param cmip_file Either "tester", or path to cmip stack within RES_Data/CMIP6
#'
#' @return Raster stack of CMIP Data, cropped to study area.
#'
import_cmip_sst_test <- function(cmip_file = "tester"){
  
  # General path to all the cmip data on Box
  cmip_path    <- shared.path(os.use = "unix", group = "RES_Data", folder = "CMIP6/")
  
  # Load the stack(s) cropped to the study area
  if(cmip_file == "tester"){
    message(paste0("Loading CMIP6 File: tos_Omon_CanESM5_historical_r1i1p2f1_gn_195501-201412.nc.1x1.nc"))
    cmip_full   <- stack(paste0(cmip_path, "TestFiles/tos_Omon_CanESM5_historical_r1i1p2f1_gn_195501-201412.nc.1x1.nc"))
    cmip_cropped <- crop(cmip_full, study_area)
  
    } else if(length(cmip_file) == 1){
        message(paste0("Loading CMIP6 File: ", cmip_file))
        cmip_full <- stack(paste0(cmip_path, cmip_file))
        cmip_cropped <- crop(cmip_full, study_area)
  
      } else if(length(cmip_file > 1)){
          cmip_cropped <- map(cmip_file, function(x){
            message(paste0("Loading Multiple CMIP Files, Returning List"))
            cmip_full   <- stack(paste0(cmip_path, x))
            cmip_cropped <- crop(cmip_full, study_area)}) %>% 
            setNames(cmip_file)
  }
  
  
  return(cmip_cropped)
  
}



#' @title Monthly Mean from Climatology
#' 
#' @details Uses the names of each daily layer in the OISST climatology pattern
#' to group the data by month and return a monthly climatology.
#'
#' @param clim_source Raster stack for climatology you want re-grouped to monthly
#' @param month_layer_key If null uses OISST key. Should be named list indicating which layers
#' belong to which group for the desired output. Ex. January is layers 1:31 in OISST Climatology,
#' month_layer_key <- list("Jan" = c(1:31), ...)
#'
#' @return Raster stack of means from original climatology source using new temporal groups.
#'
months_from_clim <- function(clim_source, month_layer_key = NULL){
  
  
  ####  Default Key to get months from Modified ordinal Day
  
  # Uses number key to match months to corresponding the day of year,
  # day of year in climatology honors the 60th as only feb 29th
  # in this system march 1st is always 61, and Dec 31st is always 366
  
  if(is.null(month_layer_key)){
    mod_months <- list(
      "Jan" = c(1:31),    "Feb" = c(32:60),   "Mar" = c(61:91),
      "Apr" = c(92:121),  "May" = c(122:152), "Jun" = c(153:182),
      "Jul" = c(183:213), "Aug" = c(214:244), "Sep" = c(245:274),
      "Oct" = c(275:305), "Nov" = c(306:335), "Dec" = c(336:366))
    
    # Put a capital X in front of each number to  match raster stack names
    mod_months <- map(mod_months, ~ str_c("X", .x))
    message("Using default key for modified ordinal days 0-366")
    month_layer_key <- mod_months} 
  
  
  # Map through the months to get mean of each one from climatology
  monthly_avgs <- map(month_layer_key, function(layer_indices){
      month_mean <- mean(clim_source[[layer_indices]])}) %>%
    stack() %>%
    setNames(names(month_layer_key))
  
  #return the monthly average stack
  return(monthly_avgs)
 }


####  SODA Processing  ####




#' @title Import SODA Monthly Climatology
#' 
#' @description Load raster stack of SODA climatology for desired variable. Choices
#' are "surf_sal", "surf_temp", "bot_sal", "bot_temp". Area is also cropped to study area.
#'
#' @param var_name variable name to use when stacking data
#' @param os.use windows mac toggle for box path
#'
#' @return Raster stack for monthly climatology, cropped to study area
#' @export
#'
#' @examples
import_soda_clim <- function(var_name = c("surf_sal", "surf_temp", "bot_sal", "bot_temp"), 
                             os.use = "unix"){
  
  # Path to OISST on Box
  soda_path <- shared.path(os.use = os.use, group = "RES_Data", folder = "SODA")
  
  # Load climatology for desired variable
  clim_stack <- stack(paste0(soda_path, "SODA_monthly_climatology1990to2019.nc"), varname = var_name)
  
  # Crop it to study area
  clim_cropped <- crop(clim_stack, study_area)
  return(clim_cropped)
}













#### CMIP Processing  ####


#' @title Get Climate Reference from CMIP6
#' 
#' @description Return monthly climatology from CMIP6 Data for specified reference period.
#'
#' @param cmip_stack Raster stack of CMIP6 Data Variable
#' @param clim_years Length 2 vector of start and end year for reference period
#'
#' @return Raster stack of monthly climatology for input CMIP6 Variable
#'
cmip_to_clim <- function(cmip_stack = cmip_cropped, clim_years = NULL){
  
  # Use 1982:2011 as default
  if(is.null(clim_years)){
    clim_years <- as.character(c(1991:2020))
    message("Using 1991-2020 as Climate Reference Period")
  
    } else {
      clim_years <- as.character(clim_years)
      message(paste0("Using ", clim_years[1], "-", clim_years[2], " as Climate Reference Period"))}
  
  
  # Pull out the names of the cmip layers for matching
  cmip_layers <- names(cmip_stack)
  
  # Pull layers of cmip data for the years of interest for climatology
  cmip_clim_years  <- cmip_stack[[which(str_sub(cmip_layers, 2,5) %in% clim_years)]]
  clim_year_layers <- names(cmip_clim_years)
  
  # Get strings to match the months, set names as abbreviations for later
  month_labels <- str_pad(c(1:12), width = 2, pad = "0", "left")
  month_labels <- setNames(month_labels, month_abbrevs)
  
  # Loop through the months, getting mean across the climatology period
  cmip_clim <- map(month_labels, function(month_index){
    
    # Indices for the month
    cmip_month_indices <- which(str_sub(clim_year_layers, 7,8) == month_index)
    
    # Mean across years
    monthly_clim <- mean(cmip_clim_years[[cmip_month_indices]])
    return(monthly_clim)
  }) %>% stack()
  
  # Return the new climatology stack
  return(cmip_clim)

}





#' @title CMIP Monthly Anomalies from Climatology
#'
#' @param cmip_data CMIP Observations
#' @param cmip_clim  Matching Climatology for a given reference period
#'
#' @return Raster stack of dimensions equal to cmip_data but with anomaly values
#' 
cmip_get_anomalies <- function(cmip_data, cmip_clim){
  
  
  # Map though the month key to pull their data, return quantile stack
  month_index <- str_pad(c(1:12), 2, "left", "0")
  
  # Map through the months
  # month_index are how they are as dates in Xyyyy.mm.dd
  # month abbrevs are how they are named in climatology
  monthly_anoms <- map2(month_index, month_abbrevs, function(month_index, month_abbrev){
    
    # Use month abbreviation to get month from climatology
    clim_month_data <- cmip_clim[[month_abbrev]]
    
    # Use month_index for CMIP layers to match the month
    cmip_layers <- names(cmip_data)
    cmip_month_indices <- which(str_sub(cmip_layers, 7, 8) == month_index)
    cmip_month_labels  <- cmip_layers[cmip_month_indices]
    
    # Pull all layers from the cropped cmip stack with that month
    month_layers <- cmip_data[[cmip_month_indices]]
    
    # Subtract climate average to get anomalies
    month_anoms <- month_layers - clim_month_data
    month_anoms <- setNames(month_anoms, cmip_month_labels)
    
    
  }) %>% stack()
  
  return(monthly_anoms)
  
  
}



####  Bias Correction  Functions  ####


#' @title Resample Resolution to Destination Grid Res
#' 
#' @description Takes two raster stacks, resamples starting_grid to resolution of
#' desired_grid. Default method is bilinear interpolation.
#'
#' @param starting_grid Raster stack whose desired resolution you wish to change
#' @param desired_grid Raaster stack with the desired output resolution to use as example grid
#' @param method Optional string indicating method to use for raster::resample(method = method)
#'
#' @return Raster stack with resolution of desired_grid and data from starting_grid
#' 
resample_grid <- function(starting_grid = cmip_anoms, 
                          desired_grid = oisst_month_avgs, 
                          method = "bilinear"){
  resample(starting_grid, desired_grid[[1]], method = method)
}





#' @title Bias Correction - Delta Method
#' 
#' @description Bias correct climate data using a reference climatology. Performs the delta-method 
#' where anomalies in climate data are applied directly to reference climate.
#'
#' @param cmip_grid Input Grid we want to bias correct the observations for
#' @param reference_climatology The climatology for real-life data to bias correct with
#'
#' @return Raster stack of bias corrected observations
#' 
delta_method_bias_correct <- function(cmip_grid = cmip_anoms_regridded, 
                                      reference_climatology = oisst_month_avgs_91){
  
  ####  Objective:
  # Add CMIP monthly anomalies directly to the climate reference monthly averages
  
  
  #### Do Some work with the names so we can find matches between the cmip anoms and the reference climatology 
  
  # change oisst names to the numeric ones to streamline the matching
  reference_climatology <- setNames(reference_climatology, months_numeric)
  
  # Get dates from cmip in non-raster format
  cmip_dates  <- as.Date(gsub("[.]", "-", gsub("X", "", names(cmip_grid))))
  
  # For each time step, add the anomalies to the
  # Alright, apply the anomalies (deltas) of the climate model to OISST climatology. The code below "loops" over each of the layers in the climate model anomalies stack and then adds those values to the matching oisst monthly climatology bsased on the month of the climate model anomalies stack
  cmip_proj_out <- lapply(seq(1:nlayers(cmip_grid)), function(layer_index) { 
    
    # Grab the first CMIP layer
    cmip_anom_layer <- cmip_grid[[layer_index]]
    
    # Get the corresponding layer number for the reference climatology
    
    # Pull month from Z, use to match reference climatology names
    # month_digits   <- format(getZ(cmip_grid)[layer_index], "%m")
    # layer_match    <- match(month_digits, gsub("X", "", names(reference_climatology)))
    
    month_digits   <- str_sub(names(cmip_anom_layer)[1], 7,8)
    layer_match    <- which(str_detect(names(reference_climatology), month_digits))
    ref_clim_layer <- reference_climatology[[layer_match]]
    
    # Add ref clim to anomalies to bias correct them, returns bias-corrected temperature
    delta_out <- cmip_anom_layer + ref_clim_layer}) %>% 
      stack()
  
  # fix those names for consistency
  names(cmip_proj_out) <- cmip_dates
  return(cmip_proj_out)
  
  
}











####  Processing Mean/5th/95th  ####

# Re-stack all the bias-corrected cmip datasets by time step
# stores years in a list


#' @title Re-Stack CMIP Delta Bias Corrected Data
#' 
#' @description Takes bias corrected data sources and re-stacks them to align on the same 
#' time steps. This sets up the stacks for assessing 5th and 95th percentile and mean data.
#'
#' @param cmip_inputs 
#' 
#'
#' @return
#' 
restack_cmip_projections <- function(cmip_inputs = cmip_delta_bias_corrected){
  
  # Get the number of total time steps from a given stack
  
  
  # map around that length and pull out the layers of each
  # Put them in stacks by time step, try to keep the names of their sources?
  
  
}







#' @title Raster Quantiles from Timestep Stacks
#' 
#' @description Get the Mean/5th/95th percentile at each time step from an ensemble of climate
#' projections. Typically used with map() to iterate through years. This function operates on
#' the name structure of the months.
#'
#' @param year_stacks Raster stack containing a full year of each CMIP model's data
#'
#' @return
#' 
timestep_stats <- function(year_stacks){
  
  # Map though the month key to pull their data, return quantile stack
  month_key <- str_pad(c(1:12), 2, "left", "0")
  month_labels <- paste0("X", month_key)
  
  
  # Get stack containing mean values
  monthly_percentiles_05 <- map(month_key, function(month_index){
    
    # What layers match the month
    raster_layers <- str_sub(names(year_stacks), 7, 8) #month letters in X2020.01.01
    which_days    <- which(str_detect(raster_layers, month_index) == TRUE)
    
    # Pull the days out
    month_subset_ras <- year_stacks[[which_days]]
    
    # Use calc to get mean + percentiles
    ras.quant <- calc(month_subset_ras, mean,  na.rm = T)
    
  }) %>% 
    setNames(month_labels) %>% 
    stack()
  
  # Get stack containing 5th percentile values
  monthly_percentiles_mean <- map(month_key, function(month_index){
    
    # What layers match the month
    raster_layers <- names(year_stacks)
    which_days <- which(str_detect(raster_layers, month_index) == TRUE)
    
    # Pull the days out
    month_subset_ras <- year_stacks[[which_days]]
    
    # Use calc to get mean + percentiles
    ras.05   <- calc(month_subset_ras, quantile, probs = 0.05, na.rm = T)
    
  }) %>% 
    setNames(month_labels) %>% 
    stack()
  
  # 95th percentile stack
  monthly_percentiles_95 <- map(month_key, function(month_index){
    
    # What layers match the month
    raster_layers <- names(year_stacks)
    which_days <- which(str_detect(raster_layers, month_index) == TRUE)
    
    # Pull the days out
    month_subset_ras <- year_stacks[[which_days]]
    
    # Use calc to get mean + percentiles
    ras.95   <- calc(month_subset_ras, quantile, probs = 0.95, na.rm = T)
  }) %>% 
    setNames(month_labels) %>% 
    stack()
  
  # For each year return a list of each month for the three quantiles
  year_quants <- list(
    "percentile_05" = monthly_percentiles_05,
    "mean"          = monthly_percentiles_mean,
    "percentile_95" = monthly_percentiles_95)
  
  
  return(year_quants)
  
  
}

# Function to re-stack the different scenarios as singular timelines by stat/quantile:
#takes the list of raster stacks, organized by year, and the corresponding names as a vector or 
# by passing the named list to imap()


#' @title Reassemble Timeseries from List of Ensemble Quantiles
#' 
#' @description 
#'
#' @param year_stack Output Raster stack or list of raster stacks from timestep_stats()
#' @param year_lab Matching string indicating the year or timestep to go with each list item
#' @param stat_group String indicating which statistic to extract ("mean", "percentile05", 
#' "percentile95")
#'
#' @return
#' 
timestep_to_full <- function(year_stack, year_lab, stat_group){
  stat_out   <- year_stack[[stat_group]] 
  orig_names <- str_replace(names(stat_out), "X", "")
  new_names  <- paste0(paste0("X", year_lab, "_"), orig_names)
  stat_out   <- setNames(stat_out, new_names)
}