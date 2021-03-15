# Beginning of _targets.R

# Load packages and set specific options for the workflow
library(targets)
library(tarchetypes)
library(sp)
library(raster)
library(tidyverse)
source("CMIP6_processing/R/sdm_workflow_funs.R")
options(tidyverse.quiet = T)


# Additional packages to load specific to target(s)
tar_option_set( packages = c("raster", "sf", "rmarkdown", "tidyverse", "gmRi") )






# Define target pipeline: Outlines high-level steps of the analysis
# Format is just a list of all the targets
# Order is not important, package sorts out connections for everything
list(
 
   #### CMIP Import  ####
  
  # Surface Temperature
  tar_target(
    name = CMIP_SurfTemp_collection,
    command = import_cmip_collection(cmip_var = "surf_temp") %>% map(get_cmip_dates)),
  tar_target(
    name = CMIP_SurfTemp_clim,
    command = cmip_to_clim(cmip_stack = CMIP_SurfTemp_collection) ),
  tar_target(
    name = CMIP_SurfTemp_anoms,
    command = cmip_get_anomalies(cmip_data = CMIP_SurfTemp_collection, 
                                 cmip_clim = CMIP_SurfTemp_clim) ),
  
  # Surface Salinity
  tar_target(
    name = CMIP_SurfSal_collection,
    command = import_cmip_collection(cmip_var = "surf_sal") %>% map(get_cmip_dates)),
  tar_target(
    name = CMIP_SurfSal_clim,
    command = cmip_to_clim(cmip_stack = CMIP_SurfSal_collection) ),
  tar_target(
    name = CMIP_SurfSal_anoms,
    command = cmip_get_anomalies(cmip_data = CMIP_SurfSal_collection, 
                                 cmip_clim = CMIP_SurfSal_clim) ),
  
  # Bottom Temperature
  tar_target(
    name = CMIP_BotTemp_collection,
    command = import_cmip_collection(cmip_var = "bot_temp") %>% map(get_cmip_dates)),
  tar_target(
    name = CMIP_BotTemp_clim,
    command = cmip_to_clim(cmip_stack = CMIP_BotTemp_collection) ),
  tar_target(
    name = CMIP_BotTemp_anoms,
    command = cmip_get_anomalies(cmip_data = CMIP_BotTemp_collection, 
                                 cmip_clim = CMIP_BotTemp_clim) ),
  
  # Bottom Salinity
  tar_target(
    name = CMIP_BotSal_collection,
    command = import_cmip_collection(cmip_var = "bot_sal") %>% map(get_cmip_dates)),
  tar_target(
    name = CMIP_BotSal_clim,
    command = cmip_to_clim(cmip_stack = CMIP_BotSal_collection) ),
  tar_target(
    name = CMIP_BotSal_anoms,
    command = cmip_get_anomalies(cmip_data = CMIP_BotSal_collection, 
                                 cmip_clim = CMIP_BotSal_clim) ),
  
  
  
  #### OISST Import  ####
  tar_target(
    name = OISST_clim_85,
    command = import_oisst_clim(climatology_period = "1985-2014", 
                                os.use = "unix") ),
  tar_target(
    name = OISST_monthly_clim,
    command = months_from_clim(clim_source = OISST_clim_85, 
                               month_layer_key = NULL) ),
  
  
  
  ####  SODA Import  ####
  
  # # load data - make climatologies
  # tar_target(
  #   name = soda_cropped,
  #   command = map(c("surf_sal", "surf_temp", "bot_sal", "bot_temp"), import_soda, os.use = "unix")),
  # tar_target(
  #   soda_climatologies,
  #   command = "Run soda_climatology.py"),
  
  # load climatologies
  tar_target(
    name = SODA_SurfSal_clim,
    command = import_soda_clim("surf_sal", os.use = "unix")),
  tar_target(
    name = SODA_BotTemp_clim,
    command = import_soda_clim("bot_temp", os.use = "unix")),
  tar_target(
    name = SODA_BotSal_clim,
    command = import_soda_clim("bot_sal", os.use = "unix")),
  
  
  #### Bias Correction  ####
  
  ####__  SST  ####
  tar_target(
    name = SurfTemp_bc,
    command =  map(CMIP_SurfTemp_anoms, function(anom_grid){
          
        # check for problem data
        if(class(anom_grid) == "character"){
          return("Problem with CMIP Naming Structure")}
        
        # run for data that passes check
        delta_method_bias_correct(cmip_grid = anom_grid, 
                                  reference_climatology = OISST_monthly_clim)})
    ),
  
  ####__  Surface Salinity  ####
  tar_target(
    name = SurfSal_bc,
    command = map(CMIP_SurfSal_anoms, function(anom_grid){
      
      # check for problem data
      if(class(anom_grid) == "character"){
        return("Problem with CMIP Naming Structure")}
      
      # run for data that passes check
      delta_method_bias_correct(cmip_grid = anom_grid, 
                                reference_climatology = SODA_SurfSal_clim)})
    ),
  
  ####__  Bottom Temperature  ####
  tar_target(
    name = BotTemp_bc,
    command = map(CMIP_BotTemp_anoms, function(anom_grid){
      
      # check for problem data
      if(class(anom_grid) == "character"){
        return("Problem with CMIP Naming Structure")}
      
      # run for data that passes check
      delta_method_bias_correct(cmip_grid = anom_grid, 
                                reference_climatology = SODA_BotTemp_clim)})
  ),
  
  ####__  Bottom Salinity  ####
  tar_target(
    name = BotSal_bc,
    command = map(CMIP_BotSal_anoms, function(anom_grid){
      
      # check for problem data
      if(class(anom_grid) == "character"){
        return("Problem with CMIP Naming Structure")}
      
      # run for data that passes check
      delta_method_bias_correct(cmip_grid = anom_grid, 
                                reference_climatology = SODA_BotSal_clim)})
  )
  
  
  #### Pull Quantiles ####
  
  
  # Up to here on the targets part
  
  # # Assemble single variable stacks of all climate projections
  # tar_target(
  #   name = unbiased_SurfTemp_ensemble,
  #   command = restack_cmip_projections(cmip_inputs = SurfTemp_bias_corrected) ),
  # tar_target(
  #   name = unbiased_SurfSal_ensemble,
  #   command = restack_cmip_projections(cmip_inputs = SurfSal_bias_corrected) ),
  # tar_target(
  #   name = unbiased_BotTemp_ensemble,
  #   command = restack_cmip_projections(cmip_inputs = BotTemp_bias_corrected) ),
  # tar_target(
  #   name = unbiased_BotSal_ensemble,
  #   command = restack_cmip_projections(cmip_inputs = BotSal_bias_corrected) ),
  
  
  
  
  
#   # Get single variable quantiles
#   tar_target(
#     name = SurfTemp_quants,
#     command = timestep_stats(year_stacks = unbiased_SurfTemp_ensemble) ),
#   tar_target(
#     name = SurfSal_quants,
#     command = timestep_stats(year_stacks = unbiased_SurfSal_ensemble) ),
#   tar_target(
#     name = BotTemp_quants,
#     command = timestep_stats(year_stacks = unbiased_BotTemp_ensemble) ),
#   tar_target(
#     name = BotSal_quants,
#     command = timestep_stats(year_stacks = unbiased_BotSal_ensemble) ),
#   
#   
#   
#   
#   # Pull out single variable means/5th/95th
#   
#   # Surface Temperature
#   tar_target(
#     name = SurfTemp_mean,
#     command = imap(SurfTemp_quants, timestep_to_full, stat_group = "mean") %>% stack() ),
#   tar_target(
#     name = SurfTemp_05,
#     command = imap(SurfTemp_quants, timestep_to_full, stat_group = "percentile_05") %>% stack() ),
#   tar_target(
#     name = SurfTemp_95,
#     command = imap(SurfTemp_quants, timestep_to_full, stat_group = "percentile_95") %>% stack()),
#   
#   # Surface Salinity
#   tar_target(
#     name = SurfSal_mean,
#     command = imap(SurfSal_quants, timestep_to_full, stat_group = "mean") %>% stack() ),
#   tar_target(
#     name = SurfSal_05,
#     command = imap(SurfSal_quants, timestep_to_full, stat_group = "percentile_05") %>% stack() ),
#   tar_target(
#     name = SurfSal_95,
#     command = imap(SurfSal_quants, timestep_to_full, stat_group = "percentile_95") %>% stack()),
#   
#   # Bottom Temperature
#   tar_target(
#     name = BotTemp_mean,
#     command = imap(BotTemp_quants, timestep_to_full, stat_group = "mean") %>% stack() ),
#   tar_target(
#     name = BotTemp_05,
#     command = imap(BotTemp_quants, timestep_to_full, stat_group = "percentile_05") %>% stack() ),
#   tar_target(
#     name = BotTemp_95,
#     command = imap(BotTemp_quants, timestep_to_full, stat_group = "percentile_95") %>% stack()),
#   
#   # Bottom Salinity
#   tar_target(
#     name = BotSal_mean,
#     command = imap(BotSal_quants, timestep_to_full, stat_group = "mean") %>% stack() ),
#   tar_target(
#     name = BotSal_05,
#     command = imap(BotSal_quants, timestep_to_full, stat_group = "percentile_05") %>% stack() ),
#   tar_target(
#     name = BotSal_95,
#     command = imap(BotSal_quants, timestep_to_full, stat_group = "percentile_95") %>% stack())
#   
#   




)
# End of _targets.R
