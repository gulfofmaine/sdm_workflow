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
  tar_target(
    name = cmip_cropped,
    command = import_cmip_sst_test(cmip_file = "tester") ),
  tar_target(
    name = oisst_clim_92,
    command = import_oisst_clim(climatology_period = "1992-2011", os.use = "unix") ),
  tar_target(
    name = oisst_month_avgs,
    command = months_from_clim(clim_source = oisst_clim_92, month_layer_key = NULL) ),
  tar_target(
    name = cmip_clim,
    command = cmip_to_clim(cmip_stack = cmip_cropped) ),
  tar_target(
    name = cmip_anoms,
    command = cmip_get_anomalies(cmip_data = cmip_cropped, cmip_clim = cmip_clim) ),
  tar_target(
    name = cmip_anoms_regridded,
    command = resample_grid(starting_grid = cmip_anoms,  desired_grid = oisst_month_avgs, method = "bilinear") ),
  tar_target(
    name = cmip_delta_bias_corrected,
    command = delta_method_bias_correct(cmip_grid = cmip_anoms_regridded, reference_climatology = oisst_month_avgs) ),
  tar_target(
    name = cmip_unbiased_ensemble,
    command = restack_cmip_projections(cmip_inputs = cmip_delta_bias_corrected) ),
  tar_target(
    name = cmip_sst_quants,
    command = timestep_stats(year_stacks = cmip_unbiased_ensemble) ),
  tar_target(
    name = cmip_sst_mean,
    command = imap(cmip_sst_quants, timestep_to_full, stat_group = "mean") %>% stack() ),
  tar_target(
    name = cmip_sst_05,
    command = imap(cmip_sst_quants, timestep_to_full, stat_group = "percentile_05") %>% stack() ),
  tar_target(
    name = cmip_sst_95,
    command = imap(cmip_sst_quants, timestep_to_full, stat_group = "percentile_95") %>% stack()),
  tar_target(
    name = soda_clim_vars,
    command = map(c("surf_sal", "surf_temp", "bot_sal", "bot_temp"), import_soda_clim, os.use = "unix")
  )
  
  
)



# End of _targets.R
