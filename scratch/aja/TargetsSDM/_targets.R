# _targets.R
# Libraries
# install.packages("devtools")
library(devtools)
# devtools::install_version("Matrix", version = "1.2.8")
library(Matrix)
# devtools::install_version("TMB", "1.7.18")
# devtools::install_github("James-Thorson-NOAA/FishStatsUtils", force = TRUE, upgrade = FALSE)
# devtools::install_github("James-Thorson-NOAA/VAST", ref = "main", force = TRUE, upgrade = FALSE)
library(VAST)
library(FishStatsUtils)
library(tidyverse)
library(lubridate)
library(sf)
library(here)
library(targets)
library(tarchetypes)
library(gifski)
library(akima)
library(splines)
library(parallel)
library(doFuture)

# Functions
source(here::here("R/dfo_functions.R"))
source(here::here("R/nmfs_functions.R"))
source(here::here("R/combo_functions.R"))
source(here::here("R/covariate_functions.R"))
source(here::here("R/vast_functions.R"))

# Targets set up
options(tidyverse.quiet = TRUE)
tar_option_set(packages = c("Matrix", "TMB", "FishStatsUtils", "VAST", "tidyverse", "lubridate", "sf", "raster", "here", "tools"))

# Model fitting stuffs
nmfs_species_code<- 105
fit_year_min<- 1985
fit_year_max<- 2017
fit_seasons<- c("DFO", "SPRING", "FALL")
pred_years <- 2025
mod_formula<- ~ Season + Year_Cov + bs(SST_seasonal, degree = 2, intercept = FALSE)

# Dynamic files
nmfs_raw_dir <- function() {
  "~/Box/RES_Data/NMFS_trawl/processed_data"
}

dfo_raw_dir <- function() {
  here::here("data/dfo/raw")
}

species_table_dir <- function() {
  here::here("data/supporting")
}

static_covariates_dir <- function() {
  here::here("data/covariates/static")
}

dynamic_covariates_dir <- function() {
  here::here("data/covariates/dynamic")
}

predict_covariates_raw_dir <- function() {
  "~/Box/RES_Data/CMIP6/BiasCorrected"
}

predict_covariates_processed_dir <- function() {
  here::here("data/predict")
}

predict_template_dir <- function() {
  here::here("data/supporting")
}

land_sf_dir <- function() {
  "~/Box/RES_Data/Shapefiles/ne_50m_land/ne_50m_land.shp"
}

# high_res_dir <- function() {
#   here::here("data/supporting")
# }

# Targets list of commands
list(
  # Species table directory
  tar_target(
    name = species_table,
    command = species_table_dir(),
    format = "file",
  ),
  
  # Species table load
  tar_target(
    name = species,
    command =  species_read_csv(species_table),
  ),
  
  # Land sf file
  tar_target(
    name = land_sf_dir,
    command =  land_sf_dir(),
  ),
  
  # Land sf file
  tar_target(
    name = land_sf,
    command =  land_read_sf(land_sf_dir),
  ),
  
  # Get DFO directory
  tar_target(
    name = dfo_raw,
    command = dfo_raw_dir(),
    format = "file",
  ),
  
  # DFO load GSINF
  tar_target(
    name = dfo_GSINF,
    command = dfo_GSINF_load(dfo_raw),
  ),

  # DFO load GSMISSIONS
  tar_target(
    name = dfo_GSMISSIONS,
    command = dfo_GSMISSIONS_load(dfo_raw),
  ),
  
  # DFO load GSCAT
  tar_target(
    name = dfo_GSCAT,
    command = dfo_GSCAT_load(dfo_raw),
  ),
  
  # DFO tows
  tar_target(
    name = dfo_tows,
    command = dfo_get_tows(dfo_GSINF = dfo_GSINF, dfo_GSMISSIONS = dfo_GSMISSIONS, out_dir = here::here("data/dfo/clean"))
  ),
  
  # DFO tidy occurrence data
  tar_target(
    name = dfo_tidy_occu,
    command = dfo_make_tidy_occu(dfo_GSCAT = dfo_GSCAT, dfo_tows = dfo_tows, species_table = species, out_dir = here::here("data/dfo/clean"))
  ),
  
  # NMFS directory
  tar_target(
    name = nmfs_raw,
    command = nmfs_raw_dir(),
    format = "file"
  ),
  
  # NMFS load
  tar_target(
    name = nmfs,
    command = nmfs_load(nmfs_raw),
  ),
  
  # NMFS tows
  tar_target(
    name = nmfs_tows,
    command = nmfs_get_tows(nmfs_raw = nmfs, out_dir = here::here("data/nmfs/clean"))
  ),
  
  # NMFS tidy occurrence data
  tar_target(
    name = nmfs_tidy_occu,
    command = nmfs_make_tidy_occu(nmfs_raw = nmfs, nmfs_tows = nmfs_tows, species_table = species, out_dir = here::here("data/nmfs/clean"))
  ),

  # Combine tow data
  tar_target(
    name = all_tows,
    command = bind_nmfs_dfo_tows(nmfs_tows = nmfs_tows, dfo_tows = dfo_tows, out_dir = here::here("data/combined"))
  ),
  
  # Get additional covariates
  # Create tows sf object
  tar_target(
    name = all_tows_sf,
    command = points_to_sf(all_tows)
  ),
  
  # Static first...
  tar_target(
    name = static_covariates,
    command = static_covariates_dir(),
    format = "file",
  ),
  
  tar_target(
    name = static_covariates_stack,
    command = static_covariates_read(static_covariates)
  ),
  
  tar_target(
    name = all_tows_with_static_covs,
    command = static_extract_wrapper(static_covariates_list = static_covariates_stack, sf_points = all_tows_sf, date_col_name = "DATE", df_sf = "sf", out_dir = here::here("data/combined"))
  ),
  
  # Now dynamic ones...
  tar_target(
    name = dynamic_covariates,
    command = dynamic_covariates_dir(),
    format = "file",
  ),
  
  tar_target(
    name = dynamic_covariates_stack,
    command = dynamic_covariates_read(dynamic_covariates)
  ),
  
  tar_target(
    name = all_tows_with_all_covs,
    command = dynamic_2d_extract_wrapper(dynamic_covariates_list = dynamic_covariates_stack, t_summ = "seasonal", t_position = NULL, sf_points = all_tows_with_static_covs, date_col_name = "DATE", df_sf = "df", out_dir = here::here("data/combined"))
  ),

  # Combine tidy occupancy data
  tar_target(
    name = all_tidy_occu,
    command = bind_nmfs_dfo_tidy_occu(nmfs_tidy_occu = nmfs_tidy_occu, dfo_tidy_occu = dfo_tidy_occu, out_dir = here::here("data/combined"))
  ),

  # Make tidy model data set
  tar_target(
    name = tidy_mod_data,
    command = make_tidy_mod_data(all_tidy_occu = all_tidy_occu, all_tows = all_tows_with_all_covs, out_dir = here::here("data/combined"))
  ),
  
  # Things get a bit weird here. Before we go into making the VAST specific data, we are going to want to gather up information for the locations and times we want to predict. These observations will be dummy observations and they won't be used in estimating the likelihood of the model, but the model will predict to them. There has got to be a better way of doing this -- but for now, just to make sure we have what we need and try to make things as flexible as possible for future use. 
  
  # First step, need to gather up our prediction covariates
  # Prediction dataset directory
  tar_target(
    name = predict_covariates_raw_dir,
    command = predict_covariates_raw_dir(),
    format = "file",
  ),

  # Read in raw covariates and summarize them
  tar_target(
    name = predict_covariates_stack_agg,
    command =  predict_covariates_stack_agg(predict_covariates_raw_dir, ensemble_stat = "mean", summarize = "seasonal", resample_template = here::here("data/supporting", "Rast0.25grid.grd"), out_dir = here::here("data/predict/")),
  ),
  
  # Processed covariate directory
  tar_target(
    name = predict_covariates_processed_dir,
    command = predict_covariates_processed_dir(),
    format = "file",
  ),

  # Process prediction covariates into a dataframe that can be eventually joined up to vast_seasonal_data, masked to focal region of interest
  # Read in shapefile
  tar_target(
    name = shapefile,
    command = read_polyshape("~/Box/RES_Data/Shapefiles/NELME_regions/NMFSandDFO_sf.shp")
  ),
  
  tar_target(
    name = vast_predict_df,
    command =  make_vast_predict_df(predict_covariates_processed_dir, extra_covariates_stack = static_covariates_stack, mask = shapefile, summarize = "seasonal", ensemble_stat = "mean", fit_year_min = fit_year_min, fit_year_max = fit_year_max, fit_seasons = fit_seasons, pred_years = pred_years, out_dir = here::here("data/predict")),
  ),

  # Make VAST seasonal data
  tar_target(
    name = vast_seasonal_data,
    command = make_vast_seasonal_data(tidy_mod_data = tidy_mod_data, fit_seasons = fit_seasons, nmfs_species_code = nmfs_species_code, fit_year_min = fit_year_min, fit_year_max = fit_year_max, pred_df = vast_predict_df, out_dir = here::here("data/combined"))
  ),

  # Make VAST sample data object
  tar_target(
    name = vast_sample_data,
    command = make_vast_sample_data(vast_seasonal_data = vast_seasonal_data, out_dir = here::here("data/combined"))
  ),

  # Make VAST covariate data object
  tar_target(
    name = vast_covariate_data,
    command = make_vast_covariate_data(vast_seasonal_data = vast_seasonal_data, out_dir = here::here("data/combined"))
  ),

  # Make extrapolation grid
  tar_target(
    name = vast_extrap_grid,
    command = vast_make_extrap_grid(shapefile = shapefile, cell_size = 50000)
  ),

  # Make settings
  tar_target(
    name = vast_settings,
    command = vast_make_settings(extrap_grid = vast_extrap_grid, FieldConfig = c("Omega1" = 1, "Epsilon1" = 1, "Omega2" = 1, "Epsilon2" = 1), RhoConfig = c("Beta1" = 3, "Beta2" = 3, "Epsilon1" = 4, "Epsilon2" = 4), bias.correct = FALSE)
  ),

  # Make covariate effect vectors
  tar_target(
    name = vast_coveff,
    command = vast_make_coveff(X1_coveff_vec = c(2, 3, 3, 2, rep(3, nlevels(vast_covariate_data$Year_Cov)-1), rep(1, 2)), X2_coveff_vec = c(2, 3, 3, 2, rep(3, nlevels(vast_covariate_data$Year_Cov)-1), rep(1, 2)))
  ),

  # Build base model
  tar_target(
    name = vast_build0,
    command = vast_build_sdm(settings = vast_settings, extrap_grid = vast_extrap_grid, sample_data = vast_sample_data, covariate_data = vast_covariate_data, X1_formula = mod_formula, X2_formula = mod_formula, X_contrasts = list(Season = contrasts(vast_covariate_data$Season, contrasts = FALSE), Year_Cov = contrasts(vast_covariate_data$Year_Cov, contrasts = FALSE)), Xconfig_list = vast_coveff)
  ),

  # Make adjustments
  tar_target(
    name = vast_adjust,
    command = vast_make_adjustments(vast_build = vast_build0, adjustments = list("log_sigmaXi1_cp" = factor(c(rep(1, 3), rep(4, nlevels(vast_covariate_data$Year_Cov)), rep(NA, 2))), "log_sigmaXi2_cp" = factor(c(rep(1, 3), rep(4, nlevels(vast_covariate_data$Year_Cov)), rep(NA, 2)))))
  ),

  # Fit model, either the base model OR by making some adjustments
  tar_target(
    name = vast_fit,
    command = vast_fit_sdm(vast_build_adjust = vast_adjust)
  ),

  # # High res template, needed for "pretty" plotting
  # tar_target(
  #   name = high_res_dir,
  #   command = high_res_dir(),
  #   format = "file",
  # ),

  # Create gif from predictions
  tar_target(
    name = plot_preds,
    command = vast_plot_density(vast_fit = vast_fit, mask = shapefile, all_times = as.character(levels(vast_seasonal_data$VAST_YEAR_SEASON)), plot_times = NULL, land_sf = land_sf, xlim = c(-80, -55), ylim = c(35, 50), panel_or_gif = "gif", out_dir = here::here("results/"))
  )
  # # Make predictions from fitted VAST model
  # # Prediction dataset directory
  # tar_target(
  #   name = predict_covariates_data,
  #   command = predict_covariates_dir(),
  #   format = "file",
  # ),
  # 
  # # Read in new covariates and aggregate them
  # tar_target(
  #   name = predict_covariates_stack_agg,
  #   command =  predict_covariates_stack_agg(predict_covariates_data, ensemble_stat = "mean", summarize = "seasonal", resample_template = here::here("data/supporting", "Rast0.25grid.grd"), out_dir = here::here("data/predict/")),
  # )
  # 
  # # Run prediction code
  # tar_target(
  #   name = vast_sdm_predictions,
  #   command =  vast_predict_sdm(x = vast_fit, new_covariate_data = predict_covariates, covariates_use, predict_variable = "D_i", out_dir = here::here("results")),
  # )
)

# End of _targets.R