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
library(raster)
library(here)
library(targets)
library(tarchetypes)
library(patchwork)
library(gifski)
library(akima)
library(splines)
library(parallel)
library(doFuture)
library(tools)

# Functions
manual<- FALSE
if(manual){
  source(here::here("scratch/aja/TargetsSDM/R/dfo_functions.R"))
  source(here::here("scratch/aja/TargetsSDM/R/nmfs_functions.R"))
  source(here::here("scratch/aja/TargetsSDM/R/combo_functions.R"))
  source(here::here("scratch/aja/TargetsSDM/R/enhance_r_funcs.R"))
  source(here::here("scratch/aja/TargetsSDM/R/vast_functions.R"))
} else {
  source(here::here("R/dfo_functions.R"))
  source(here::here("R/nmfs_functions.R"))
  source(here::here("R/combo_functions.R"))
  source(here::here("R/enhance_r_funcs.R"))
  source(here::here("R/vast_functions.R"))
}


# Targets set up
options(tidyverse.quiet = TRUE)
tar_option_set(packages = c("Matrix", "TMB", "FishStatsUtils", "VAST", "tidyverse", "lubridate", "sf", "raster", "here", "tools"))

# Model fitting stuffs
nmfs_species_code<- 301
nice_category_names<- "American lobster"
depth_cut<- 400
fit_year_min<- 1985
fit_year_max<- 2015
fit_seasons<- c("SPRING", "SUMMER", "FALL")
#fit_seasons<- c("FALL")
pred_years <- 2100
gam_degree<- 3
hab_formula<- ~ Season + Year_Cov + bs(Depth, degree = 3, intercept = FALSE) + bs(SST_seasonal, degree = 3, intercept = FALSE) + bs(BT_seasonal, degree = 3, intercept = FALSE) # Seasonal
#hab_formula<- ~ bs(Depth, degree = 2, intercept = FALSE) + bs(SST_seasonal, degree = 2, intercept = FALSE) # Annnual
hab_env_coeffs_n<- length(attributes(terms.formula(hab_formula))$term.labels)-2 # Seasonal
#hab_env_coeffs_n<- length(attributes(terms.formula(hab_formula))$term.labels)-1 # Annual
# hab_env_coeffs_n<- length(attributes(terms.formula(hab_formula))$term.labels) # Annual habitat covs only
field_config<- c("Omega1" = 1, "Epsilon1" = 1, "Omega2" = 1, "Epsilon2" = 1)
rho_config<- c("Beta1" = 2, "Beta2" = 2, "Epsilon1" = 2, "Epsilon2" = 2)

catch_formula<- ~ factor(Survey)
#catch_formula<- ~0
strata_use<- data.frame("STRATA" = c("NMFS_and_DFO", "DFO", "NMFS", "Scotian_Shelf", "Georges_Bank", "Gulf_of_Maine", "Southern_New_England", "Mid_Atlantic_Bight"))
#strata_use<- data.frame("STRATA" = c("NMFS", "DFO"))
#strata_use<- data.frame("STRATA" = c("NMFS"))


# Dynamic files
nmfs_raw_mk_dir <- function() {
  "~/Box/RES_Data/NMFS_trawl/processed_data"
}

dfo_raw_mk_dir <- function() {
  here::here("data/dfo/raw")
}

species_table_mk_dir <- function() {
  here::here("data/supporting")
}

static_covariates_mk_dir <- function() {
  here::here("data/covariates/static")
}

dynamic_covariates_mk_dir <- function() {
  here::here("data/covariates/dynamic")
}

predict_covariates_raw_mk_dir <- function() {
  "~/Box/RES_Data/CMIP6/BiasCorrected"
}

# predict_covariates_processed_dir <- function() {
#   here::here("data/predict")
# }

predict_template_mk_dir <- function() {
  here::here("data/supporting")
}

land_sf_mk_dir <- function() {
  "~/Box/RES_Data/Shapefiles/ne_50m_land/ne_50m_land.shp"
}

index_shapefiles_mk_dir <- function() {
  here::here("data/supporting/index_shapefiles")
}

region_shapefile_mk_dir <- function() {
  here::here("data/supporting/region_shapefile")
}

# high_res_dir <- function() {
#   here::here("data/supporting")
# }

# Targets list of commands
list(
  # Species table directory
  tar_target(
    name = species_table,
    command = species_table_mk_dir(),
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
    command =  land_sf_mk_dir(),
  ),

  # Land sf file
  tar_target(
    name = land_sf,
    command =  land_read_sf(land_sf_dir),
  ),

  # Get DFO directory
  tar_target(
    name = dfo_raw,
    command = dfo_raw_mk_dir(),
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
    command = nmfs_raw_mk_dir(),
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
    command = static_covariates_mk_dir(),
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
    command = dynamic_covariates_mk_dir(),
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
  
  tar_target(
    name = all_tows_with_all_covs_rescale,
    command = rescale_all_covs(all_tows_with_all_covs, depth_cut = depth_cut, type = "AJA", center = TRUE, scale = TRUE, out_dir = here::here("data/combined"))
  ),
  
  tar_target(
    name = rescale_params,
    command = get_rescale_params(all_tows_with_all_covs, depth_cut = depth_cut, out_dir = here::here("data/covariates/"))
  ),

  # Combine tidy occupancy data
  tar_target(
    name = all_tidy_occu,
    command = bind_nmfs_dfo_tidy_occu(nmfs_tidy_occu = nmfs_tidy_occu, dfo_tidy_occu = dfo_tidy_occu, out_dir = here::here("data/combined"))
  ),

  # Make tidy model data set
  tar_target(
    name = tidy_mod_data,
    command = make_tidy_mod_data(all_tidy_occu = all_tidy_occu, all_tows = all_tows_with_all_covs_rescale, out_dir = here::here("data/combined"))
  ),

  # Things get a bit weird here. Before we go into making the VAST specific data, we are going to want to gather up information for the locations and times we want to predict. These observations will be dummy observations and they won't be used in estimating the likelihood of the model, but the model will predict to them. There has got to be a better way of doing this -- but for now, just to make sure we have what we need and try to make things as flexible as possible for future use.

  # First step, need to gather up our prediction covariates
  # Prediction dataset directory
  tar_target(
    name = predict_covariates_raw_dir,
    command = predict_covariates_raw_mk_dir(),
    format = "file",
  ),

  # Read in raw covariates and summarize them
  tar_target(
    name = predict_covariates_stack_agg_out,
    command =  predict_covariates_stack_agg(predict_covariates_raw_dir, ensemble_stat = "mean", summarize = "seasonal", resample_template = here::here("data/supporting", "Rast0.25grid.grd"), out_dir = here::here("data/predict")),
    format = "file",
  ),

  # # Processed covariate directory
  # tar_target(
  #   name = predict_covariates_processed_dir,
  #   command = predict_covariates_processed_dir(),
  #   format = "file",
  # ),

  # Process prediction covariates into a dataframe that can be eventually joined up to vast_seasonal_data, masked to focal region of interest
  tar_target(
    name = region_shapefile_dir,
    command = region_shapefile_mk_dir(),
    format = "file"
  ),

  # Read in region shapefile
  tar_target(
    name = region_shapefile,
    command = vast_read_region_shape(region_shapefile_dir)
  ),

  # Read in index shapefiles
  tar_target(
    name = index_shapefiles_dir,
    command = index_shapefiles_mk_dir(),
    format = "file"
  ),

  tar_target(
    name = index_shapefiles,
    command = vast_read_index_shapes(index_shapefiles_dir)
  ),

  tar_target(
    name = vast_predict_df,
    command =  make_vast_predict_df(predict_covariates_stack_agg = predict_covariates_stack_agg_out, extra_covariates_stack = static_covariates_stack, covs_rescale = c("Depth", "BS_seasonal", "BT_seasonal", "SS_seasonal", "SST_seasonal"), rescale_params = rescale_params, depth_cut = depth_cut, mask = region_shapefile, summarize = "seasonal", ensemble_stat = "mean", fit_year_min = fit_year_min, fit_year_max = fit_year_max, fit_seasons = fit_seasons, pred_years = pred_years, out_dir = here::here("data/predict")),
  ),
  
  # Make extrapolation grid
  tar_target(
    name = vast_extrap_grid,
    command = vast_make_extrap_grid(region_shapefile = region_shapefile, index_shapes = index_shapefiles, strata.limits = strata_use, cell_size = 25000)
  ),
  
  # Make settings 
  tar_target(
    name = vast_settings,
    command = vast_make_settings(extrap_grid = vast_extrap_grid, n_knots = 400, FieldConfig = field_config, RhoConfig = rho_config, OverdispersionConfig = c(0, 0), bias.correct = TRUE, knot_method = "samples", inla_method = "Barrier", Options = c("Calculate_Range"=TRUE), strata.limits = strata_use)
  ),
  
  # Make spatial lists 
  tar_target(
    name = vast_spatial_lists,
    command = vast_make_spatial_lists(extrap_grid = vast_extrap_grid, vast_settings = vast_settings, tidy_mod_data = tidy_mod_data, out_dir = paste0(here::here(""), "/"))
  ),
  
  # Reduce the prediction dataframe from regular grid to have size based on number of knots
  tar_target(
    name = vast_predict_df_reduced,
    command = reduce_vast_predict_df(vast_predict_df = vast_predict_df, vast_spatial_lists = vast_spatial_lists, out_dir = here::here("data/predict"))
  ), 
  
  # Make VAST seasonal data
  tar_target(
    name = vast_seasonal_data,
    command = make_vast_seasonal_data(tidy_mod_data = tidy_mod_data, fit_seasons = fit_seasons, nmfs_species_code = nmfs_species_code, fit_year_min = fit_year_min, fit_year_max = fit_year_max, pred_years = pred_years, pred_df = vast_predict_df_reduced, out_dir = here::here("data/combined"))
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

  # Make VAST catachability data object
  tar_target(
    name = vast_catchability_data,
    command = make_vast_catchability_data(vast_seasonal_data = vast_seasonal_data, out_dir = here::here("data/combined"))
  ),

  # # Make covariate effect vectors -- just habitat covariates
  # tar_target(
  #   name = vast_coveff,
  #   command = vast_make_coveff(X1_coveff_vec = rep(rep(1, gam_degree), hab_env_coeffs_n), X2_coveff_vec = rep(rep(1, gam_degree), hab_env_coeffs_n), Q1_coveff_vec = NULL, Q2_coveff_vec = NULL)
  # ),
  
  # # Make covariate effect vectors -- annual intercept
  # tar_target(
  #   name = vast_coveff,
  #   command = vast_make_coveff(X1_coveff_vec = c(2, rep(3, nlevels(vast_covariate_data$Year_Cov)-1), rep(rep(1, gam_degree), hab_env_coeffs_n)), X2_coveff_vec = c(2, rep(3, nlevels(vast_covariate_data$Year_Cov)-1), rep(rep(1, gam_degree), hab_env_coeffs_n)), Q1_coveff_vec = NULL, Q2_coveff_vec = NULL)
  # ),

  # # Make covariate effect vectors -- multiple seasons
  tar_target(
    name = vast_coveff,
    command = vast_make_coveff(X1_coveff_vec = c(2, rep(3, length(unique(fit_seasons))-1), 2, rep(3, nlevels(vast_covariate_data$Year_Cov)-1), rep(rep(1, gam_degree), hab_env_coeffs_n)), X2_coveff_vec = c(2, rep(3, length(unique(fit_seasons))-1), 2, rep(3, nlevels(vast_covariate_data$Year_Cov)-1), rep(rep(1, gam_degree), hab_env_coeffs_n)), Q1_coveff_vec = NULL, Q2_coveff_vec = NULL)
  ),
  # # 
  # # # # Build base model -- annual
  # tar_target(
  #   name = vast_build0,
  #   command = vast_build_sdm(settings = vast_settings, extrap_grid = vast_extrap_grid, Method = "Barrier", sample_data = vast_sample_data, covariate_data = vast_covariate_data, X1_formula = hab_formula, X2_formula = hab_formula, X_contrasts = list(Year_Cov = contrasts(vast_covariate_data$Year_Cov, contrasts = FALSE)), Xconfig_list = vast_coveff, catchability_data = vast_catchability_data, Q1_formula = catch_formula, Q2_formula = catch_formula, index_shapes = index_shapefiles)
  # ),
  #
  # # Build base model -- seasonal
  tar_target(
    name = vast_build0,
    command = vast_build_sdm(settings = vast_settings, extrap_grid = vast_extrap_grid, sample_data = vast_sample_data, covariate_data = vast_covariate_data, X1_formula = hab_formula, X2_formula = hab_formula, X_contrasts = list(Season = contrasts(vast_covariate_data$Season, contrasts = FALSE), Year_Cov = contrasts(vast_covariate_data$Year_Cov, contrasts = FALSE)), Xconfig_list = vast_coveff, catchability_data = vast_catchability_data, Q1_formula = catch_formula, Q2_formula = catch_formula, index_shapes = index_shapefiles)
  ),

  # # # Make adjustments -- annual just habitat covariates, no adjustments needed
  # tar_target(
  #   name = vast_adjust,
  #   command = vast_make_adjustments(vast_build = vast_build0, adjustments = NULL, index_shapes = index_shapefiles)
  # ),

  # # # Make adjustments -- annual
  # tar_target(
  #   name = vast_adjust,
  #   command = vast_make_adjustments(vast_build = vast_build0, adjustments = list("log_sigmaXi1_cp" = factor(c(rep(4, nlevels(vast_covariate_data$Year_Cov)), rep(NA, gam_degree*hab_env_coeffs_n))), "log_sigmaXi2_cp" = factor(c(rep(4, nlevels(vast_covariate_data$Year_Cov)), rep(NA, gam_degree*hab_env_coeffs_n)))), index_shapes = index_shapefiles)
  # ),

  # Make adjustments -- seasonal
  tar_target(
    name = vast_adjust,
    command = vast_make_adjustments(vast_build = vast_build0, adjustments = list("log_sigmaXi1_cp" = factor(c(rep(1, length(unique(fit_seasons))), rep(4, nlevels(vast_covariate_data$Year_Cov)), rep(NA, gam_degree*hab_env_coeffs_n))), "log_sigmaXi2_cp" = factor(c(rep(1, length(unique(fit_seasons))), rep(4, nlevels(vast_covariate_data$Year_Cov)), rep(NA, gam_degree*hab_env_coeffs_n))), "lambda1_k" = factor(c(1, NA)), "lambda2_k" = factor(c(1, NA))), index_shapes = index_shapefiles)
  ),

  # Fit model, either the base model OR by making some adjustments
  tar_target(
    name = vast_fit,
    command = vast_fit_sdm(vast_build_adjust = vast_adjust, nmfs_species_code = nmfs_species_code, out_dir = here::here("results/mod_fits"), index_shapes = index_shapefiles)
  ),

  # Plot samples and knot locations
  tar_target(
    name = vast_sample_knot_plot,
    command = vast_plot_design(vast_fit = vast_fit, land = land_sf, spat_grid = here::here("data/predict/predict_stack_SST_seasonal_mean.grd"), xlim = c(-80, -55), ylim = c(35, 50), land_color = "#d9d9d9", out_dir = here::here("results/plots_maps/"))
  ),
  # 
  # # Get covariate effects --
  tar_target(
    name = vast_covariate_effects,
    command = get_vast_covariate_effects(vast_fit = vast_fit, params_plot = c("Depth", "SST_seasonal", "BT_seasonal"), params_plot_levels = 100, effects_pad_values = c(1), nice_category_names = nice_category_names, out_dir = here::here("results/tables/"))
  ),
  # 
  # tar_target(
  #   name = vast_covariate_effects,
  #   command = get_vast_covariate_effects(vast_fit = vast_fit, params_plot = c("Depth", "SST_seasonal"), params_plot_levels = 100, effects_pad_values = c(), nice_category_names = nice_category_names, out_dir = here::here("results/tables/"))
  # ),
  
  # tar_target(
  #   name = vast_covariate_effects,
  #   command = get_vast_covariate_effects(vast_fit = vast_fit, params_plot = c("Depth", "SST_seasonal"), params_plot_levels = 100, effects_pad_values = c(1), nice_category_names = nice_category_names, out_dir = here::here("results/tables/"))
  # ),
  
  # Plot covariate effects
  tar_target(
    name = vast_covariate_effects_plot,
    command = plot_vast_covariate_effects(vast_covariate_effects = vast_covariate_effects, vast_fit = vast_fit, nice_category_names = nice_category_names, out_dir = here::here("results/plots_maps/"))
  ),

  # Predict with fitted model
  tar_target(
    name = vast_predictions,
    command = predict_vast(vast_fitted_sdm = vast_fit, nmfs_species_code = nmfs_species_code, predict_variable = "D_i", predict_category = 0, predict_vessel = 0, predict_covariates_df_all = vast_predict_df, out_dir = here::here("results/prediction_df"))
  ),

  # Create gif from predictions
  tar_target(
    name = plot_preds,
    command = vast_fit_plot_density(vast_fit = vast_fit, nice_category_names = nice_category_names, mask = region_shapefile, all_times = as.character(levels(vast_seasonal_data$VAST_YEAR_SEASON)), plot_times = NULL, land_sf = land_sf, xlim = c(-80, -55), ylim = c(35, 50), panel_or_gif = "gif", out_dir = here::here("results/plots_maps"))
  ),

  # tar_target(
  #   name = plot_preds,
  #   comman = vast_df_plot_density(pred_df = vast_predictions, nice_category_names = nice_category_names, mask = region_shapefile, all_times = as.character(levels(vast_seasonal_data$VAST_YEAR_SEASON)), plot_times = NULL, land_sf = land_sf, xlim = c(-80, -55), ylim = c(35, 50), panel_or_gif = "gif", out_dir = here::here("results/plots_maps"))
  # ),
  # 
  # # Calculate biomass index values
  tar_target(
    name = biomass_indices,
    command = get_vast_index_timeseries(vast_fit = vast_fit, nice_category_names = nice_category_names, index_scale = c("raw"), out_dir = here::here("results/tables"))
  ),
  # 
  # # Plot biomass index time series
  tar_target(
    name = plot_biomass_timeseries,
    command = plot_vast_index_timeseries(index_res_df = biomass_indices, index_scale = "raw", nice_category_names = nice_category_names, nice_xlab = "Year-Season", nice_ylab= "Biomass index (metric tons)", paneling = "none", color_pal = NULL, out_dir = here::here("results/plots_maps"))
  )
)

# End of _targets.R