# _targets.R
# Libraries
# install.packages("devtools")
# library(devtools)
# devtools::install_version("Matrix", version = "1.2.8")
# library(Matrix)
# devtools::install_version("TMB", "1.7.18")
# devtools::install_github("James-Thorson-NOAA/FishStatsUtils", ref = "development", force = TRUE, upgrade = FALSE)
# devtools::install_github("James-Thorson-NOAA/VAST", ref = "development", force = TRUE, upgrade = FALSE)
# library(VAST)
# library(FishStatsUtils)
# library(tidyverse)
# library(lubridate)
# library(sf)
# library(here)
library(targets)

# Functions
source(here::here("scratch/aja/targets_flow/R/dfo_functions.R"))
source(here::here("scratch/aja/targets_flow/R/nmfs_functions.R"))
source(here::here("scratch/aja/targets_flow/R/combo_functions.R"))
source(here::here("scratch/aja/targets_flow/R/covariate_functions.R"))
#source(here::here("scratch/aja/targets_flow/R/vast_functions.R"))

# Targets set up
options(tidyverse.quiet = TRUE)
tar_option_set(packages = c("Matrix", "TMB", "FishStatsUtils", "VAST", "tidyverse", "lubridate", "sf", "here", "tools"))

# Dynamic files
nmfs_raw_dir <- function() {
  "~/Box/RES_Data/NMFS_trawl/processed_data"
}

dfo_raw_dir <- function() {
  here::here("scratch/aja/targets_flow/data/dfo/raw")
}

species_table_dir <- function() {
  here::here("scratch/aja/targets_flow/data/supporting")
}

static_covariates_dir <- function() {
  here::here("scratch/aja/targets_flow/data/covariates/static")
}

dynamic_covariates_dir <- function() {
  here::here("scratch/aja/targets_flow/data/covariates/dynamic")
}

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
    command = dfo_get_tows(dfo_GSINF = dfo_GSINF, dfo_GSMISSIONS = dfo_GSMISSIONS, out_dir = here::here("scratch/aja/targets_flow/data/dfo/clean"))
  ),
  
  # DFO tidy occurrence data
  tar_target(
    name = dfo_tidy_occu,
    command = dfo_make_tidy_occu(dfo_GSCAT = dfo_GSCAT, dfo_tows = dfo_tows, species_table = species, out_dir = here::here("scratch/aja/targets_flow/data/dfo/clean"))
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
    command = nmfs_get_tows(nmfs_raw = nmfs, out_dir = here::here("scratch/aja/targets_flow/data/nmfs/clean"))
  ),
  
  # NMFS tidy occurrence data
  tar_target(
    name = nmfs_tidy_occu,
    command = nmfs_make_tidy_occu(nmfs_raw = nmfs, nmfs_tows = nmfs_tows, species_table = species, out_dir = here::here("scratch/aja/targets_flow/data/nmfs/clean"))
  ),

  # Combine tow data
  tar_target(
    name = all_tows,
    command = bind_nmfs_dfo_tows(nmfs_tows = nmfs_tows, dfo_tows = dfo_tows, out_dir = here::here("scratch/aja/targets_flow/data/combined"))
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
    command = static_extract_wrapper(static_covariates_list = static_covariates_stack, sf_points = all_tows_sf, date_col_name = "EST_DATE", df_sf = "sf", out_dir = here::here("scratch/aja/targets_flow/data/combined"))
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
    command = dynamic_2d_extract_wrapper(dynamic_covariates_list = dynamic_covariates_stack, sf_points = all_tows_with_static_covs, date_col_name = "DATE", df_sf = "df", out_dir = here::here("scratch/aja/targets_flow/data/combined"))
  ),
  
  # Combine tidy occupancy data
  tar_target(
    name = all_tidy_occu,
    command = bind_nmfs_dfo_tidy_occu(nmfs_tidy_occu = nmfs_tidy_occu, dfo_tidy_occu = dfo_tidy_occu, out_dir = here::here("scratch/aja/targets_flow/data/combined"))
  ),
  
  # Make tidy model data set
  tar_target(
    name = tidy_mod_data,
    command = make_tidy_mod_data(all_tidy_occu = all_tidy_occu, all_tows = all_tows_with_all_covs, out_dir = here::here("scratch/aja/targets_flow/data/combined"))
  )
   
  # # Make VAST covariate data object
  # tar_target(
  #   name = vast_cov_dat,
  #   command = vast_make_cov_data(all_cov_dir = here::here("scratch/aja/targets_flow/data/combined"), year_min = 1985, year_max = 2015, out_dir = here::here("scratch/aja/targets_flow/data/combined"))
  # ),
  # 
  # # Make VAST sample data object
  # tar_target(
  #   name = vast_sample_dat,
  #   command = vast_make_sample_data(all_sample_dir = here::here("scratch/aja/targets_flow/data/combined"), year_min = 1985, year_max = 2015, nmfs_species_code = 73, out_dir = here::here("scratch/aja/targets_flow/data/combined"))
  # ),
  # 
  # # Read in shapefile
  # tar_target(
  #   name = shapefile,
  #   command = read_polyshape("~/Box/RES_Data/Shapefiles/NELME_regions/NELME_sf.shp")
  # ),
  # 
  # # Make extrapolation grid
  # tar_target(
  #   name = vast_extrap_grid,
  #   command = vast_make_extrap_grid(shapefile = shapefile, cell_size = 50000)
  # ),
  # 
  # # Make settings
  # tar_target(
  #   name = vast_settings,
  #   command = vast_make_settings(extrap_grid = vast_extrap_grid, FieldConfig = c("Omega1" = 1, "Epsilon1" = 1, "Omega2" = 1, "Epsilon2" = 1), RhoConfig = c("Beta1" = 3, "Beta2" = 3, "Epsilon1" = 4, "Epsilon2" = 4), bias.correct = FALSE)
  # ),
  # 
  # # Make covariate effect vectors
  # tar_target(
  #   name = vast_coveff,
  #   command = vast_make_coveff(X1_coveff_vec = c(2, 3, 3, 2, rep(3, nlevels(vast_cov_dat$Year_Cov)-1)), X2_coveff_vec = c(2, 3, 3, 2, rep(3, nlevels(vast_cov_dat$Year_Cov)-1)))
  # ),
  # 
  # # Build base model
  # tar_target(
  #   name = vast_build,
  #   command = vast_build_sdm(settings = vast_settings, extrap_grid = vast_extrap_grid, samp_dat = vast_sample_dat, cov_dat = vast_cov_dat, X1_formula = ~ Season + Year_Cov, X2_formula = ~ Season + Year_Cov, X_contrasts = list(Season = contrasts(vast_cov_dat$Season, contrasts = FALSE), Year_Cov = contrasts(vast_cov_dat$Year_Cov, contrasts = FALSE)), Xconfig_list = vast_coveff)
  # ),
  # 
  # # Make adjustments
  # tar_target(
  #   name = vast_adjust,
  #   command = vast_make_adjustments(vast_build, adjustments = list("log_sigmaXi1_cp" = factor(c(rep(1, 3), rep(4, nlevels(vast_cov_dat$Year_Cov)))), "log_sigmaXi2_cp" = factor(c(rep(1, 3), rep(4, nlevels(vast_cov_dat$Year_Cov))))))
  # ),
  # 
  # # Fit model, either the base model OR by making some adjustments
  # tar_target(
  #   name = vast_fit,
  #   command = vast_fit_sdm(vast_adjust)
  # )
)

# End of _targets.R