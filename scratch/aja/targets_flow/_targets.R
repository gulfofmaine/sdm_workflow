# _targets.R
# Libraries
install.packages("devtools")
library(devtools)
devtools::install_version("Matrix", version = "1.2.8")
library(Matrix)
devtools::install_version("TMB", "1.7.18")
devtools::install_github("James-Thorson-NOAA/FishStatsUtils", ref = "development", force = TRUE, upgrade = FALSE)
devtools::install_github("James-Thorson-NOAA/VAST", ref = "development", force = TRUE, upgrade = FALSE)
library(VAST)
library(FishStatsUtils)
library(tidyverse)
library(sf)
library(here)
library(targets)

# Functions
source(here::here("scratch/aja/targets_flow/R/vast_functions.R"))

# Targets set up
options(tidyverse.quiet = TRUE)
tar_option_set(packages = c("Matrix", "TMB", "FishStatsUtils", "VAST", "tidyverse", "sf", "here"))

# Targets list of commands
list(
  # Read in sample data
  tar_target(
    name = samp_dat,
    command = read_samp_dat_csv(here::here("scratch/aja/targets_flow/data", "sample_data.csv"))
  ),
  
  # Read in covariate data
  tar_target(
    name = cov_dat,
    command = read_cov_dat_csv(here::here("scratch/aja/targets_flow/data", "covariate_data.csv"))
  ),
  
  # Read in shapefile
  tar_target(
    name = shapefile,
    command = read_polyshape("~/Box/RES_Data/Shapefiles/NELME_regions/NELME_sf.shp")
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
    command = vast_make_coveff(X1_coveff_vec = c(2, 3, 3, 2, rep(3, 32)), X2_coveff_vec = c(2, 3, 3, 2, rep(3, 32)))
  ),
  
  # Build base model
  tar_target(
    name = vast_build,
    command = vast_build_sdm(settings = vast_settings, extrap_grid = vast_extrap_grid, samp_dat = samp_dat, cov_dat = cov_dat, X1_formula = ~ Season + Year_Cov, X2_formula = ~ Season + Year_Cov, X_contrasts = list(Season = contrasts(cov_dat$Season, contrasts = FALSE), Year_Cov = contrasts(cov_dat$Year_Cov, contrasts = FALSE)), Xconfig_list = vast_coveff)
  ),
  
  # Make adjustments
  tar_target(
    name = vast_adjust,
    command = vast_make_adjustments(vast_build, adjustments = list("log_sigmaXi1_cp" = factor(c(rep(1, 3), rep(4, 33))), "log_sigmaXi2_cp" = factor(c(rep(1, 3), rep(4, 33)))))
  ),
  
  # Fit model, either the base model OR by making some adjustments
  tar_target(
    name = vast_fit,
    command = vast_fit_sdm(vast_adjust)
  )
)

# End of _targets.R