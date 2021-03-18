##########
##### Walking through steps of _targets.R
##########
#install.packages("usethis")
#install.packages("devtools")
library(devtools)
#install_version("Matrix", version = "1.2.8")
library(Matrix)
#install_version("TMB", "1.7.18")
# Now VAST and FishStatsUtils, do NOT update TMB
#devtools::install_github("James-Thorson-NOAA/FishStatsUtils", ref = "development", force = TRUE)
#devtools::install_github("James-Thorson-NOAA/VAST", ref = "development", force = TRUE)
library(VAST)
library(FishStatsUtils)
library(profvis)
library(future)
library(sf)
library(tidyverse)
library(gmRi)

source(here::here("scratch/aja/targets_flow/R/vast_functions.R"))

shape_path<- shared.path(os.use = "unix", group = "root", folder = "RES_Data/Shapefiles/")

# Data prep
samp_dat<- read_samp_dat_csv(here::here("scratch/aja/targets_flow/data", "sample_data.csv"))
cov_dat<- read_cov_dat_csv(here::here("scratch/aja/targets_flow/data", "covariate_data.csv"))
shapefile<- read_polyshape(paste0(shape_path, "NELME_regions/NELME_sf.shp"))

# Input grid and settings
vast_extrap_grid<- vast_make_extrap_grid(shapefile = shapefile, cell_size = 50000)
vast_settings<- vast_make_settings(extrap_grid = vast_extrap_grid, FieldConfig = c("Omega1" = 1, "Epsilon1" = 1, "Omega2" = 1, "Epsilon2" = 1), RhoConfig = c("Beta1" = 3, "Beta2" = 3, "Epsilon1" = 4, "Epsilon2" = 4), bias.correct = FALSE)
vast_coveff<- vast_make_coveff(X1_coveff_vec = c(2, 3, 3, 2, rep(3, 32)), X2_coveff_vec = c(2, 3, 3, 2, rep(3, 32)))

# Initial build
vast_build<- vast_build_sdm(settings = vast_settings, extrap_grid = vast_extrap_grid, samp_dat = samp_dat, cov_dat = cov_dat, X1_formula = ~ Season + Year_Cov, X2_formula = ~ Season + Year_Cov, X_contrasts = list(Season = contrasts(cov_dat$Season, contrasts = FALSE), Year_Cov = contrasts(cov_dat$Year_Cov, contrasts = FALSE)), Xconfig_list = vast_coveff)

# Adjustments
vast_build_adjust<- vast_make_adjustments(vast_build, adjustments = list("log_sigmaXi1_cp" = factor(c(rep(1, 3), rep(4, 33))), "log_sigmaXi2_cp" = factor(c(rep(1, 3), rep(4, 33)))))

# Final fit
vast_fit<- vast_fit_sdm(vast_build_adjust)
