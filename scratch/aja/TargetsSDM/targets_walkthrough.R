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

# Function files
func_files<- list.files(here::here("scratch/aja/targets_flow/R"), full.names = TRUE)
invisible(lapply(func_files, source))

# Shared paths
shape_path<- shared.path(os.use = "unix", group = "root", folder = "RES_Data/Shapefiles/")

# Data prep
dfo_tows<- dfo_get_tows(here::here("scratch/aja/targets_flow/data/dfo/raw"), out_dir = here::here("scratch/aja/targets_flow/data/dfo/clean"))

dfo_tidy_occu<- dfo_make_tidy_occu(dfo_raw_dir = here::here("scratch/aja/targets_flow/data/dfo/raw"), dfo_tows = dfo_tows, spp_table_dir = here::here("scratch/aja/targets_flow/data/supporting"), out_dir = here::here("scratch/aja/targets_flow/data/dfo/clean"))

nmfs_tows<- nmfs_get_tows(nmfs_raw_dir = here::here("scratch/aja/targets_flow/data/nmfs/raw"), out_dir = here::here("scratch/aja/targets_flow/data/nmfs/clean"))

nmfs_tidy_occu<- nmfs_make_tidy_occu(nmfs_raw_dir = here::here("scratch/aja/targets_flow/data/nmfs/raw"), nmfs_tows = nmfs_tows, spp_table_dir = here::here("scratch/aja/targets_flow/data/supporting"), out_dir = here::here("scratch/aja/targets_flow/data/nmfs/clean"))

all_cov<- bind_nmfs_dfo_cov(nmfs_cov_dir = here::here("scratch/aja/targets_flow/data/nmfs/clean"), dfo_cov_dir = here::here("scratch/aja/targets_flow/data/dfo/clean"), out_dir = here::here("scratch/aja/targets_flow/data/combined"))
all_sample<- bind_nmfs_dfo_sample(nmfs_sample_dir = here::here("scratch/aja/targets_flow/data/nmfs/clean"), dfo_sample_dir = here::here("scratch/aja/targets_flow/data/dfo/clean"), out_dir = here::here("scratch/aja/targets_flow/data/combined"))

year_min_use<- 1985
year_max_use<- 2015
species_code_use<- 73 # Atlantic cod
cov_dat<- vast_make_cov_data(all_cov_dir = here::here("scratch/aja/targets_flow/data/combined"), year_min = year_min_use, year_max = year_max_use, out_dir = here::here("scratch/aja/targets_flow/data/combined"))
samp_dat<- vast_make_sample_data(all_sample_dir = here::here("scratch/aja/targets_flow/data/combined"), year_min = year_min_use, year_max = year_max_use, nmfs_species_code = species_code_use, out_dir = here::here("scratch/aja/targets_flow/data/combined"))

shapefile<- read_polyshape(paste0(shape_path, "NELME_regions/NELME_sf.shp"))

# Input grid and settings
vast_extrap_grid<- vast_make_extrap_grid(shapefile = shapefile, cell_size = 50000)
vast_settings<- vast_make_settings(extrap_grid = vast_extrap_grid, FieldConfig = c("Omega1" = 1, "Epsilon1" = 1, "Omega2" = 1, "Epsilon2" = 1), RhoConfig = c("Beta1" = 3, "Beta2" = 3, "Epsilon1" = 4, "Epsilon2" = 4), bias.correct = FALSE)
vast_coveff<- vast_make_coveff(X1_coveff_vec = c(2, 3, 3, 2, rep(3, nlevels(cov_dat$Year_Cov)-1)), X2_coveff_vec = c(2, 3, 3, 2, rep(3, nlevels(cov_dat$Year_Cov)-1)))

# Initial build
vast_build<- vast_build_sdm(settings = vast_settings, extrap_grid = vast_extrap_grid, samp_dat = samp_dat, cov_dat = cov_dat, X1_formula = ~ Season + Year_Cov, X2_formula = ~ Season + Year_Cov, X_contrasts = list(Season = contrasts(cov_dat$Season, contrasts = FALSE), Year_Cov = contrasts(cov_dat$Year_Cov, contrasts = FALSE)), Xconfig_list = vast_coveff)

# Adjustments
vast_build_adjust<- vast_make_adjustments(vast_build, adjustments = list("log_sigmaXi1_cp" = factor(c(rep(1, 3), rep(4, nlevels(cov_dat$Year_Cov)))), "log_sigmaXi2_cp" = factor(c(rep(1, 3), rep(4, nlevels(cov_dat$Year_Cov))))))

# Final fit
vast_fit<- vast_fit_sdm(vast_build_adjust)
