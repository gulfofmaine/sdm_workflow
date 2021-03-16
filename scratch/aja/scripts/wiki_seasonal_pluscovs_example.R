#####
## Basic example, season and year only
#####
library(VAST)
library(tidyverse)
library(lubridate)
library(raster)
library(splines)
example = load_example( data_set="NWA_yellowtail_seasons" )

# Load data and quick exploration of structure
# Set of years and seasons. The DFO spring survey usually occurs before the NOAA NEFSC spring survey, so ordering accordingly.
year_set = sort(unique(example$sampling_data[,'year']))
season_set = c("DFO", "SPRING", "FALL")

# Create a grid with all unique combinations of seasons and years and then combine these into one "year_season" variable
yearseason_grid = expand.grid("season" = season_set, "year" = year_set)
yearseason_levels = apply(yearseason_grid[,2:1], MARGIN = 1, FUN = paste, collapse = "_")
yearseason_labels = round(yearseason_grid[,'year'] + (as.numeric(factor(yearseason_grid[,'season'], levels = season_set))-1)/length(season_set), digits=1)

# Similar process, but for the observations
yearseason_i = apply(example$sampling_data[,c("year","season")], MARGIN = 1, FUN = paste, collapse = "_")
yearseason_i = factor(yearseason_i, levels = yearseason_levels)

# Add the year_season factor column to our sampling_data data set
example$sampling_data$year_season = yearseason_i
example$sampling_data$season = factor(example$sampling_data$season, levels = season_set)

# Some last processing steps
example$sampling_data = example$sampling_data[, c("year", "season", "year_season", "latitude", "longitude", "swept", "weight")]

# Make dummy observation for each season-year combination
dummy_data = data.frame(
  year = yearseason_grid[,'year'],
  season = yearseason_grid[,'season'],
  year_season = yearseason_levels,
  latitude = mean(example$sampling_data[,'latitude']),
  longitude = mean(example$sampling_data[,'longitude']),
  swept = mean(example$sampling_data[,'swept']),
  weight = 0,
  dummy = TRUE)

# Combine with sampling data
full_data = rbind(cbind(example$sampling_data, dummy = FALSE), dummy_data)

# Create sample data
samp_dat = data.frame(
  "year_season" = as.numeric(full_data$year_season)-1,
  "Lat" = full_data$latitude,
  "Lon" = full_data$longitude,
  "weight" = full_data$weight,
  "Swept" = full_data$swept,
  "Dummy" = full_data$dummy )

# Covariate data. Note here, case sensitive!
cov_dat = data.frame(
  "Year" = as.numeric(full_data$year_season)-1,
  "Year_Cov" = factor(full_data$year, levels = year_set),
  "Season" = full_data$season,
  "Lat" = full_data$latitude,
  "Lon" = full_data$longitude )
# Inspect
table("year_season"=cov_dat$Year, "Actual_year"=cov_dat$Year_Cov)
table("year_season"=cov_dat$Year, "Actual_season"=cov_dat$Season)

#####
## Model settings
#####

# Make settings
settings = make_settings(n_x = 100,
                         Region = example$Region,
                         strata.limits = example$strata.limits,
                         purpose = "index2",
                         FieldConfig = c("Omega1" = 1, "Epsilon1" = 1, "Omega2" = 1, "Epsilon2" = 1),
                         RhoConfig = c("Beta1" = 3, "Beta2" = 3, "Epsilon1" = 4, "Epsilon2" = 4),
                         ObsModel = c(1, 1),
                         bias.correct = FALSE,
                         Options = c('treat_nonencounter_as_zero' = TRUE) )

# Creating model formula
formula_use = ~ Season + Year_Cov

# Implement corner constraint for linear effect but not spatially varying effect:
# * one level for each term is 2 (just spatially varying)
# * all other levels for each term is 3 (spatialy varying plus linear effect)
X1config_cp_use = matrix( c(2, rep(3,nlevels(cov_dat$Season)-1), 2, rep(3,nlevels(cov_dat$Year_Cov)-1) ), nrow=1 )
X2config_cp_use = matrix( c(2, rep(3,nlevels(cov_dat$Season)-1), 2, rep(3,nlevels(cov_dat$Year_Cov)-1) ), nrow=1 )

#####
## Model fit -- make sure to use new functions
#####

fit_orig = fit_model("settings" = settings,
                     "Lat_i" = samp_dat[, 'Lat'],
                     "Lon_i" = samp_dat[, 'Lon'],
                     "t_i" = samp_dat[, 'year_season'],
                     "b_i" = samp_dat[, 'weight'],
                     "a_i" = samp_dat[, 'Swept'],
                     "X1config_cp" = X1config_cp_use,
                     "X2config_cp" = X2config_cp_use,
                     "covariate_data" = cov_dat,
                     "X1_formula" = formula_use,
                     "X2_formula" = formula_use,
                     "X_contrasts" = list(Season = contrasts(cov_dat$Season, contrasts = FALSE), Year_Cov = contrasts(cov_dat$Year_Cov, contrasts = FALSE)),
                     "run_model" = FALSE,
                     "PredTF_i" = samp_dat[, 'Dummy'] )

# Adjust mapping for log_sigmaXi and fitting final model -- pool variance for all seasons and then set year's to NA
Map_adjust = fit_orig$tmb_list$Map

# Pool variances for each term to a single value
Map_adjust$log_sigmaXi1_cp = factor(c(rep(as.numeric(Map_adjust$log_sigmaXi1_cp[1]), nlevels(cov_dat$Season)),
                                      rep(as.numeric(Map_adjust$log_sigmaXi1_cp[nlevels(cov_dat$Season)+1]), nlevels(cov_dat$Year_Cov))))
Map_adjust$log_sigmaXi2_cp = factor(c(rep(as.numeric(Map_adjust$log_sigmaXi2_cp[1]), nlevels(cov_dat$Season)),
                                      rep(as.numeric(Map_adjust$log_sigmaXi2_cp[nlevels(cov_dat$Season)+1]), nlevels(cov_dat$Year_Cov))))

# Fit final model with new mapping
fit  = fit_model("settings" = settings,
                 "Lat_i" = samp_dat[, 'Lat'],
                 "Lon_i" = samp_dat[, 'Lon'],
                 "t_i" = samp_dat[, 'year_season'],
                 "b_i" = samp_dat[, 'weight'],
                 "a_i" = samp_dat[, 'Swept'],
                 "X1config_cp" = X1config_cp_use,
                 "X2config_cp" = X2config_cp_use,
                 "covariate_data" = cov_dat,
                 "X1_formula" = formula_use,
                 "X2_formula" = formula_use,
                 "X_contrasts" = list(Season = contrasts(cov_dat$Season, contrasts = FALSE), Year_Cov = contrasts(cov_dat$Year_Cov, contrasts = FALSE)),
                 "newtonsteps" = 1,
                 "PredTF_i" = samp_dat[, 'Dummy'],
                 "Map" = Map_adjust,
                 "run_model" = FALSE)

######
## Adding other covariates...
######
# Going to want to add a column to the covariate data set.
source(here::here("scratch/aja/scripts/dynamic_2d_extract_func.R"))
dynamic_2d_extract

# OISST stack...
oisst_stack<- raster::stack(paste(shared.path(os.use = "unix", group = "root", folder = "RES_Data/OISST/"), "ThroughFeb2020.grd", sep = ""))
# Annoyingly, dates not preserved...
oisst_ts_data<- seq(from = ymd('1981-09-01'), to = ymd('1981-09-01') + nlayers(oisst_stack)-1, by = 'day')
names(oisst_stack)<- oisst_ts_data

# Need a spatial points object for the tows and to make some column modifications to line up with the existing dynamic_2d_extract function...
cov_sf<- st_as_sf(cov_dat, coords = c("Lon", "Lat"), remove = FALSE, crs = crs(oisst_stack))

# Some name futzing and creation of a date column
cov_sf<- cov_sf %>%
  mutate(., "EST_YEAR" = Year_Cov,
         "EST_MONTH" = case_when(
           Season == "SPRING" ~ 05,
           Season == "DFO" ~ 04,
           Season == "FALL" ~ 09),
         "EST_DAY" = rep(12, nrow(.)),
         "EST_DATE" = as.Date(paste(EST_YEAR, EST_MONTH, EST_DAY, sep = "-")))

cov_dat_new<- dynamic_2d_extract(rast_ts_stack = oisst_stack, stack_name = "sst", t_summ = "seasonal", t_position = NULL, sf_points = cov_sf, df_sf = "df", out_path = here::here("scratch/aja/data"))

# Need to adjust our covariate matrices...
# Creating model formula
formula_use = ~ Season + Year_Cov + bs(sst_seasonal, degree = 2, intercept = FALSE)
#formula_use = ~ Season + Year_Cov + sst_seasonal

# For the season and year bits, we have the same implemented corner constraint for linear effect but not spatially varying effect:
# * one level for each term is 2 (just spatially varying)
# * all other levels for each term is 3 (spatialy varying plus linear effect)

# Then, we need to have something for SST. It looks like for this we just repeat 1 times degrees in the bs splines
X1config_cp_use = matrix( c(2, rep(3,nlevels(cov_dat_new$Season)-1), 2, rep(3,nlevels(cov_dat_new$Year_Cov)-1), rep(1, 2) ), nrow=1 )
X2config_cp_use = matrix( c(2, rep(3,nlevels(cov_dat_new$Season)-1), 2, rep(3,nlevels(cov_dat_new$Year_Cov)-1), rep(1, 2) ), nrow=1 )
# 
# X1config_cp_use = matrix( c(2, rep(3,nlevels(cov_dat_new$Season)-1), 2, rep(3,nlevels(cov_dat_new$Year_Cov)-1), rep(1, 1) ), nrow=1 )
# X2config_cp_use = matrix( c(2, rep(3,nlevels(cov_dat_new$Season)-1), 2, rep(3,nlevels(cov_dat_new$Year_Cov)-1), rep(1, 1) ), nrow=1 )

# Refit...
fit_orig = fit_model("settings" = settings,
                     "Lat_i" = samp_dat[, 'Lat'],
                     "Lon_i" = samp_dat[, 'Lon'],
                     "t_i" = samp_dat[, 'year_season'],
                     "b_i" = samp_dat[, 'weight'],
                     "a_i" = samp_dat[, 'Swept'],
                     "X1config_cp" = X1config_cp_use,
                     "X2config_cp" = X2config_cp_use,
                     "covariate_data" = cov_dat_new,
                     "X1_formula" = formula_use,
                     "X2_formula" = formula_use,
                     "X_contrasts" = list(Season = contrasts(cov_dat$Season, contrasts = FALSE), Year_Cov = contrasts(cov_dat$Year_Cov, contrasts = FALSE)),
                     "run_model" = FALSE,
                     "PredTF_i" = samp_dat[, 'Dummy'] )

# Adjust mapping for log_sigmaXi and fitting final model -- pool variance for all seasons and then set year's to NA
Map_adjust = fit_orig$tmb_list$Map

# Pool variances for each term to a single value
Map_adjust$log_sigmaXi1_cp = factor(c(rep(as.numeric(Map_adjust$log_sigmaXi1_cp[1]), nlevels(cov_dat_new$Season)),
                                      rep(as.numeric(Map_adjust$log_sigmaXi1_cp[nlevels(cov_dat_new$Season)+1]), nlevels(cov_dat$Year_Cov)),
                                      rep(NA, 2)))
Map_adjust$log_sigmaXi2_cp = factor(c(rep(as.numeric(Map_adjust$log_sigmaXi2_cp[1]), nlevels(cov_dat_new$Season)),
                                      rep(as.numeric(Map_adjust$log_sigmaXi2_cp[nlevels(cov_dat_new$Season)+1]), nlevels(cov_dat_new$Year_Cov)),
                                      rep(NA, 2)))
# Map_adjust$log_sigmaXi1_cp = factor(c(rep(as.numeric(Map_adjust$log_sigmaXi1_cp[1]), nlevels(cov_dat_new$Season)),
#                                       rep(as.numeric(Map_adjust$log_sigmaXi1_cp[nlevels(cov_dat_new$Season)+1]), nlevels(cov_dat$Year_Cov)),
#                                       rep(NA, 1)))
# Map_adjust$log_sigmaXi2_cp = factor(c(rep(as.numeric(Map_adjust$log_sigmaXi2_cp[1]), nlevels(cov_dat_new$Season)),
#                                       rep(as.numeric(Map_adjust$log_sigmaXi2_cp[nlevels(cov_dat_new$Season)+1]), nlevels(cov_dat_new$Year_Cov)),
#                                       rep(NA, 1)))

# Fit final model with new mapping
fit_new = fit_model("settings" = settings,
                     "Lat_i" = samp_dat[, 'Lat'],
                     "Lon_i" = samp_dat[, 'Lon'],
                     "t_i" = samp_dat[, 'year_season'],
                     "b_i" = samp_dat[, 'weight'],
                     "a_i" = samp_dat[, 'Swept'],
                     "X1config_cp" = X1config_cp_use,
                     "X2config_cp" = X2config_cp_use,
                     "covariate_data" = cov_dat_new,
                     "X1_formula" = formula_use,
                     "X2_formula" = formula_use,
                     "X_contrasts" = list(Season = contrasts(cov_dat_new$Season, contrasts = FALSE), Year_Cov = contrasts(cov_dat_new$Year_Cov, contrasts = FALSE)),
                    "Map" = Map_adjust,
                     "run_model" = TRUE,
                     "PredTF_i" = samp_dat[, 'Dummy'] )

