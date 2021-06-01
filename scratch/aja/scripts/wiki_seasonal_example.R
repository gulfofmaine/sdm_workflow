#####
## VAST seasonal example and predictions
#####

### Libraries
library(VAST)
library(sf)


### Data loading and processing
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
VAST_Testing
setwd("~/Desktop/VAST_Testing")
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
# This will prevent refitting the model and just load it if its been fit already
first = FALSE 
if(first){
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
                   "run_model" = TRUE)
  saveRDS(fit, here::here("", "VASTseasonalexample_fit.rds"))
} else {
  fit = readRDS("~/GitHub/sdm_workflow/VASTseasonalexample_fit.rds")
}

#####
## Model predictions using predict.fit_model
#####
# Make predictions -- pretend like the last year of obs is the "next" year we want to predict...
pred_dat = cov_dat
pred_dat = pred_dat[pred_dat$Year_Cov == "2017",]

# Adjusting "Year" and "Year_Cov"
pred_dat$Year = pred_dat$Year + 1
pred_dat$Year_Cov = "2018"
pred_dat$Year_Cov = factor(pred_dat$Year_Cov, levels = c(levels(cov_dat$Year_Cov), "2018"))

# Check it out
summary(pred_dat)

# Projected coordinates...
pred_dat_sf<- st_as_sf(pred_dat, coords = c("Lon", "Lat"), remove = FALSE, crs = 4326)

# Test to see if predict crashes?
test_crash = FALSE
# Prediction....
if(test_crash){
  vast_pred = predict(x = fit, what = "D_i", Lat_i = pred_dat$Lat, Lon_i = pred_dat$Lon, t_i = pred_dat$Year, a_i = rep(0.03, nrow(pred_dat)), new_covariate_data = pred_dat, new_catchability_data = NULL, do_checks = TRUE, working_dir = paste0(getwd(),"/"))
  
  # Contrast error. Can `predict` accept this as an additional argument?
  vast_pred = predict(x = fit, what = "D_i", Lat_i = pred_dat$Lat, Lon_i = pred_dat$Lon, t_i = pred_dat$Year, a_i = rep(0.03, nrow(pred_dat)), new_covariate_data = pred_dat, X_contrasts = list(Season = contrasts(pred_dat$Season, contrasts = FALSE), Year_Cov = contrasts(pred_dat$Year_Cov, contrasts = FALSE)), new_catchability_data = NULL, do_checks = TRUE, working_dir = paste0(getwd(),"/"))
  
  # Not as it stands now. Edited function to accept this as well as our "CP" matrix bits that can then pass through to `make_data`
  source("~/GitHub/sdm_workflow/scratch/aja/TargetsSDM/R/vast_functions.R")
  vast_pred = predict.fit_model_aja(x = fit, what = "D_i", Lat_i = pred_dat$Lat, Lon_i = pred_dat$Lon, t_i = pred_dat$Year, a_i = rep(0.03, nrow(pred_dat)), new_covariate_data = pred_dat, Xcontrasts_pred = list(Season = contrasts(pred_dat$Season, contrasts = FALSE), Year_Cov = contrasts(pred_dat$Year_Cov, contrasts = FALSE)), X1config_cp_pred = matrix(c(2, rep(3, nlevels(pred_dat$Season)-1), 2, rep(3, nlevels(pred_dat$Year_Cov)-1)), nrow = 1), X2config_cp_pred = matrix(c(2, rep(3,nlevels(pred_dat$Season)-1), 2, rep(3, nlevels(pred_dat$Year_Cov)-1) ), nrow = 1))
  
  # Crash. A little detective work and this crash occurs at the `dyn.load` line in `make_model`. I have not really found a great way around this yet...
}


