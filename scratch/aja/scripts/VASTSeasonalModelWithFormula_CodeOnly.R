#####
## Fitting VAST seasonal model through formula interface
install.packages("devtools")
library(devtools)
install_version("Matrix", version = "1.2.8")
library(Matrix)
install_version("TMB", "1.7.18")
devtools::install_github("James-Thorson-NOAA/FishStatsUtils", ref = "development", force = TRUE)
devtools::install_github("James-Thorson-NOAA/VAST", ref = "development", force = TRUE)
library(VAST)
library(FishStatsUtils)
library(tidyverse)
library(here)

# # Older version of TMB to avoid rtweedie issues (https://github.com/James-Thorson-NOAA/VAST/issues/276#issue-804990507)
# install_version("TMB", version = "1.7.18", repos = "http://cran.us.r-project.org")
# library(TMB)

# Making sure new functions are being used
source("~/Github/FishStatsUtils/R/make_covariates.R")
source("~/Github/FishStatsUtils/R/fit_model.R")
source("~/Github/VAST/R/make_data.R")

#####
## Data processing
#####
# Load data and quick exploration of structure
load("~/GitHub/VAST/R/VAST_SeasonalExampleData.RData")

# Set of seasons and years. The DFO spring survey usually occurs before the NOAA NEFSC spring survey, so ordering accordingly.
season_set<- c("DFO", "SPRING", "FALL")
year_set<- sort(unique(svdata[,'year']))

# Create a grid with all unique combinations of seasons and years and then combine these into one "season_year" variable
seasonyear_grid<- expand.grid("season" = season_set, "year" = year_set)
t_levels<- apply(seasonyear_grid, MARGIN = 1, FUN = paste, collapse = "_")
t_labels<- round(seasonyear_grid[,'year'] + (as.numeric(factor(seasonyear_grid[,'season'], levels = season_set))-1)/length(season_set), digits=1)

# Similar process, but for the observations
seasonyear_i<- apply(svdata[,c("season","year")], MARGIN = 1, FUN = paste, collapse = "_")
seasonyear_i<- factor(seasonyear_i, levels = t_levels)

# Add the season_year factor column to our sampling_data data set
svdata <- cbind(svdata, "season_year" = seasonyear_i)
svdata$season<- factor(svdata$season, levels = season_set)

# Some last processing steps
svdata<- svdata[, c("year", "season", "season_year", "latitude", "longitude", "swept", "weight")]
colnames(svdata)<- c("year", "season", "season_year", "latitude", "longitude", "swept", "response")

# Adding in dummy data for missing season_years
# Make dummy observation for each season-year combination
dummy_data<- data.frame(year = seasonyear_grid[,'year'], season = seasonyear_grid[,'season'], season_year = t_levels, latitude = mean(svdata[,'latitude']), longitude = mean(svdata[,'longitude']), swept = mean(svdata[,'swept']), response = 0, dummy = TRUE)

# Combine with sampling data
svdata<- rbind(cbind(svdata, dummy = FALSE), dummy_data)

# Create sample data
samp_dat<- data.frame("spp" = rep("Yellowtail", nrow(svdata)), "Year" = as.numeric(svdata$season_year)-1, "Season" = svdata$season, "Season_Year" = svdata$season_year, "Lat" = svdata$latitude, "Lon" = svdata$longitude, "Response" = svdata$response, "Swept" = svdata$swept, "Dummy" = svdata$dummy)

# Covariate data. Note here, case sensitive!
cov_dat<- data.frame("spp" = rep("Yellowtail", nrow(svdata)), "Year" = as.numeric(svdata$season_year)-1, "Year_Cov" = factor(svdata$year, levels = year_set), "Season" = svdata$season, "Lat" = svdata$latitude, "Lon" = svdata$longitude)

#####
## Model settings
#####
# Extrapolation grid
maximum_distance_from_sample<- 10
grid_dim_km<- c(5, 5)
region_code<- "northwest_atlantic"
strata.limits<- list("Yellowtail" = c(1130, 1140, 1150, 1160, 1170, 1180, 1190, 1200, 1210))

# Observation model -- "Poisson-link delta model" with lognormal positive catch rate distribution
ObsModel<- c(1, 1)

# Number of spatial and spatio-temporal factors to use for each linear predictor. Here, using a single species model, and turning spatial and spatio-temporal variability effects for both linear predictors "on".
FieldConfig<- c("Omega1" = 1, "Epsilon1" = 1, "Omega2" = 1, "Epsilon2" = 1)

# Setting the structure of parameters across season_year surveys. Here, setting intercepts for the first and second linear predictors (Beta) to be estimated as random effects, which are constant across season_years, and then estimating a first order auto-regressive structure for the spatio-temporal variability (Epsilon) of both linear predictors. This implementation facilitates estimating occurrence at unsampled times/areas.
RhoConfig<- c("Beta1" = 3, "Beta2" = 3, "Epsilon1" = 4, "Epsilon2" = 4)
Options<- c('treat_nonencounter_as_zero' = TRUE)

# OutFile
OutFile<- "~/Desktop/VASTIntraAnnualModel/"
if(!file.exists(OutFile)){
  dir.create(OutFile)
}

# Set the run directory
RunDir<- OutFile

# Make settings
settings<- make_settings(n_x = 100, Region = region_code, purpose = "index2", FieldConfig = FieldConfig, RhoConfig = RhoConfig, ObsModel = ObsModel, bias.correct = FALSE, strata.limits = strata.limits, Options = Options)

# Creating model formula
formula_use<- ~ Season + Year_Cov

# Making these to match what was in original analysis
# Empty at first
season_effect<- rep(NA, length(levels(cov_dat$Season)))
year_cov_effect<- rep(NA, length(unique(cov_dat$Year_Cov)))

X1config_cp_use<- matrix(data = c(season_effect, year_cov_effect), nrow = 1)
X2config_cp_use<- matrix(data = c(season_effect, year_cov_effect), nrow = 1)

# Populate manually
X1config_cp_use[1,]<- c(2, 3, 3, 0, rep(1, 32))
X2config_cp_use[1,]<- c(0, 3, 3, 0, rep(1, 32))

#####
## Model fit -- make sure to use new functions
#####

# Fit model -- Within fit_model and make_data, the make_covariates call is currently attributed to FishStatsUtils as: "FishStatsUtils::make_covariates". This means the call is not using the edited "make_covariates" function. Not sure the best way around this, but potentially dealing with it through `assignInNamespace`?
make_covariates_aja<- source("~/Github/FishStatsUtils/R/make_covariates.R")[[1]]
assignInNamespace("make_covariates", make_covariates_aja, ns = "FishStatsUtils")

Use_REML<- FALSE
run_model_use<- FALSE
fit_seas_form<- fit_model("settings" = settings, "Lat_i" = samp_dat[, 'Lat'], "Lon_i" = samp_dat[, 'Lon'], "t_i" = samp_dat[, 'Year'], "c_i" = rep(0, nrow(samp_dat)), "b_i" = samp_dat[, 'Response'], "a_i" = samp_dat[, 'Swept'], "X1config_cp" = X1config_cp_use, "X2config_cp" = X2config_cp_use, "covariate_data" = cov_dat, "X1_formula" = formula_use, "X2_formula" = formula_use, Xcontrasts = list(Season = contrasts(cov_dat$Season, contrasts = FALSE), Year_Cov = contrasts(cov_dat$Year_Cov, contrasts = FALSE)), "newtonsteps" = 1, "getsd" = TRUE, "getReportCovariance" = TRUE, "observations_LL" = cbind("Lat" = samp_dat[, 'Lat'], "Lon" = samp_dat[, 'Lon']), "maximum_distance_from_sample" = maximum_distance_from_sample, "grid_dim_km" = grid_dim_km, "run_model" = run_model_use, "test_fit" = FALSE, working_dir = RunDir, "PredTF_i" = samp_dat[, 'Dummy'], "Use_REML" = Use_REML, "getJointPrecision" = FALSE)

# Adjust mapping for log_sigmaXi and fitting final model -- pool variance for all seasons and then set year's to NA
Map_adjust<- fit_seas_form$tmb_list$Map

# Need to think about how to make this more robust and less hard-wirey...
Map_adjust$log_sigmaXi1_cp<- factor(c(rep(as.numeric(Map_adjust$log_sigmaXi1_cp[1]), nlevels(cov_dat$Season)), rep(NA, nlevels(cov_dat$Year_Cov))))
Map_adjust$log_sigmaXi2_cp<- factor(c(rep(as.numeric(Map_adjust$log_sigmaXi2_cp[1]), nlevels(cov_dat$Season)), rep(NA, nlevels(cov_dat$Year_Cov))))

# Fit final model with new mapping
run_model_use<- TRUE
fit_seas_form<- fit_model("settings" = settings, "Lat_i" = samp_dat[, 'Lat'], "Lon_i" = samp_dat[, 'Lon'], "t_i" = samp_dat[, 'Year'], "c_i" = rep(0, nrow(samp_dat)), "b_i" = samp_dat[, 'Response'], "a_i" = samp_dat[, 'Swept'], "X1config_cp" = X1config_cp_use, "X2config_cp" = X2config_cp_use, "covariate_data" = cov_dat, "X1_formula" = formula_use, "X2_formula" = formula_use, Xcontrasts = list(Season = contrasts(cov_dat$Season, contrasts = FALSE), Year_Cov = contrasts(cov_dat$Year_Cov, contrasts = FALSE)), Map = Map_adjust, "newtonsteps" = 1, "getsd" = TRUE, "getReportCovariance" = TRUE, "observations_LL" = cbind("Lat" = samp_dat[, 'Lat'], "Lon" = samp_dat[, 'Lon']), "maximum_distance_from_sample" = maximum_distance_from_sample, "grid_dim_km" = grid_dim_km, "run_model" = run_model_use, "test_fit" = FALSE, working_dir = RunDir, "PredTF_i" = samp_dat[, 'Dummy'], "Use_REML" = Use_REML, "getJointPrecision" = FALSE)

