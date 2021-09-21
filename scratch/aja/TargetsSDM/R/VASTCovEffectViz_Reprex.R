#####
## Visualizing covariate effects when using a mix of factor and continuous, spatially-varying covariates
#####

library(VAST)
library(splines)  # Used to include basis-splines
library(effects)  # Used to visualize covariate effects
library(pdp)

#####
## 1. Continuous covariate only -- this essentially mirrors the "Visualizing covariate effects" wiki 
#####

# load data set
# see `?load_example` for list of stocks with example data
# that are installed automatically with `FishStatsUtils`.
example = load_example( data_set="covariate_example" )

# Make settings (turning off bias.correct to save time for example)
settings = make_settings( n_x=100,
                          Region=example$Region,
                          purpose="index2",
                          use_anisotropy=FALSE,
                          bias.correct=FALSE,
                          fine_scale=TRUE )

# Define formula, here using a b-spline basis with the splines library. For the cubic regression spline, we want to set "degree = 3"
X1_formula = ~ bs( BOT_DEPTH, degree=3, intercept=FALSE) 
X2_formula = ~ bs( BOT_DEPTH, degree=3, intercept=FALSE)

# Bottom depth is "static" (not changing among years), so setting Year = NA to cause values to be duplicated internally for all values of Year
example$covariate_data[,'Year'] = NA

# Scale/center bottom depth (for stability and original comparison with mgcv::gam)
example$covariate_data[,'BOT_DEPTH']<- as.numeric(scale(example$covariate_data[,'BOT_DEPTH']))

# Adjusting ObsModel from default for "purpose = index2", which uses "ObsModel = c(2, 1)" and specifies the "Poisson" link delta model, to use "ObsModel = c(1, 0)". Done originally for comparison with mgcv::gam
settings$ObsModel<- c(1, 0)

# Run model
fit_base_example = fit_model( "settings" = settings,
                              Lat_i = example$sampling_data[,'Lat'],
                              Lon_i = example$sampling_data[,'Lon'],
                              t_i = example$sampling_data[,'Year'],
                              b_i = example$sampling_data[,'Catch_KG'],
                              a_i = example$sampling_data[,'AreaSwept_km2'],
                              X1_formula = X1_formula,
                              X2_formula = X2_formula,
                              covariate_data = example$covariate_data )

## Visualizing with `effects`
# Must add dataframes to global environment
covariate_data_full = fit_base_example$effects$covariate_data_full
catchability_data_full = fit_base_example$effects$catchability_data_full

# Plot 1st linear predictor
pred = Effect.fit_model( fit_base_example,
                         focal.predictors = c("BOT_DEPTH"),
                         which_formula = "X1", 
                         xlevels = 100)
plot(pred) # Looks good

## Visualizing with 'pdp'
# Just a note here, the pdp implementation following the visualizing covariates wiki code will not work for this example. I've done some digging, and plan to come back to this at another time. 

#####
## 2. Factor and a continuous, spatially varying covariate
#####
# load data set
# see `?load_example` for list of stocks with example data
# that are installed automatically with `FishStatsUtils`.
example = load_example( data_set="covariate_example" )

# What do we have available?
summary(example$sampling_data)

# Year would be the logical thing to include. There are a variety of issues with this.
# Isssue 1: During the model set up, `make_covariates` will complain that there are missing years because the survey is done once every twice year. So, some adjusting to force the years to be successive annually.
# Reloading dataset
example = load_example( data_set="covariate_example" )
example$sampling_data<- data.frame(example$sampling_data)

# Converting years to be sequential and creating a "Year_Covariate" factor to include in the model formula
example$sampling_data[,'Year']<- as.numeric(factor(example$sampling_data[,'Year']))
example$sampling_data<- data.frame(example$sampling_data, "Year_Cov" = factor(example$sampling_data[,'Year'], levels = unique(example$sampling_data[,'Year'])))

example$covariate_data[,'Year']<- as.numeric(factor(example$covariate_data[,'Year']))
example$covariate_data<- data.frame(example$covariate_data, "Year_Cov" = factor(example$covariate_data[,'Year'], levels = unique(example$covariate_data[,'Year'])))

example$covariate_data[,'BOT_DEPTH']<- as.numeric(scale(example$covariate_data[,'BOT_DEPTH']))

# Make settings (turning off bias.correct to save time for example)
settings = make_settings( n_x=100,
                          Region=example$Region,
                          purpose="index2",
                          use_anisotropy=FALSE,
                          bias.correct=FALSE,
                          fine_scale=TRUE )
settings$ObsModel<- c(1, 0)

# Model formula
X1_formula<- ~ Year_Cov + bs( BOT_DEPTH, degree=3, intercept=FALSE)
X2_formula<- ~ Year_Cov + bs( BOT_DEPTH, degree=3, intercept=FALSE)

# Issue 2: If we just proceed as normal, VAST will generate a beta1/2_ft parameter for each year, as well as gamma1/2_cp parameters for each year (-1, given null factor contrasts) and then the three bs::spline parameters. Of course, that is going to be an issue. I *think* we would be after a model that had one beta1/2_ft parameter, which is a global intercept to capture average occurrence across entire domain/time series, then have a gamma1/2_cp parameter for each year and the three bs::splines parameters. I originally thought we would adjust the 'X_contrasts' argument to FALSE to ensure the first year isn't set as a reference level. This results in a singularity issue, which I guess isn't all that surprising as it is unable to estimate the beta1/2_fts with the each level of the Year_Cov factor? So, leaving the contrasts as is, which means beta1/2_ft is really the first level in Year_Cov. 
# Build base model that we can then adjust as needed. Had originally thought about mapping off the beta1/2_fts, which isn't possible.
fit_base = fit_model( "settings" = settings,
                      Lat_i = example$sampling_data[,'Lat'],
                      Lon_i = example$sampling_data[,'Lon'],
                      t_i = example$sampling_data[,'Year'],
                      b_i = example$sampling_data[,'Catch_KG'],
                      a_i = example$sampling_data[,'AreaSwept_km2'],
                      X1_formula = X1_formula,
                      X2_formula = X2_formula,
                      covariate_data = example$covariate_data,
                      run_model = FALSE)

# Since the beta1/2_fts are generated based on the time vector (t_i), what if we fix all of the years to the same value, while retaining the "Year_Cov" to soak up the year to year variability? 
example$sampling_data[,'Year']<- rep(1, nrow(example$sampling_data))
example$covariate_data[,'Year']<- rep(1, nrow(example$covariate_data))

fit_adjust1 = fit_model( "settings" = settings,
                         Lat_i = example$sampling_data[,'Lat'],
                         Lon_i = example$sampling_data[,'Lon'],
                         t_i = example$sampling_data[,'Year'],
                         b_i = example$sampling_data[,'Catch_KG'],
                         a_i = example$sampling_data[,'AreaSwept_km2'],
                         X1_formula = X1_formula,
                         X2_formula = X2_formula,
                         covariate_data = example$covariate_data,
                         run_model = FALSE)

# Looks good...run it 
fit_adjust1 = fit_model( "settings" = settings,
                         Lat_i = example$sampling_data[,'Lat'],
                         Lon_i = example$sampling_data[,'Lon'],
                         t_i = example$sampling_data[,'Year'],
                         b_i = example$sampling_data[,'Catch_KG'],
                         a_i = example$sampling_data[,'AreaSwept_km2'],
                         X1_formula = X1_formula,
                         X2_formula = X2_formula,
                         covariate_data = example$covariate_data,
                         run_model = TRUE)

# Warning about the Hessian. For future reference, the Hessian matrix is central to calculating the standard errors (matrix of second derivitives at the MLE). This warning is signaling that at the MLE, the curvature of the negative log-likelihood surface is problematic --  downward-curving, or flat, in some direction -- rather than what we would hope for, which would look something like a well defined trough. I was really hoping to just ignore that to get to the actual plotting error. Of course, no such luck as the object needed in the Effects function is empty when Hessian is not positive definite. This can occur because a parameter is hitting a bound, or maybe there is not enough data to estimate a parameter.
# Checking estimatability
fit_adjust1_check = fit_model( "settings" = settings,
                         Lat_i = example$sampling_data[,'Lat'],
                         Lon_i = example$sampling_data[,'Lon'],
                         t_i = example$sampling_data[,'Year'],
                         b_i = example$sampling_data[,'Catch_KG'],
                         a_i = example$sampling_data[,'AreaSwept_km2'],
                         X1_formula = X1_formula,
                         X2_formula = X2_formula,
                         covariate_data = example$covariate_data,
                         run_model = TRUE,
                         get_sd = FALSE,
                         newtonsteps = 0)
TMBhelper::check_estimability(fit_adjust1_check$tmb_list$Obj)

# Turn off spatio and spatio-temporal variation. I originally just switched these off for the second linear predictor, but then had issues estimating them for the first linear predictor, too. So, turning them all off.
settings_adjust<- settings
settings_adjust[['FieldConfig']]<- c("Omega1" = 0, "Epsilon1" = 0, "Omega2" = 0, "Epsilon2" = 0)
fit_adjust2_check = fit_model( "settings" = settings_adjust,
                               Lat_i = example$sampling_data[,'Lat'],
                               Lon_i = example$sampling_data[,'Lon'],
                               t_i = example$sampling_data[,'Year'],
                               b_i = example$sampling_data[,'Catch_KG'],
                               a_i = example$sampling_data[,'AreaSwept_km2'],
                               X1_formula = X1_formula,
                               X2_formula = X2_formula,
                               covariate_data = example$covariate_data,
                               run_model = TRUE,
                               get_sd = FALSE,
                               newtonsteps = 0)
TMBhelper::check_estimability(fit_adjust2_check$tmb_list$Obj)

# Looks good, run it
fit_adjust2_run = fit_model( "settings" = settings_adjust,
                               Lat_i = example$sampling_data[,'Lat'],
                               Lon_i = example$sampling_data[,'Lon'],
                               t_i = example$sampling_data[,'Year'],
                               b_i = example$sampling_data[,'Catch_KG'],
                               a_i = example$sampling_data[,'AreaSwept_km2'],
                               X1_formula = X1_formula,
                               X2_formula = X2_formula,
                               covariate_data = example$covariate_data,
                               run_model = TRUE)

## Visualizing effects with 'effects' package
# Must add data-frames to global environment
covariate_data_full = fit_adjust2_run$effects$covariate_data_full
catchability_data_full = fit_adjust2_run$effects$catchability_data_full

# Plot 1st linear predictor
pred = Effect.fit_model( fit_adjust2_run,
                         focal.predictors = c("BOT_DEPTH"),
                         which_formula = "X1", 
                         xlevels = 100)

# Error in mod.matrix %*% scoef : non-conformable arguments -- same as seasonal model error.
