#####
## Visualizing covariate effects when using a mix of factor and continuous covariates
#####

library(VAST)
library(splines)  # Used to include basis-splines
library(effects)  # Used to visualize covariate effects
library(pdp)

#####
## 1. Continuous covariate only. 
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

# Adjusting ObsModel from default for "purpose = index2", which uses "ObsModel = c(2, 1)" and specifies the "Poisson" link delta model, to use "ObsModel = c(1, 0)". Again, done originally for comparison with mgcv::gam
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
# Must add data-frames to global environment
covariate_data_full = fit_base_example$effects$covariate_data_full
catchability_data_full = fit_base_example$effects$catchability_data_full

# Plot 1st linear predictor
pred = Effect.fit_model( fit_base_example,
                         focal.predictors = c("BOT_DEPTH"),
                         which_formula = "X1", 
                         xlevels = 100)
plot(pred) # Looks good


## Visualizing with `pdp`
# Make function to interface with pdp
pred.fun = function( object, newdata ){
  predict( x=object,
           Lat_i = object$data_frame$Lat_i,
           Lon_i = object$data_frame$Lon_i,
           t_i = object$data_frame$t_i,
           a_i = object$data_frame$a_i,
           what = "P1_iz",
           new_covariate_data = newdata,
           do_checks = FALSE )
}

# Run partial - this returns an error. Fairly certain it is arising because of the NAs for year in our covariate_data.
Partial = pdp::partial( object = fit_base_example,
                        pred.var = "BOT_DEPTH",
                        pred.fun = pred.fun,
                        train = fit_example$covariate_data) 

# Reloading dataset
example = load_example( data_set="covariate_example" )
example$sampling_data<- data.frame(example$sampling_data)

# Converting years to be sequential...
example$sampling_data[,'Year']<- as.numeric(factor(example$sampling_data[,'Year']))
example$covariate_data[,'Year']<- as.numeric(factor(example$covariate_data[,'Year']))

# Refitting

# In the Wiki example, covariates are rescaled to have an SD >0.1 and <10 (for numerical stability) with the following line. The mgcv::gam function I believe automatically scales/centers any continuous covariates. So, rather than doing the rescaling as in the wiki example, going to scale/center instead for direct comparison between the two approaches.
# example$covariate_data[,'BOT_DEPTH'] = example$covariate_data[,'BOT_DEPTH'] / 100
example$covariate_data[,'BOT_DEPTH']<- as.numeric(scale(example$covariate_data[,'BOT_DEPTH']))

# Another thing we need to adjust is the default settings and the ObsModel default for "purpose = index2", which uses "ObsModel = c(2, 1)" and specifies the "Poisson" link delta model. I'm not entirely sure how to implement that in mgcv::gam. So, instead going to use "ObsModel = c(1, 0)", which is the "traditional" two-stage delta log normal GAM.
settings$ObsModel<- c(1, 0)

# Run model
fit_example = fit_model( "settings" = settings,
                         Lat_i = example$sampling_data[,'Lat'],
                         Lon_i = example$sampling_data[,'Lon'],
                         t_i = example$sampling_data[,'Year'],
                         b_i = example$sampling_data[,'Catch_KG'],
                         a_i = example$sampling_data[,'AreaSwept_km2'],
                         X1_formula = X1_formula,
                         X2_formula = X2_formula,
                         covariate_data = example$covariate_data )

# Redoing with effects...
# Must add data-frames to global environment
covariate_data_full = fit_example$effects$covariate_data_full
catchability_data_full = fit_example$effects$catchability_data_full

# Plot 1st linear predictor
pred = Effect.fit_model( fit_example,
                         focal.predictors = c("BOT_DEPTH"),
                         which_formula = "X1", 
                         xlevels = 100)
plot(pred)

# PDP?
pred.fun = function( object, newdata ){
  predict( x=object,
           Lat_i = object$data_frame$Lat_i,
           Lon_i = object$data_frame$Lon_i,
           t_i = object$data_frame$t_i,
           a_i = object$data_frame$a_i,
           what = "P1_iz",
           new_covariate_data = newdata,
           do_checks = FALSE )
}

# Run partial - this returns the same NA error. 
Partial = pdp::partial( object = fit_example,
                        pred.var = "BOT_DEPTH",
                        pred.fun = pred.fun,
                        train = fit_example$covariate_data)
str(Partial)
autoplot(Partial) # Yikes...

partial_plot_dat<- data.frame(Partial) %>%
  arrange(., BOT_DEPTH)
str(partial_plot_dat)

# A whole lot of predictions for the same BOT_DEPTH value. 6158 to be exact -- this is the number of observations in fit_example$covariate_data. 
t<- partial_plot_dat[partial_plot_dat$BOT_DEPTH == max(partial_plot_dat$BOT_DEPTH),]

# Average at each depth value and then plot??
partial_plot_dat_avg<- partial_plot_dat %>%
  group_by(., BOT_DEPTH) %>%
  summarize(., "AvgPred" = mean(yhat))

ggplot() +
  geom_line(data = partial_plot_dat_avg, aes(x = BOT_DEPTH, y = AvgPred))

# Stepping aside from that for a second..back to just the effects one that was working.

#####
## 2. Factor and a continuous, spatially varying covariate
#####
# load data set
# see `?load_example` for list of stocks with example data
# that are installed automatically with `FishStatsUtils`.
example = load_example( data_set="covariate_example" )

# What do we have available?
summary(example$sampling_data)

# Year would be the logical thing to include. There are a variety of issues with this...
  # During the model set up, `make_covariates` will complain that there are missing years because the survey is done once every twice year. So, we will need to do the same trick we did above and just make sure the years are successive.
# Reloading dataset
example = load_example( data_set="covariate_example" )
example$sampling_data<- data.frame(example$sampling_data)

# Converting years to be sequential and creating a "Year_Covariate" factor to include in the model formula
example$sampling_data[,'Year']<- as.numeric(factor(example$sampling_data[,'Year']))
example$sampling_data<- bind_cols(example$sampling_data, "Year_Cov" = factor(example$sampling_data[,'Year'], levels = unique(example$sampling_data[,'Year'])))

example$covariate_data[,'Year']<- as.numeric(factor(example$covariate_data[,'Year']))
example$covariate_data<- bind_cols(example$covariate_data, "Year_Cov" = factor(example$covariate_data[,'Year'], levels = unique(example$covariate_data[,'Year'])))

example$covariate_data[,'BOT_DEPTH']<- as.numeric(scale(example$covariate_data[,'BOT_DEPTH']))

# Make settings (turning off bias.correct to save time for example)
settings = make_settings( n_x=100,
                          Region=example$Region,
                          purpose="index2",
                          use_anisotropy=FALSE,
                          bias.correct=FALSE,
                          fine_scale=TRUE )
# Another thing we need to adjust is the default settings and the ObsModel default for "purpose = index2", which uses "ObsModel = c(2, 1)" and specifies the "Poisson" link delta model. I'm not entirely sure how to implement that in mgcv::gam. So, instead going to use "ObsModel = c(1, 0)", which is the "traditional" two-stage delta log normal GAM.
settings$ObsModel<- c(1, 0)

# Model formula
X1_formula<- ~ Year_Cov + bs( BOT_DEPTH, degree=3, intercept=FALSE)
X2_formula<- ~ Year_Cov + bs( BOT_DEPTH, degree=3, intercept=FALSE)

  # During the model building if we include year as a factor covariate (e.g., X1_formula = Year_Cov + bs(BOT_DEPTH, degree=3, intercept = FALSE)), we are also going to run into issues as internally, VAST will assume we want to estimate a parameter for each Year_Cov levels. We do want to do this BUT without any other modifications the Year_Cov parameters will be completely confounded by the `beta1/2_fts` (or yearly intercepts) that VAST internally generates. This will cause a computationally singular error. 
show_problems<- TRUE
if(show_problems){
  
  
  # Run model -- set run_model to FALSE, otherwise singular system error, which makes sense as Year_Cov is completely confounded with the betas..
  fit_base = fit_model( "settings" = settings,
                        Lat_i = example$sampling_data[,'Lat'],
                        Lon_i = example$sampling_data[,'Lon'],
                        t_i = example$sampling_data[,'Year'],
                        b_i = example$sampling_data[,'Catch_KG'],
                        a_i = example$sampling_data[,'AreaSwept_km2'],
                        X1_formula = X1_formula,
                        X2_formula = X2_formula,
                        covariate_data = example$covariate_data,
                        run_model = TRUE)
  
  # We can see why this is happening by looking at the beta1/2_fts and the gamma1/2_cps. Recall that during the first example, we had 9 beta1/2_fts and 3 gamma1/2_cps. Now, we still have the 9 beta1/2_fts, but we also have 11 of the gamma1/2cps, which corresponds to the levels in our year factor (minus 1, as the first year is set to the reference level within the default contrasts behavior) and then the 3 bs::splines parameters. 
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
}

# How do we move forward? Thinking about the options. 
  # Option 1: One beta1_ft parameter and then each year as a gamma1/2_cp parameter, plus the three for the GAM? What would this imply? Intercept would capture the average occurrence of a species across the entire spatial domain and timeseries. That seems...logical? With the years in gamma, though, the default behavior of contrasts would drop the first level as the reference level. We could adjust that with `X_contrasts`?
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
fit_adjust_map<- fit_base$tmb_list$Map
names(fit_adjust_map)

# Doesn't look like we can map off the intercept terms easily. Maybe there is an easier way to do this. If the beta1/2_fts are generated based on the time vector (t_i), we could just fix all of the years to the same value, while retaining the "Year_Cov" to soak up the year to year variability? I tried this with setting year to one numeric value, it didn't like that -- hessian error. Then tried year = NA, didn't like that either. 
example$sampling_data[,'Year']<- rep(1, nrow(example$sampling_data))
example$covariate_data[,'Year']<- rep(1, nrow(example$covariate_data))

# Settings
settings$ObsModel<- c(1, 0)

# Model formula
X1_formula<- ~ Year_Cov + bs( BOT_DEPTH, degree=3, intercept=FALSE)
X2_formula<- ~ Year_Cov + bs( BOT_DEPTH, degree=3, intercept=FALSE)

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

# Almost there -- we are still missing a year effect for each year as evidenced by the number of gamma1_cp parameters. To adjust this, need to adjust the Xcontrasts bit, which by default set the first level of a factor as the base level. 
fit_adjust2 = fit_model( "settings" = settings,
                         Lat_i = example$sampling_data[,'Lat'],
                         Lon_i = example$sampling_data[,'Lon'],
                         t_i = example$sampling_data[,'Year'],
                         b_i = example$sampling_data[,'Catch_KG'],
                         a_i = example$sampling_data[,'AreaSwept_km2'],
                         X1_formula = X1_formula,
                         X2_formula = X2_formula,
                         X_contrasts = list(Year_Cov = contrasts(example$covariate_data$Year_Cov, contrasts = FALSE)),
                         covariate_data = example$covariate_data,
                         run_model = FALSE)

# Alright, all looks good with now 12 gamma1/2_cps. Switch run_model to TRUE and re run...
fit_adjust2 = fit_model( "settings" = settings,
                         Lat_i = example$sampling_data[,'Lat'],
                         Lon_i = example$sampling_data[,'Lon'],
                         t_i = example$sampling_data[,'Year'],
                         b_i = example$sampling_data[,'Catch_KG'],
                         a_i = example$sampling_data[,'AreaSwept_km2'],
                         X1_formula = X1_formula,
                         X2_formula = X2_formula,
                         X_contrasts = list(Year_Cov = contrasts(example$covariate_data$Year_Cov, contrasts = FALSE)),
                         covariate_data = example$covariate_data,
                         run_model = TRUE)

# Ugh, warning about the Hessian. Really hoping to ignore that right now and see if we can plot covariate effects. Of course, no such luck as the object needed in the Effects function is empty when Hessian is not positive definite. This can occur because a parameter is hitting a bound, or maybe there is not enough data to estimate a parameter...
table(example$covariate_data$Year_Cov)
samp_dat_check<- example$sampling_data
samp_dat_check$Presence<- ifelse(example$sampling_data$Catch_KG > 0, "Present", "Absent")
table(samp_dat_check$Year_Cov, samp_dat_check$Presence)

fit_base_example$parameter_estimates$diagnostics
fit_adjust2$parameter_estimates$opt$diagnostics

# Maybe we need the intercept to actually be the first year? And go back to regular contrasts?
fit_adjust3_map<- fit_adjust2$tmb_list$Map
fit_adjust3_map$gamma1_cp[1]<- factor(c(NA, seq(1:)
fit_adjust3_map$gamma2_cp[1]<- NA

fit_adjust3_params<- fit_adjust3$tmb_list$Parameters
fit_adjust3_params$gamma1_cp<- matrix(fit_adjust3_params$gamma1_cp[,-1], nrow = 1)
fit_adjust3_params$gamma2_cp<- matrix(fit_adjust3_params$gamma2_cp[,-1], nrow = 1)


fit_adjust3<- fit_model( "settings" = settings,
                         Lat_i = example$sampling_data[,'Lat'],
                         Lon_i = example$sampling_data[,'Lon'],
                         t_i = example$sampling_data[,'Year'],
                         b_i = example$sampling_data[,'Catch_KG'],
                         a_i = example$sampling_data[,'AreaSwept_km2'],
                         X1_formula = X1_formula,
                         X2_formula = X2_formula,
                         #X_contrasts = list(Year_Cov = contrasts(example$covariate_data$Year_Cov, contrasts = FALSE)),
                         #Map = fit_adjust3_map,
                         parameters = fit_adjust3_params,
                         covariate_data = example$covariate_data,
                         run_model = FALSE)

# Conceptually a bit lost here, but that's alright. Intercept is now the average effect of the first year and not the average effect across all years. I guess that's okay? Thought process being it would be hard to estimate the other factors without a reference level, which is what happens with contrasts = FALSE.
fit_adjust3<- fit_model( "settings" = settings,
                         Lat_i = example$sampling_data[,'Lat'],
                         Lon_i = example$sampling_data[,'Lon'],
                         t_i = example$sampling_data[,'Year'],
                         b_i = example$sampling_data[,'Catch_KG'],
                         a_i = example$sampling_data[,'AreaSwept_km2'],
                         X1_formula = X1_formula,
                         X2_formula = X2_formula,
                         #X_contrasts = list(Year_Cov = contrasts(example$covariate_data$Year_Cov, contrasts = FALSE)),
                         #Map = fit_adjust3_map,
                         parameters = fit_adjust3_params,
                         covariate_data = example$covariate_data,
                         run_model = TRUE)









Just as an aside, the parameter estimates seem to aling with the `fit_example` ones for all of the parameters.
# Must add data-frames to global environment
covariate_data_full = fit_adjust2$effects$covariate_data_full
catchability_data_full = fit_adjust2$effects$catchability_data_full

# Plot 1st linear predictor
pred = Effect.fit_model( fit_adjust2,
                         focal.predictors = c("BOT_DEPTH"),
                         which_formula = "X1", 
                         xlevels = 100)


# Error -- this isn't the error I was getting (of course). Digging into the Effect.fit_model call and quickly see that the error from above IS going to be an issue as no parameter estimates are pulled. The saga continues...parameter must be at a bound...
fit_adjust2$parameter_estimates$opt$diagnostics

# Potential culprits -- logSigmaM? 
fit_adjust3_map<- fit_adjust2$tmb_list$Map
fit_adjust3_map['logSigmaM']<- factor(rep(NA, 3))

# Refitting
fit_adjust3 = fit_model( "settings" = settings,
                         Lat_i = example$sampling_data[,'Lat'],
                         Lon_i = example$sampling_data[,'Lon'],
                         t_i = example$sampling_data[,'Year'],
                         b_i = example$sampling_data[,'Catch_KG'],
                         a_i = example$sampling_data[,'AreaSwept_km2'],
                         X1_formula = X1_formula,
                         X2_formula = X2_formula,
                         X_contrasts = list(Year_Cov = contrasts(example$covariate_data$Year_Cov, contrasts = FALSE)),
                         covariate_data = example$covariate_data,
                         Map = fit_adjust3_map,
                         run_model = TRUE)









mod = fit_adjust2
focal.predictors = c("BOT_DEPTH")
which_formula = "X1"
xlevels = 100

# Function checks, no issues here...
# if (mod$data_list$n_c > 1) {
#   stop("`Effect.fit_model` is not currently designed for multivariate models")
# }
# if (!all(c("covariate_data_full", "catchability_data_full") %in%  ls(.GlobalEnv))) {
#   stop("Please load `covariate_data_full` and `catchability_data_full` into global memory")
# }
# if (!requireNamespace("effects")) {
#   stop("please install the effects package")
# }
# if (!("effects" %in% names(mod))) {
#   stop("`effects` slot not detected in input to `Effects.fit_model`. Please update model using later package version.")
# }

# Formula work...
if (which_formula == "X1") {
  formula_orig = mod$X1_formula
  parname = "gamma1_cp"
  mod$call = mod$effects$call_X1
} else if (which_formula == "X2") {
  formula_orig = mod$X2_formula
  parname = "gamma2_cp"
  mod$call = mod$effects$call_X2
} else if (which_formula == "Q1") {
  formula_orig = mod$Q1_formula
  parname = "lambda1_k"
  mod$call = mod$effects$call_Q1
} else if (which_formula == "Q2") {
  formula_orig = mod$Q2_formula
  parname = "lambda2_k"
  mod$call = mod$effects$call_Q2
} else {
  stop("Check `which_formula` input")
}
whichnum = which(names(mod$parameter_estimates$par) == parname)
mod$parhat = mod$parameter_estimates$par[whichnum]
mod$covhat = mod$parameter_estimates$SD$cov.fixed[whichnum, 
                                                  whichnum, drop = FALSE]
if (parname %in% names(mod$tmb_list$Obj$env$map)) {
  mod$parhat = mod$parhat[mod$tmb_list$Obj$env$map[[parname]]]
  mod$covhat = mod$covhat[mod$tmb_list$Obj$env$map[[parname]], 
                          mod$tmb_list$Obj$env$map[[parname]], drop = FALSE]
  mod$parhat = ifelse(is.na(mod$parhat), 0, mod$parhat)
  mod$covhat = ifelse(is.na(mod$covhat), 0, mod$covhat)
}
names(mod$parhat)[] = parname
rownames(mod$covhat) = colnames(mod$covhat) = names(mod$parhat)
formula_full = stats::update.formula(formula_orig, linear_predictor ~ 
                                       . + 0)
mod$coefficients = mod$parhat
mod$vcov = mod$covhat
mod$formula = formula_full
mod$family = stats::gaussian(link = "identity")
family.fit_model = function(x, ...) x$family
vcov.fit_model = function(x, ...) x$vcov
dummyfuns = list(variance = function(mu) mu, initialize = expression(mustart = y + 
                                                                       0.1), dev.resids = function(...) stats::poisson()$dev.res(...))
fam = mod$family
for (i in names(dummyfuns)) {
  if (is.null(fam[[i]])) 
    fam[[i]] = dummyfuns[[i]]
}
if (length(formals(fam$variance)) > 1) {
  warning("overriding variance function for effects: computed variances may be incorrect")
  fam$variance = dummyfuns$variance
}
args = list(call = mod$call, coefficients = mod$coefficients, 
            vcov = mod$vcov, family = fam, formula = formula_full)
effects::Effect.default(focal.predictors, mod, ..., sources = args)
sources <- if (missing(sources)) 
  effSources(mod)
else sources
formula <- if (is.null(sources$formula)) 
  insight::find_formula(mod)$conditional
else sources$formula
if (is.null(focal.predictors)) 
  return(formula)
cl <- if (is.null(sources$call)) {
  if (isS4(mod)) 
    mod@call
  else mod$call
}
else sources$call
cl$formula <- formula
type <- if (is.null(sources$type)) 
  "glm"
else sources$type
fam <- try(family(mod), silent = TRUE)
if (inherits(fam, "try-error")) 
  fam <- NULL
if (!is.null(sources$family)) {
  fam <- sources$family
}
if (!is.null(fam)) {
  fam$aic <- function(...) NULL
  if (!is.null(fam$variance)) {
    if (length(formals(fam$variance)) > 1) 
      stop("Effect plots are not implemented for families with more than\n             one parameter in the variance function (e.g., negitave binomials).")
  }
}
cl$family <- fam
coefficients <- if (is.null(sources$coefficients)) 
  effCoef(mod)
else sources$coefficients
vcov <- if (is.null(sources$vcov)) 
  as.matrix(vcov(mod, complete = TRUE))
else sources$vcov
zeta <- if (is.null(sources$zeta)) 
  NULL
else sources$zeta
cl$control <- switch(type, glm = glm.control(epsilon = Inf, 
                                             maxit = 1), polr = list(maxit = 1), multinom = c(maxit = 1))
cl$method <- sources$method
.m <- switch(type, glm = match(c("formula", "data", "family", 
                                 "contrasts", "subset", "control", "offset"), names(cl), 
                               0L), polr = match(c("formula", "data", "family", "contrasts", 
                                                   "subset", "control", "method"), names(cl), 0L), multinom = match(c("formula", 
                                                                                                                      "data", "family", "contrasts", "subset", "family", "maxit", 
                                                                                                                      "offset"), names(cl), 0L))
cl <- cl[c(1L, .m)]
cl[[1L]] <- as.name(type)
mod2 <- eval(cl)
mod2$coefficients <- coefficients
mod2$vcov <- vcov
if (!is.null(zeta)) 
  mod2$zeta <- zeta
if (type == "glm") {
  mod2$weights <- as.vector(with(mod2, prior.weights * 
                                   (family$mu.eta(linear.predictors)^2/family$variance(fitted.values))))
}
class(mod2) <- c("fakeeffmod", class(mod2))
Effect(focal.predictors, mod2, ...)

