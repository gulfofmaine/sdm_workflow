###### DEVELOPMENT VERSION OF VAST SHOULD BE USED
library(VAST)
library(splines)  # Used to include basis-splines
library(effects)  # Used to visualize covariate effects
library(mgcv)

#####
## For comparison with GAMs, leveraging "covariate example" wiki
####
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

# Define formula, here using a b-spline basis with the splines library. For the cubic regression spline, we want to set "degree = 3".
X1_formula = ~ bs( BOT_DEPTH, degree=3, intercept=FALSE)
X2_formula = ~ bs( BOT_DEPTH, degree=3, intercept=FALSE)

# If all covariates as "static" (not changing among years), then set Year = NA to cause values to be duplicated internally for all values of Year
example$covariate_data[,'Year'] = NA

# In the Wiki example, covariates are rescaled to have an SD >0.1 and <10 (for numerical stability) with the following line. The mgcv::gam function I believe automatically scales/centers any continuous covariates. So, rather than doing the rescaling as in the wiki example, going to scale/center instead for direct comparison between the two approaches.
# example$covariate_data[,'BOT_DEPTH'] = example$covariate_data[,'BOT_DEPTH'] / 100
example$covariate_data[,'BOT_DEPTH']<- scale(example$covariate_data[,'BOT_DEPTH'])

# Another thing we need to adjust is the default settings and the ObsModel default for "purpose = index2", which uses "ObsModel = c(2, 1)" and specifies the "Poisson" link delta model. I'm not entirely sure how to implement that in mgcv::gam. So, instead going to use "ObsModel = c(2, 0)", which is the "traditional" two-stage delta log normal GAM.
settings$ObsModel<- c(2, 0)

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

library(effects)  # Used to visualize covariate effects

# Must add data-frames to global environment
covariate_data_full = fit$effects$covariate_data_full
catchability_data_full = fit$effects$catchability_data_full

# Plot 1st linear predictor
pred = Effect.fit_model( fit,
                         focal.predictors = c("BOT_DEPTH"),
                         which_formula = "X1", 
                         xlevels = 100)
plot(pred)


# Now, trying to do the same thing with mgcv. To do this, relies heavily on Gavin Simpson's excellent tutorial here: https://fromthebottomoftheheap.net/2020/06/03/extrapolating-with-gams/. First, convert "response" to 0 and 1 and then also get the bottom depth covariate
gam_df<- data.frame("Response" = ifelse(example$sampling_data[,'Catch_KG'] > 0, 1, 0), "BOT_DEPTH" = example$covariate_data[,'BOT_DEPTH'])

# Specify and fit the model. According to Gavin's blog post and the R help file, a gam with s(variable, bs = "bs", m = c(3, 2)) is the "garden variety cubic B spline with second order penalty" -- so penalty placed on the second derivative of the spline, or its curvature (as opposed to m[2] = 1, which would penalize the first derivative and departures from a linear function). With splines::bs(), I wasn't sure if there was actually any penalty being applied and tried "m = c(3,0)" first. This gave something similar, but not identical. So, went ahead and tried "m = c(3, 2)" instead. 
gam_fit_2pen<- gam(Response ~ s(BOT_DEPTH, bs = "bs", m = c(3, 2)), data = gam_df, family = "binomial")
plot(gam_fit_2pen, main = "MGCV::GAM, m = c(3, 2)")
plot(pred, main = "VAST SPLINES::BS()")


######
## Model with a factor variable AND continuous variable
######
# Load packages
library(VAST)
library(splines)  # Used to include basis-splines
library(effects)  # Used to visualize covariate effects

# load data set
# see `?load_example` for list of stocks with example data
# that are installed automatically with `FishStatsUtils`.
example = load_example( data_set="covariate_example" )

# Make settings (turning off bias.correct to save time for example)
settings = make_settings( n_x=200,
                          Region=example$Region,
                          purpose="index2",
                          use_anisotropy=FALSE,
                          bias.correct=FALSE,
                          fine_scale=TRUE )

# Define formula.
X1_formula = ~ Year_Cov + bs( log(BOT_DEPTH), degree=2, intercept=FALSE)
X2_formula = ~ Year_Cov + bs( log(BOT_DEPTH), degree=2, intercept=FALSE)

# Create year covariate as a factor
example$covariate_data[,'Year_Cov']<- factor(example$covariate_data[,'Year'], levels = unique(example$covariate_data[,'Year']))
str(example$sampling_data)
# Get as a dataframe...
example$sampling_data<- data.frame(example$sampling_data)
example$sampling_data[,'Year_Cov']<- factor(example$sampling_data[,'Year'], levels = unique(example$sampling_data[,'Year']))

# Rescale covariates being used to have an SD >0.1 and <10 (for numerical stability)
example$covariate_data[,'BOT_DEPTH'] = example$covariate_data[,'BOT_DEPTH'] / 100

# Run model
fit = fit_model( "settings" = settings,
                 Lat_i = example$sampling_data[,'Lat'],
                 Lon_i = example$sampling_data[,'Lon'],
                 t_i = example$sampling_data[,'Year'],
                 b_i = example$sampling_data[,'Catch_KG'],
                 a_i = example$sampling_data[,'AreaSwept_km2'],
                 X1_formula = X1_formula,
                 X2_formula = X2_formula,
                 covariate_data = example$covariate_data )

# Nope, built in sequencing for years issue in make_covariates -- Year_Set = min(Year_i):max(Year_i). Just manually adjust years...
example$sampling_data[, 'Year']<- as.numeric(example$sampling_data[,'Year_Cov']) + 2000
example$sampling_data[, 'Year_Cov']<- factor(example$sampling_data[,'Year'], levels = unique(example$sampling_data[,'Year']))
example$covariate_data[, 'Year']<- as.numeric(example$covariate_data[,'Year_Cov']) + 2000
example$covariate_data[, 'Year_Cov']<- factor(example$covariate_data[,'Year'], levels = unique(example$covariate_data[,'Year']))

# Build model
fit = fit_model( "settings" = settings,
                 Lat_i = example$sampling_data[,'Lat'],
                 Lon_i = example$sampling_data[,'Lon'],
                 t_i = example$sampling_data[,'Year'],
                 b_i = example$sampling_data[,'Catch_KG'],
                 a_i = example$sampling_data[,'AreaSwept_km2'],
                 X1_formula = X1_formula,
                 X2_formula = X2_formula,
                 covariate_data = example$covariate_data,
                 run_model = FALSE)

# Make adjustments -- right now 10 gamma1_cp terms, should only be two, one for each of the parameters of the spline (not one for each years)
map_orig<- fit$tmb_list$Map
map_new<- map_orig
map_new$gamma1_cp<- factor(rep(c(1, 2), times = 5))
map_new$gamma2_cp<- factor(rep(c(1, 2), times = 5))

# Build model
fit = fit_model( "settings" = settings,
                 Lat_i = example$sampling_data[,'Lat'],
                 Lon_i = example$sampling_data[,'Lon'],
                 t_i = example$sampling_data[,'Year'],
                 b_i = example$sampling_data[,'Catch_KG'],
                 a_i = example$sampling_data[,'AreaSwept_km2'],
                 X1_formula = X1_formula,
                 X2_formula = X2_formula,
                 covariate_data = example$covariate_data,
                 Map = map_new,
                 run_model = FALSE)

# All looks good, fit
fit = fit_model( "settings" = settings,
                 Lat_i = example$sampling_data[,'Lat'],
                 Lon_i = example$sampling_data[,'Lon'],
                 t_i = example$sampling_data[,'Year'],
                 b_i = example$sampling_data[,'Catch_KG'],
                 a_i = example$sampling_data[,'AreaSwept_km2'],
                 X1_formula = X1_formula,
                 X2_formula = X2_formula,
                 covariate_data = example$covariate_data,
                 Map = map_new,
                 run_model = TRUE)

#####
## Recreating effects plots -- this doesn't work!
#####
# Must add data-frames to global environment
covariate_data_full = fit$effects$covariate_data_full
catchability_data_full = fit$effects$catchability_data_full

# Plot 1st linear predictor
pred = Effect.fit_model( fit,
                         focal.predictors = c("BOT_DEPTH"),
                         which_formula = "X1", 
                         xlevels = 100)
plot(pred)
# Same error I was getting. Digging into the function...

#####
## Function sleuthing...
#####

mod = fit #Fit example should work, fit gives the error
# Define formula.
X1_formula_full = X1_formula
X2_formula_full = X2_formula
# Must add data-frames to global environment (hope to fix in future)
covariate_data_full = mod$effects$covariate_data_full
catchability_data_full = mod$effects$catchability_data_full
focal.predictors = c("BOT_DEPTH")
which_formula = "X1"

# Error checks
if( mod$data_list$n_c>1 ){
  stop("`Effect.fit_model` is not currently designed for multivariate models")
}
if( !all(c("covariate_data_full","catchability_data_full") %in% ls(.GlobalEnv)) ){
  stop("Please load `covariate_data_full` and `catchability_data_full` into global memory")
}
if( !requireNamespace("effects") ){
  stop("please install the effects package")
}

# Check for effects slot in mod, if it isn't there, try to make it -- seems to happen if using older version of FishStatsUtils::fit_model. Not applicable here...
# if(is.null(mod$effects)){
#   mod$effects<- list()
#   if(!is.null(mod$catchability_data)) {
#     catchability_data_full = data.frame(mod$catchability_data, linear_predictor = 0)
#     Q1_formula_full = update.formula(mod$Q1_formula, linear_predictor ~ . + 0)
#     call_Q1 = lm(Q1_formula_full, data = catchability_data_full)$call
#     Q2_formula_full = update.formula(mod$Q2_formula, linear_predictor ~ . + 0)
#     call_Q2 = lm(Q2_formula_full, data = catchability_data_full)$call
#     mod$effects = c(Return$effects, list(call_Q1 = call_Q1, call_Q2 = call_Q2, catchability_data_full = catchability_data_full))
#   }
#   if(!is.null(mod$covariate_data)) {
#     covariate_data_full = data.frame(mod$covariate_data, linear_predictor = 0)
#     X1_formula_full = update.formula(mod$X1_formula, linear_predictor ~ . + 0)
#     call_X1 = lm(X1_formula_full, data = covariate_data_full)$call
#     X2_formula_full = update.formula(mod$X2_formula, linear_predictor ~ . + 0)
#     call_X2 = lm(X2_formula_full, data = covariate_data_full)$call
#     mod$effects = c(mod$effects, list(call_X1 = call_X1, call_X2 = call_X2, covariate_data_full = covariate_data_full))
#   }
# }

# Identify formula-specific stuff -- uses formula to grab the parameter name and corresponding formula
if( which_formula=="X1" ){
  formula_orig = mod$X1_formula
  parname = "gamma1_cp"
  mod$call = mod$effects$call_X1
} else if( which_formula=="X2" ){
  formula_orig = mod$X2_formula
  parname = "gamma2_cp"
  mod$call = mod$effects$call_X2
} else if( which_formula=="Q1" ){
  formula_orig = mod$Q1_formula
  parname = "lambda1_k"
  mod$call = mod$effects$call_Q1
} else if( which_formula=="Q2" ){
  formula_orig = mod$Q2_formula
  parname = "lambda2_k"
  mod$call = mod$effects$call_Q2
} else{
  stop("Check `which_formula` input")
}

# Extract parameters / covariance
# Get index of the gamma1_cp parameters...
whichnum = which(names(mod$parameter_estimates$par)==parname)

# Extract parameter estimate
mod$parhat = mod$parameter_estimates$par[whichnum]

# Extract covariance
mod$covhat = mod$parameter_estimates$SD$cov.fixed[whichnum,whichnum,drop=FALSE]

# Fill in values that are mapped off -- I'm not sure if this makes sense...
# if(parname %in% names(mod$tmb_list$Obj$env$map) ){
#   mod$parhat = mod$parhat[ mod$tmb_list$Obj$env$map[[parname]] ]
#   mod$covhat = mod$covhat[ mod$tmb_list$Obj$env$map[[parname]], mod$tmb_list$Obj$env$map[[parname]], drop=FALSE ]
#   mod$parhat = ifelse( is.na(mod$parhat), 0, mod$parhat)
#   mod$covhat = ifelse( is.na(mod$covhat), 0, mod$covhat)
# }
# add names
names(mod$parhat)[] = parname
rownames(mod$covhat) = colnames(mod$covhat) = names(mod$parhat)

# Augment stuff
formula_full = update.formula(formula_orig, linear_predictor~.+0)
mod$coefficients = mod$parhat
mod$vcov = mod$covhat
mod$formula = formula_full
mod$family = gaussian(link = "identity")

# Functions for package
family.fit_model = function(x,...) x$family
vcov.fit_model = function(x,...) x$vcov

# dummy functions to make Effect.default work
dummyfuns = list(variance = function(mu) mu,
                 initialize = expression(mustart = y + 0.1),
                 dev.resids = function(...) poisson()$dev.res(...) )

# Replace family (for reasons I don't really understand)
fam = mod$family
for( i in names(dummyfuns) ){
  if( is.null(fam[[i]]) ) fam[[i]] = dummyfuns[[i]]
}

# allow calculation of effects ...
if (length(formals(fam$variance))>1) {
  warning("overriding variance function for effects: computed variances may be incorrect")
  fam$variance = dummyfuns$variance
}

# Bundle arguments
args = list(call = mod$call,
            coefficients = mod$coefficients,
            vcov = mod$vcov,
            family = fam,
            formula = formula_full)

# Do call
effects::Effect.default(focal.predictors,
                        mod,
                        sources = args)

# Error...going into effects::Effect.default...
focal.predictors = focal.predictors
mod = mod
sources = args

# Creating objects given null or missing inputs. 
sources <- if (missing(sources))
  effSources(mod) else sources
formula <- if (is.null(sources$formula))
  insight::find_formula(mod)$conditional else sources$formula
if(is.null(focal.predictors))
  return(formula)
cl<- if(is.null(sources$call)) {
  if(isS4(mod))
    mod@call else mod$call
} else sources$call
cl$formula <- formula

type <- if (is.null(sources$type)) 
  "glm" else sources$type
fam <- try(family(mod), silent = TRUE)
if(inherits(fam, "try-error")) 
  fam <- NULL
if(!is.null(sources$family)) {
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
coefficients <- if(is.null(sources$coefficients)) 
  effCoef(mod) else sources$coefficients
vcov <- if (is.null(sources$vcov)) 
  as.matrix(vcov(mod, complete = TRUE)) else sources$vcov
zeta <- if (is.null(sources$zeta)) 
  NULL else sources$zeta
cl$control <- switch(type, glm = glm.control(epsilon = Inf, maxit = 1), polr = list(maxit = 1), multinom = c(maxit = 1))
cl$method <- sources$method
.m <- switch(type, glm = match(c("formula", "data", "family", "contrasts", "subset", "control", "offset"), names(cl), 0L), polr = match(c("formula", "data", "family", "contrasts", "subset", "control", "method"), names(cl), 0L), multinom = match(c("formula", "data", "family", "contrasts", "subset", "family", "maxit", "offset"), names(cl), 0L))
cl <- cl[c(1L, .m)]
cl[[1L]] <- as.name(type)
mod2 <- eval(cl)
mod2$coefficients <- coefficients
mod2$vcov <- vcov
if (!is.null(zeta)) 
  mod2$zeta <- zeta
if (type == "glm") {
  mod2$weights <- as.vector(with(mod2, prior.weights * (family$mu.eta(linear.predictors)^2/family$variance(fitted.values))))
}
class(mod2) <- c("fakeeffmod", class(mod2))

# Weird issue -- so when I do this with the Year_Cov + bs(DEPTH), I get a weird output from mod2 that again shows: glm(formula = linear_predictor ~ Year_Cov + bs(log(BOT_DEPTH), degree = 2, intercept = FALSE) - 1 and again the 10 gamma1 coefficients! In contrast, with the example that works (just bs(DEPTH), I get the correct 3 coefficients...where is that happening above??). I can solve that if I comment out the mapping parameters bits 259-265. That still doesn't solve the error though...
Effect(focal.predictors, mod2)

# Digging into Effect.lm
focal.predictors = focal.predictors
mod = mod 
xlevels = 100
vcov. = vcov
se = TRUE
residuals = FALSE
quantiles = seq(0.2, 0.8, by = 0.2)
x.var = NULL

partial.residuals <- residuals
if (missing(transformation)) 
  transformation <- list(link = family(mod)$linkfun, inverse = family(mod)$linkinv)
if (missing(fixed.predictors)) 
  fixed.predictors <- NULL
fixed.predictors <- effects:::applyDefaults(fixed.predictors, list(given.values = NULL, typical = mean, apply.typical.to.factors = FALSE, offset = mean), arg = "fixed.predictors")
given.values <- fixed.predictors$given.values
typical <- fixed.predictors$typical
offset <- fixed.predictors$offset
apply.typical.to.factors <- fixed.predictors$apply.typical.to.factors
confint <- effects:::applyDefaults(se, list(compute = TRUE, level = 0.95, type = "pointwise"), onFALSE = list(compute = FALSE, level = 0.95, type = "pointwise"), arg = "se")
se <- confint$compute
confidence.level <- confint$level
confidence.type <- match.arg(confint$type, c("pointwise", "Scheffe", "scheffe"))
default.levels <- NULL

data <- if(partial.residuals) {
  all.vars <- all.vars(formula(mod))
  expand.model.frame(mod, all.vars)[, all.vars]
} else NULL
if (!is.null(given.values) && !all(which <- names(given.values) %in%  names(coef(mod)))) 
  stop("given.values (", names(given.values[!which]), ") not in the model")
off<- if (is.numeric(offset) && length(offset) == 1) 
  offset else if (is.function(offset)) {
  mod.off <- model.offset(model.frame(mod))
  if (is.null(mod.off)) 
    0 else offset(mod.off)
} else stop("offset must be a function or a number")
formula.rhs <- formula(mod)[[3]]
if (!missing(x.var)) {
  if (!is.numeric(x.var)) {
    x.var.name <- x.var
    x.var <- which(x.var == focal.predictors)
  }
  if (length(x.var) == 0) 
    stop("'", x.var.name, "' is not among the focal predictors")
  if (length(x.var) > 1) 
    stop("x.var argument must be of length 1")
}
model.components <- Analyze.model(focal.predictors, mod, 
                                  xlevels, default.levels, formula.rhs, partial.residuals = partial.residuals, 
                                  quantiles = quantiles, x.var = x.var, data = data, typical = typical)
excluded.predictors <- model.components$excluded.predictors
predict.data <- model.components$predict.data
predict.data.all.rounded <- predict.data.all <- if (partial.residuals) 
  na.omit(data[, all.vars(formula(mod))])
else NULL
factor.levels <- model.components$factor.levels
factor.cols <- model.components$factor.cols
n.focal <- model.components$n.focal
x <- model.components$x
X.mod <- model.components$X.mod
cnames <- model.components$cnames
X <- model.components$X
x.var <- model.components$x.var
formula.rhs <- formula(mod)[c(1, 3)]
Terms <- delete.response(terms(mod))
mf <- model.frame(Terms, predict.data, xlev = factor.levels, 
                  na.action = NULL)
mod.matrix <- model.matrix(formula.rhs, data = mf, contrasts.arg = mod$contrasts)
if (is.null(x.var)) 
  partial.residuals <- FALSE
factors <- sapply(predict.data, is.factor)
if (partial.residuals) {
  for (predictor in focal.predictors[-x.var]) {
    if (!factors[predictor]) {
      values <- unique(predict.data[, predictor])
      predict.data.all.rounded[, predictor] <- values[apply(outer(predict.data.all[, 
                                                                                   predictor], values, function(x, y) (x - y)^2), 
                                                            1, which.min)]
    }
  }
}
mod.matrix.all <- model.matrix(mod)
wts <- weights(mod)
if (is.null(wts)) 
  wts <- rep(1, length(residuals(mod)))
mod.matrix <- Fixup.model.matrix(mod, mod.matrix, mod.matrix.all, 
                                 X.mod, factor.cols, cnames, focal.predictors, excluded.predictors, 
                                 typical, given.values, apply.typical.to.factors)
null.basis <- estimability::nonest.basis(mod)
is.estimable <- estimability::is.estble(mod.matrix, null.basis)
scoef <- ifelse(is.na(mod$coefficients), 0L, mod$coefficients)
effect <- off + mod.matrix %*% scoef
effect[!is.estimable] <- NA
if (partial.residuals) {
  res <- na.omit(residuals(mod, type = "working"))
  fitted <- na.omit(if (inherits(mod, "glm")) 
    predict(mod, type = "link")
    else predict(mod))
  partial.residuals.range <- range(fitted + res)
}
else {
  res <- partial.residuals.range <- NULL
}
result <- list(term = paste(focal.predictors, collapse = "*"), 
               formula = formula(mod), response = response.name(mod), 
               variables = x, fit = effect, x = predict.data[, 1:n.focal, 
                                                             drop = FALSE], x.all = predict.data.all.rounded[, 
                                                                                                             focal.predictors, drop = FALSE], model.matrix = mod.matrix, 
               data = X, discrepancy = 0, offset = off, residuals = res, 
               partial.residuals.range = partial.residuals.range, x.var = x.var)
if (se) {
  if (any(family(mod)$family == c("binomial", "poisson"))) {
    z <- if (confidence.type == "pointwise") {
      qnorm(1 - (1 - confidence.level)/2)
    }
    else {
      p <- length(na.omit(coef(mod)))
      scheffe(confidence.level, p)
    }
  }
  else {
    z <- if (confidence.type == "pointwise") {
      qt(1 - (1 - confidence.level)/2, df = mod$df.residual)
    }
    else {
      p <- length(na.omit(coef(mod)))
      scheffe(confidence.level, p, mod$df.residual)
    }
  }
  V <- vcov.(mod, complete = FALSE)
  mmat <- mod.matrix[, !is.na(mod$coefficients)]
  eff.vcov <- mmat %*% V %*% t(mmat)
  rownames(eff.vcov) <- colnames(eff.vcov) <- NULL
  var <- diag(eff.vcov)
  result$vcov <- eff.vcov
  result$se <- sqrt(var)
  result$se[!is.estimable] <- NA
  result$lower <- effect - z * result$se
  result$upper <- effect + z * result$se
  result$confidence.level <- confidence.level
}
if (is.null(transformation$link) && is.null(transformation$inverse)) {
  transformation$link <- I
  transformation$inverse <- I
}
result$transformation <- transformation
result$family <- family(mod)$family
result$link <- family(mod)
class(result) <- "eff"
result
