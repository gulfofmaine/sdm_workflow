######
## Exploring problematic parameter and singularity
######

# First, general advice is to use getSD = FALSE and newtonsteps = 0, then fit the model. After doing that, we can run TMBhelper::check_estimatability(vast_model$tmb_list$Obj) to see if all parameters are estimatable. In our case, clearly getting some issue that is arising with the suvey catchability coefficient (lambda). We can see that through:
names(vast_fit$tmb_list$Lower)
TMBhelper::check_estimability(vast_fit$tmb_list$Obj)
# What's going on internally? Checking the make_data call
# Annual
tar_load(vast_settings)
settings = vast_settings
tar_load(vast_extrap_grid)
extrap_grid = vast_extrap_grid
tar_load(vast_sample_data)
sample_data = vast_sample_data
tar_load(vast_covariate_data)
covariate_data = vast_covariate_data
X1_formula = hab_formula
X2_formula = hab_formula
hab_env_coeffs_n = hab_env_coeffs_n
tar_load(vast_catchability_data)
catchability_data = vast_catchability_data
catch_formula<- ~ Survey
Q1_formula = catch_formula
Q2_formula = catch_formula
tar_load(vast_coveff)
Xconfig_list = vast_coveff
working_dir = paste0(getwd(), "/")
input_grid = as.matrix(extrap_grid, ncol = 3)
"Lat_i" = sample_data[, 'Lat']
"Lon_i" = sample_data[, 'Lon']
"t_i" = sample_data[, 'Year']
"c_i" = rep(0, nrow(sample_data))
"b_i" = sample_data[, 'Biomass']
"a_i" = sample_data[, 'Swept']
"v_i" = rep(0, length(b_i))
"c_iz" = rep(0, length(b_i))
"PredTF_i" = sample_data[, 'Pred_TF']
"X1config_cp" = Xconfig_list[['X1config_cp']]
"X2config_cp" = Xconfig_list[['X2config_cp']]
"covariate_data" = covariate_data
"X1_formula" = X1_formula
"X2_formula" = X2_formula
"catchability_data" = catchability_data
"Q1_formula" = Q1_formula
"Q2_formula" = Q2_formula
"Q1config_k" = Xconfig_list[['Q1config_k']]
"Q2config_k" = Xconfig_list[['Q2config_k']]

# Capture extra arguments to function
extra_args = list("X_contrasts" = list(Year_Cov = contrasts(vast_covariate_data$Year_Cov, contrasts = FALSE)), "Xconfig_list" = vast_coveff, "input_grid" = as.matrix(extrap_grid, ncol = 3))
# Backwards-compatible way to capture previous format to input extra arguments for each function via specific input-lists
extra_args = c( extra_args,
                extra_args$extrapolation_args,
                extra_args$spatial_args,
                extra_args$optimize_args,
                extra_args$model_args )

# Assemble inputs
data_frame = data.frame( "Lat_i"=Lat_i, "Lon_i"=Lon_i, "a_i"=a_i, "v_i"=v_i, "b_i"=b_i, "t_i"=t_i, "c_iz"=c_iz )
# Decide which years to plot
year_labels = seq( min(t_i), max(t_i) )
years_to_plot = which( year_labels %in% t_i )

# Save record
message("\n### Writing output from `fit_model` in directory: ", working_dir)
dir.create(working_dir, showWarnings=FALSE, recursive=TRUE)
#save( settings, file=file.path(working_dir,"Record.RData"))
capture.output( settings, file=file.path(working_dir,"settings.txt"))

# Build extrapolation grid
message("\n### Making extrapolation-grid")
extrapolation_args_default = list(Region = settings$Region,
                                  strata.limits = settings$strata.limits,
                                  zone = settings$zone,
                                  max_cells = settings$max_cells,
                                  DirPath = working_dir)
extrapolation_args_input = combine_lists( input = extra_args,
                                          default = extrapolation_args_default,
                                          args_to_use = formalArgs(make_extrapolation_info) )
extrapolation_list = do.call( what=make_extrapolation_info, args=extrapolation_args_input )

# Build information regarding spatial location and correlation
message("\n### Making spatial information")
spatial_args_default = list( grid_size_km = settings$grid_size_km,
                             n_x = settings$n_x,
                             Method = settings$Method,
                             Lon_i = Lon_i,
                             Lat_i = Lat_i,
                             Extrapolation_List = extrapolation_list,
                             DirPath = working_dir,
                             Save_Results = TRUE,
                             fine_scale = settings$fine_scale,
                             knot_method = settings$knot_method)
spatial_args_input = combine_lists( input=extra_args, default=spatial_args_default, args_to_use=c(formalArgs(make_spatial_info),formalArgs(INLA::inla.mesh.create)) )
spatial_list = do.call( what=make_spatial_info, args=spatial_args_input )

# Build data
# Do *not* restrict inputs to formalArgs(make_data) because other potential inputs are still parsed by make_data for backwards compatibility
message("\n### Making data object") # VAST::
if(missing(covariate_data)) covariate_data = NULL
if(missing(catchability_data)) catchability_data = NULL
"Version" = settings$Version 
"FieldConfig" = settings$FieldConfig
"OverdispersionConfig" = settings$OverdispersionConfig
"RhoConfig" = settings$RhoConfig
"VamConfig" = settings$VamConfig
"ObsModel" = settings$ObsModel
"c_iz" = c_iz
"b_i" = b_i
"a_i" = a_i
"v_i" = v_i
"s_i" = spatial_list$knot_i-1
"t_i" = t_i
"spatial_list" = spatial_list
"Options" = settings$Options
"Aniso" = settings$use_anisotropy
"X1config_cp" = X1config_cp
"X2config_cp" = X2config_cp
"covariate_data" = covariate_data
"X1_formula" = X1_formula
"X2_formula" = X2_formula
"Q1config_k" = Q1config_k
"Q2config_k" = Q2config_k
"catchability_data" = catchability_data
"Q1_formula" = Q1_formula
"Q2_formula" = Q2_formula

alternate_inputs = list()
if (missing(spatial_list)) {
  warning("Consider changing use of `make_data` to include `spatial_list` as input")
  a_xl = a_gl = alternate_inputs[["a_xl"]]
  MeshList = alternate_inputs[["MeshList"]]
  GridList = alternate_inputs[["GridList"]]
  Method = alternate_inputs[["Method"]]
  s_i = alternate_inputs[["s_i"]]
} else {
  MeshList = spatial_list[["MeshList"]]
  GridList = spatial_list[["GridList"]]
  Method = spatial_list[["Method"]]
  a_xl = a_gl = spatial_list[["a_gl"]]
  s_i = spatial_list[["knot_i"]] - 1
}
if (missing(t_i) & !is.null(alternate_inputs[["t_iz"]])) {
  message("Detecting deprecated input `t_iz` and coercing this to expected input `t_i`")
  t_i = alternate_inputs$t_iz[, 1]
}
if ("X_xj" %in% names(alternate_inputs)) 
  stop("`X_xj` is fully deprecated; please check inputs")
Options2use = c(SD_site_density = FALSE, SD_site_logdensity = FALSE, 
                Calculate_Range = FALSE, SD_observation_density = FALSE, 
                Calculate_effective_area = FALSE, Calculate_Cov_SE = FALSE, 
                Calculate_Synchrony = FALSE, Calculate_Coherence = FALSE, 
                Calculate_proportion = FALSE, normalize_GMRF_in_CPP = TRUE, 
                Calculate_Fratio = FALSE, Estimate_B0 = FALSE, Project_factors = FALSE, 
                treat_nonencounter_as_zero = FALSE, simulate_random_effects = TRUE, 
                observation_error_as_CV = TRUE, report_additional_variables = FALSE, 
                zerosum_penalty = 0)
for (i in seq_along(Options)) {
  if (tolower(names(Options)[i]) %in% tolower(names(Options2use))) {
    Options2use[[match(tolower(names(Options)[i]), tolower(names(Options2use)))]] = Options[[i]]
  }
}  
if (is.vector(FieldConfig) && length(FieldConfig) == 4) {
  FieldConfig = rbind(matrix(FieldConfig, ncol = 2, dimnames = list(c("Omega", "Epsilon"), c("Component_1", "Component_2"))), Beta = c("IID", "IID"))
} 
if (is.matrix(FieldConfig) & all(dim(FieldConfig) == c(3, 2))) {
  FieldConfig = rbind(FieldConfig, Epsilon_time = c("Identity", "Identity"))
} 
if (!is.matrix(FieldConfig) || !all(dim(FieldConfig) == c(4, 2))) {
  stop("`FieldConfig` has the wrong dimensions in `make_data`")
}  
dimnames(FieldConfig) = list(c("Omega", "Epsilon", "Beta", "Epsilon_time"), c("Component_1", "Component_2"))
tprime_i = t_i - min(t_i, na.rm = TRUE)
if (Options2use[12] == 1) {
  tprime_i = tprime_i + 1
  F_ct = cbind(0, F_ct)
}
if (!is.matrix(c_iz)) 
  c_iz = matrix(c_iz, ncol = 1)
n_t = max(tprime_i, na.rm = TRUE) + 1
n_c = max(c_iz, na.rm = TRUE) + 1
n_e = max(e_i) + 1
n_v = length(unique(v_i))
n_i = length(b_i)
n_x = nrow(a_gl)
n_l = ncol(a_gl)
n_g = ifelse(is.null(spatial_list), 1, spatial_list$n_g)  
if (!is.matrix(ObsModel_ez)) 
  ObsModel_ez = matrix(ObsModel_ez, ncol = 2, nrow = n_e, byrow = TRUE)
if (Options2use["treat_nonencounter_as_zero"] == TRUE) {
  Index = list(factor(c_iz[, 1], levels = 0:max(c_iz[, 1])), factor(tprime_i, levels = 0:max(tprime_i)))
  Num_ct = tapply(b_i, INDEX = Index, FUN = function(vec) {
    sum(vec > 0, na.rm = TRUE)
  })
  Num_ct = ifelse(is.na(Num_ct), 0, Num_ct)
  b_i = ifelse(Num_ct[cbind(as.numeric(Index[[1]]), as.numeric(Index[[2]]))] ==  0, NA, b_i)
}
X_xtp = alternate_inputs[["X_xtp"]]
X_gctp = alternate_inputs[["X_gctp"]]
X1_gctp = alternate_inputs[["X1_gctp"]]
X2_gctp = alternate_inputs[["X2_gctp"]]
X_itp = alternate_inputs[["X_itp"]]
X1_itp = alternate_inputs[["X1_itp"]]
X2_itp = alternate_inputs[["X2_itp"]]
if (missing(X1_formula) & missing(X2_formula) & "formula" %in% names(alternate_inputs)) {
  X1_formula = X2_formula = alternate_inputs[["formula"]]
}
if (is.null(X1config_cp) & is.null(X2config_cp) & "Xconfig_zcp" %in% names(alternate_inputs)) {
  X1config_cp = array(alternate_inputs[["Xconfig_zcp"]][1, , ], dim = dim(alternate_inputs[["Xconfig_zcp"]])[2:3])
  X2config_cp = array(alternate_inputs[["Xconfig_zcp"]][2, , ], dim = dim(alternate_inputs[["Xconfig_zcp"]])[2:3])
}
if (is.null(X_itp) & "X_ip" %in% names(alternate_inputs)) {
  X_itp = aperm(outer(alternate_inputs[["X_ip"]], rep(1, n_t)), c(1, 3, 2))
}
if (is.null(X_gctp) & "X_gtp" %in% names(alternate_inputs)) {
  X_gctp = aperm(outer(alternate_inputs[["X_gtp"]], rep(1, n_c)), c(1, 4, 2, 3))
}
Covariates_created = FALSE
if (Covariates_created == FALSE) {
  if (!is.null(X1_gctp) & !is.null(X1_itp) & !is.null(X2_gctp) &  !is.null(X2_itp)) {
    Covariates_created = TRUE
    if (!is.array(X1_gctp) || !(all(dim(X1_gctp)[1:3] == c(n_g, n_c, n_t)))) {
      stop("`X1_gctp` has wrong dimensions")
    }
    if (!is.array(X1_itp) || !(all(dim(X1_itp)[1:2] == c(n_i, n_t)))) {
      stop("`X1_itp` has wrong dimensions")
    }
    if (!is.array(X2_gctp) || !(all(dim(X2_gctp)[1:3] == c(n_g, n_c, n_t)))) {
      stop("`X2_gctp` has wrong dimensions")
    }
    if (!is.array(X2_itp) || !(all(dim(X2_itp)[1:2] == c(n_i, n_t)))) {
      stop("`X2_itp` has wrong dimensions")
    }
  }
}
if (Covariates_created == FALSE) {
  if (!is.null(X_gctp) & !is.null(X_itp)) {
    Covariates_created = TRUE
    X1_gctp = X2_gctp = X_gctp
    X1_itp = X2_itp = X_itp
    if (!is.array(X_gctp) || !(all(dim(X_gctp)[1:3] == c(n_g, n_c, n_t)))) {
      stop("`X_gctp` has wrong dimensions")
    }
    if (!is.array(X_itp) || !(all(dim(X_itp)[1:2] == c(n_i, n_t)))) {
      stop("`X_itp` has wrong dimensions")
    }
  }
}
if (Covariates_created == FALSE) {
  if (!is.null(covariate_data)) {
    if ("formula" %in% names(alternate_inputs)) {
      Covariates_created = TRUE
      warning("Using input `formula` to generate covariates. This interface is soft-deprecated but still available for backwards compatibility; please switch to using `X1_formula` and `X2_formula`")
      covariate_list = FishStatsUtils::make_covariates(formula = alternate_inputs[["formula"]], 
                                                       covariate_data = covariate_data, Year_i = t_i, 
                                                       spatial_list = spatial_list, contrasts.arg = X_contrasts)
      X1_gtp = X2_gtp = covariate_list$X_gtp
      X1_itp = X2_itp = covariate_list$X_itp
      X1_gctp = X2_gctp = aperm(outer(X1_gtp, rep(1, n_c)), c(1, 4, 2, 3))
    }
  }
}
if (Covariates_created == FALSE) {
  if (!is.null(covariate_data)) {
    Covariates_created = TRUE
    if (FishStatsUtils::convert_version_name(Version) <= 
        FishStatsUtils::convert_version_name("VAST_v9_4_0")) {
      stop("To use separate formula interface for linear predictors, please use version >= CPP 10.0.0")
    }
    covariate_list = FishStatsUtils::make_covariates(formula = X1_formula, 
                                                     covariate_data = covariate_data, Year_i = t_i, 
                                                     spatial_list = spatial_list, contrasts.arg = X_contrasts)
    X1_gtp = covariate_list$X_gtp
    X1_itp = covariate_list$X_itp
    X1_gctp = aperm(outer(X1_gtp, rep(1, n_c)), c(1, 4, 2, 3))
    covariate_list = FishStatsUtils::make_covariates(formula = X2_formula, 
                                                     covariate_data = covariate_data, Year_i = t_i, 
                                                     spatial_list = spatial_list, contrasts.arg = X_contrasts)
    X2_gtp = covariate_list$X_gtp
    X2_itp = covariate_list$X_itp
    X2_gctp = aperm(outer(X2_gtp, rep(1, n_c)), c(1, 4, 2, 3))
  }
}
if (Covariates_created == FALSE) {
  if (is.null(covariate_data)) {
    Covariates_created = TRUE
    X1_gctp = X2_gctp = array(0, dim = c(n_g, n_c, n_t, 0))
    X1_itp = X2_itp = array(0, dim = c(n_i, n_t, 0))
  }
}
create_Xconfig = function(Xconfig_cp, n_c, n_p) {
  if (is.null(Xconfig_cp)) {
    Xconfig_cp = array(1, dim = c(n_c, n_p))
  } else {
    if (!is.array(Xconfig_cp) || !(all(dim(Xconfig_cp) == c(n_c, n_p)))) {
      stop("`Xconfig_cp` has wrong dimensions")
    }
    if (!all(Xconfig_cp %in% c(-1, 0, 1, 2, 3))) {
      stop("`Xconfig_cp` has some wrong element(s)")
    }
    if (any(Xconfig_cp %in% -1)) {
      warning("Using `Xconfig_cp[] = -1` is unconventional and warrents caution")
    }
  }
  return(Xconfig_cp)
}
X1config_cp = create_Xconfig(Xconfig_cp = X1config_cp, n_c = n_c, n_p = dim(X1_gctp)[4])
X2config_cp = create_Xconfig(Xconfig_cp = X2config_cp, n_c = n_c, n_p = dim(X2_gctp)[4])

Catchability_created = FALSE
if (Catchability_created == FALSE) {
  if (!is.null(catchability_data)) {
    if (nrow(catchability_data) != n_i) 
      stop("`catchability_data` has the wrong number of rows; please supply one row for each observation `i`")
    Catchability_created = TRUE
    Model_matrix1 = stats::model.matrix(stats::update.formula(Q1_formula, ~. + 1), data = catchability_data)
    Columns_to_keep = which(attr(Model_matrix1, "assign") != 0)
    coefficient_names_Q1 = attr(Model_matrix1, "dimnames")[[2]][Columns_to_keep]
    Q1_ik = Model_matrix1[, Columns_to_keep, drop = FALSE]
    Model_matrix2 = stats::model.matrix(stats::update.formula(Q2_formula, ~. + 1), data = catchability_data)
    Columns_to_keep = which(attr(Model_matrix2, "assign") != 0)
    coefficient_names_Q2 = attr(Model_matrix2, "dimnames")[[2]][Columns_to_keep]
    Q2_ik = Model_matrix2[, Columns_to_keep, drop = FALSE]
  }
}
if (Catchability_created == FALSE) {
  Q1_ik = Q2_ik = matrix(0, nrow = n_i, ncol = 0)
}
create_Qconfig = function(Qconfig_k, n_k) {
  if (is.null(Qconfig_k)) {
    Qconfig_k = array(1, dim = c(n_k))
  } else {
    if (!is.vector(Qconfig_k) || length(Qconfig_k) != n_k) {
      stop("`Qconfig_k` has wrong length")
    }
    if (!all(Qconfig_k %in% c(0, 1, 2, 3))) {
      stop("`Qconfig_k` has some wrong element(s)")
    }
  }
  return(Qconfig_k)
}
Q1config_k = create_Qconfig(Qconfig_k = Q1config_k, n_k = ncol(Q1_ik))
Q2config_k = create_Qconfig(Qconfig_k = Q2config_k, n_k = ncol(Q2_ik))

str(Q1_ik)
str(Q1config_k)

# Adding a column for the DUMMY survey, which isn't what we want. Just make this 
