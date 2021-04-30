function (formula, covariate_data, Year_i, spatial_list, contrasts.arg = NULL) 
{
  formula = X1_formula
  covariate_data = covariate_data[!is.na(covariate_data$SST_seasonal),]
  Year_i = t_i
  spatial_list = spatial_list
  contrasts.arg = X_contrasts
  if (!is.data.frame(covariate_data)) 
    stop("Please ensure that `covariate_data` is a data frame")
  if (!all(c("Lat", "Lon", "Year") %in% names(covariate_data))) {
    stop("`data` in `make_covariates(.)` must include columns `Lat`, `Lon`, and `Year`")
  }
  Year_Set = min(Year_i):max(Year_i)
  covariate_df = covariate_data[which(!is.na(covariate_data[, "Year"])), ]
  for (tI in seq_along(Year_Set)) {
    newrows = covariate_data[which(is.na(covariate_data[, "Year"])), ]
    newrows[, "Year"] = rep(Year_Set[tI], nrow(newrows))
    covariate_df = rbind(covariate_df, newrows)
    print(tI)
  }
  Model_matrix = model.matrix(update.formula(formula, ~. + 1), data = covariate_df, contrasts.arg = contrasts.arg)
  Columns_to_keep = which(attr(Model_matrix, "assign") != 0)
  coefficient_names = attr(Model_matrix, "dimnames")[[2]][Columns_to_keep]
  X = Model_matrix[, Columns_to_keep, drop = FALSE]
  dimnames(X) = list(NULL, coefficient_names)
  sample_i = data.frame(Year = Year_i, Lat = spatial_list$latlon_i[, "Lat"], Lon = spatial_list$latlon_i[, "Lon"])
  latlon_g = spatial_list$latlon_g
  X_gtp = array(NA, dim = c(nrow(latlon_g), length(Year_Set), ncol(X)), dimnames = list(NULL, Year_Set, colnames(X)))
  X_ip = array(NA, dim = c(nrow(sample_i), ncol(X)), dimnames = list(NULL, colnames(X)))
  for (tI in seq_along(Year_Set)) {
  #for (tI in 1:100) {
    print(tI)
    tmp_covariate_df = covariate_df[which(Year_Set[tI] ==  covariate_df[, "Year"]), , drop = FALSE]
    tmp_X = X[which(Year_Set[tI] == covariate_df[, "Year"]),  , drop = FALSE]
    if (nrow(tmp_covariate_df) == 0) {
      stop("Year ", Year_Set[tI], " not found in `covariate_data` please specify covariate values for all years")
    }
    Which = which(Year_Set[tI] == sample_i[, "Year"])
    if (length(Which) > 0) {
      NN = RANN::nn2(data = tmp_covariate_df[, c("Lat", "Lon")], query = sample_i[Which, c("Lat", "Lon")], k = 1)
      X_ip[Which, ] = tmp_X[NN$nn.idx[, 1], , drop = FALSE]
    }
    NN = RANN::nn2(data = tmp_covariate_df[, c("Lat", "Lon")], 
                   query = latlon_g[, c("Lat", "Lon")], k = 1)
    X_gtp[, tI, ] = tmp_X[NN$nn.idx[, 1], , drop = FALSE]
  }
  X_itp = aperm(X_ip %o% rep(1, length(Year_Set)), perm = c(1, 3, 2))
  if (any(is.na(X_itp))) 
    stop("Problem with `X_itp` in `make_covariates(.)")
  if (any(is.na(X_gtp))) 
    stop("Problem with `X_gtp` in `make_covariates(.)")
  if (any(apply(X_gtp, MARGIN = 2:3, FUN = sd) > 10 | apply(X_itp, 
                                                            MARGIN = 2:3, FUN = sd) > 10)) {
    warning("The package author recommends that you rescale covariates in `covariate_data` to have mean 0 and standard deviation 1.0")
  }
  Return = list(X_gtp = X_gtp, X_itp = X_itp, coefficient_names = coefficient_names)
  return(Return)
}
