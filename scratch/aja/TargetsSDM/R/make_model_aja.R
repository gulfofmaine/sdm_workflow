make_model_aja<- function (TmbData, Version, RhoConfig = c(Beta1 = 0, Beta2 = 0, Epsilon1 = 0, Epsilon2 = 0), Method = "Mesh", Npool = 0, ConvergeTol = 1, Use_REML = FALSE, loc_x = NULL, Parameters = "generate", Random = "generate", Map = "generate", DiagnosticDir = NULL, TmbDir = system.file("executables", package = "VAST"), RunDir = getwd(), CompileDir = TmbDir, build_model = TRUE){
  
  if(FALSE){
    TmbData = model_args_input$TmbData
    Version = model_args_input$Version
    RhoConfig = model_args_input$RhoConfig
    Method = model_args_input$Method
    Npool = 0
    ConvergeTol = 1
    Use_REML = FALSE
    loc_x = model_args_input$loc_x
    Parameters = model_args_input$Parameters
    Random = model_args_input$Random
    Map = model_args_input$Map
    DiagnosticDir = NULL
    TmbDir = here::here("")
    RunDir = here::here("")
    CompileDir = TmbDir
    build_model = TRUE
  }
  
  if (all(c("Options", "Options_vec") %in% names(TmbData))) {
    Options_vec = TmbData$Options_vec
    Options = TmbData$Options
  }
  if ("Options_list" %in% names(TmbData)) {
    Options_vec = TmbData$Options_list$Options_vec
    Options = TmbData$Options_list$Options
  }
  if (!("n_e" %in% names(TmbData))) {
    TmbData[["n_e"]] = TmbData$n_c
  }
  if (!("ObsModel_ez" %in% names(TmbData))) {
    TmbData[["ObsModel_ez"]] = rep(1, TmbData[["n_e"]]) %o% 
      TmbData$ObsModel
  }
  if (!("c_iz" %in% names(TmbData))) {
    TmbData[["c_iz"]] = matrix(TmbData$c_i, ncol = 1)
  }
  if (!("e_i" %in% names(TmbData))) {
    TmbData[["e_i"]] = TmbData$c_iz[, 1]
  }
  if (!("t_iz" %in% names(TmbData))) {
    TmbData[["t_iz"]] = matrix(TmbData$t_i, ncol = 1)
  }
  capture.output(packageDescription("VAST"), file = paste0(RunDir, "/packageDescription.txt"))
  capture.output(packageDescription("FishStatsUtils"), file = paste0(RunDir, "/packageDescription.txt"), append = TRUE)
  if (length(Parameters) == 1 && Parameters == "generate") 
    Parameters = make_parameters(Version = Version, DataList = TmbData, RhoConfig = RhoConfig)
  if (length(Map) == 1 && Map == "generate") { RhoConfig = RhoConfig, Npool = Npool)
  } else {
    warning("Please carefully check starting values for all parameters to ensure that mapping off parameters will work as expected.")
  }
  if (length(Random) == 1 && Random == "generate") {
    Random = c("Epsiloninput1_sct", "Omegainput1_sc", "Epsiloninput1_sft", 
               "Omegainput1_sf", "eta1_vf", "Xiinput1_scp", "Phiinput1_sk", 
               "Epsiloninput1_sff", "Epsiloninput2_sct", "Omegainput2_sc", 
               "Epsiloninput2_sft", "Omegainput2_sf", "eta2_vf", 
               "Xiinput2_scp", "Phiinput2_sk", "Epsiloninput2_sff", 
               "delta_i")
    if (RhoConfig[["Beta1"]] %in% c(1, 2, 4)) 
      Random = c(Random, "beta1_ct", "beta1_ft")
    if (RhoConfig[["Beta2"]] %in% c(1, 2, 4)) 
      Random = c(Random, "beta2_ct", "beta2_ft")
    if (Use_REML == TRUE) {
      Random = union(Random, c("beta1_ct", "beta1_ft", 
                               "gamma1_j", "gamma1_tp", "gamma1_ctp", "lambda1_k", 
                               "gamma1_cp", "beta2_ct", "beta2_ft", "gamma2_j", 
                               "gamma2_tp", "gamma2_ctp", "lambda2_k", "gamma2_cp"))
    }
    Random = Random[which(Random %in% names(Parameters))]
    if (length(Random) == 0) 
      Random = NULL
  }
  if (build_model == FALSE) {
    Return = list(Map = Map, Data = TmbData, Parameters = Parameters, Random = Random)
    return(Return)
  }
  file.copy(from = paste0(TmbDir, "/", Version, ".cpp"), to = paste0(CompileDir, "/", Version, ".cpp"), overwrite = FALSE)
  origwd = getwd()
  on.exit(setwd(origwd), add = TRUE)
  setwd(here::here())
  TMB::compile(paste0(Version, ".cpp"), flags = "-Wno-ignored-attributes -O2 -mfpmath=sse -msse2 -mstackrealign")
  dyn.load(paste0(here::here(), "/", TMB::dynlib(Version)))
  Obj<- MakeADFun(data = TmbData, parameters = Parameters, hessian = FALSE, map = Map, random = Random, inner.method = "newton", DLL = "VAST_v13_0_0")
  Obj$control <- list(parscale = 1, REPORT = 1, reltol = 1e-12, maxit = 100)
  if (FishStatsUtils::convert_version_name(Version) >= FishStatsUtils::convert_version_name("VAST_v4_1_0")) {
    if (Options["normalize_GMRF_in_CPP"] == FALSE) {
      message("Normalizing GMRF in R using `TMB::normalize` feature")
      Obj = TMB::normalize(Obj, flag = "include_data", value = FALSE)
    }
  }
  if (!is.null(DiagnosticDir)) {
    Obj$gr_orig = Obj$gr
    Obj$fn_orig = Obj$fn
    Obj$fn = function(vec) {
      utils::capture.output(matrix(vec, ncol = 1, dimnames = list(names(Obj$par), 
                                                                  NULL)), file = paste0(DiagnosticDir, "fn.txt"))
      utils::write.table(matrix(vec, nrow = 1), row.names = FALSE, 
                         sep = ",", col.names = FALSE, append = TRUE, 
                         file = paste0(DiagnosticDir, "trace.csv"))
      return(Obj$fn_orig(vec))
    }
    Obj$gr = function(vec) {
      utils::capture.output(matrix(vec, ncol = 1, dimnames = list(names(Obj$par), 
                                                                  NULL)), file = paste0(DiagnosticDir, "gr.txt"))
      return(Obj$gr_orig(vec))
    }
    utils::write.table(matrix(Obj$par, nrow = 1), row.names = FALSE, 
                       sep = ",", col.names = FALSE, file = paste0(DiagnosticDir, 
                                                                   "trace.csv"))
  }
  boundsifpresent_fn = function(par, map, name, lower, upper, 
                                bounds) {
    if (name %in% names(par)) {
      bounds[grep(name, names(par)), c("Lower", "Upper")] = rep(1, 
                                                                length(grep(name, names(par)))) %o% c(lower, 
                                                                                                      upper)
    }
    return(bounds)
  }
  Bounds = matrix(NA, ncol = 2, nrow = length(Obj$par), dimnames = list(names(Obj$par), 
                                                                        c("Lower", "Upper")))
  Bounds[, "Lower"] = rep(-Inf, length(Obj$par))
  Bounds[, "Upper"] = rep(Inf, length(Obj$par))
  Bounds[grep("SigmaM", names(Obj$par)), "Upper"] = 10
  if (any(TmbData$ObsModel_ez[1, ] == 8)) 
    Bounds[grep("SigmaM", names(Obj$par)), "Upper"] = 3
  if (!is.null(loc_x) && !is.na(Options_vec["Method"]) && 
      Options_vec["Method"] == 0 && Method != "Spherical_mesh") {
    Dist = stats::dist(loc_x)
    Bounds[grep("logkappa", names(Obj$par)), "Lower"] = log(sqrt(8)/max(Dist))
    Bounds[grep("logkappa", names(Obj$par)), "Upper"] = log(sqrt(8)/min(Dist))
  }
  if (!is.na(Options_vec["Method"]) && Options_vec["Method"] == 
      1 && Method != "Spherical_mesh") {
    Bounds[grep("logkappa", names(Obj$par)), "Upper"] = log(0.9999)
  }
  Bounds = boundsifpresent_fn(par = Obj$par, name = "ln_H_input", 
                              lower = -5, upper = 5, bounds = Bounds)
  Bounds = boundsifpresent_fn(par = Obj$par, name = "gamma1", 
                              lower = -20, upper = 20, bounds = Bounds)
  Bounds = boundsifpresent_fn(par = Obj$par, name = "gamma2", 
                              lower = -20, upper = 20, bounds = Bounds)
  Bounds = boundsifpresent_fn(par = Obj$par, name = "lambda1", 
                              lower = -20, upper = 20, bounds = Bounds)
  Bounds = boundsifpresent_fn(par = Obj$par, name = "lambda2", 
                              lower = -20, upper = 20, bounds = Bounds)
  Bounds = boundsifpresent_fn(par = Obj$par, name = "Beta_rho1", 
                              lower = -0.99, upper = 0.99, bounds = Bounds)
  Bounds = boundsifpresent_fn(par = Obj$par, name = "Beta_rho2", 
                              lower = -0.99, upper = 0.99, bounds = Bounds)
  Bounds = boundsifpresent_fn(par = Obj$par, name = "Beta_rho1_f", 
                              lower = -0.99, upper = 0.99, bounds = Bounds)
  Bounds = boundsifpresent_fn(par = Obj$par, name = "Beta_rho2_f", 
                              lower = -0.99, upper = 0.99, bounds = Bounds)
  Bounds = boundsifpresent_fn(par = Obj$par, name = "Epsilon_rho1", 
                              lower = -0.99, upper = 0.99, bounds = Bounds)
  Bounds = boundsifpresent_fn(par = Obj$par, name = "Epsilon_rho2", 
                              lower = -0.99, upper = 0.99, bounds = Bounds)
  Bounds = boundsifpresent_fn(par = Obj$par, name = "Epsilon_rho1_f", 
                              lower = -0.99, upper = 0.99, bounds = Bounds)
  Bounds = boundsifpresent_fn(par = Obj$par, name = "Epsilon_rho2_f", 
                              lower = -0.99, upper = 0.99, bounds = Bounds)
  Bounds = boundsifpresent_fn(par = Obj$par, name = "rho_c1", 
                              lower = -0.99, upper = 0.99, bounds = Bounds)
  Bounds = boundsifpresent_fn(par = Obj$par, name = "rho_c2", 
                              lower = -0.99, upper = 0.99, bounds = Bounds)
  if (("n_f_input" %in% names(TmbData)) && TmbData[["n_f_input"]] == 
      0) {
    Bounds = boundsifpresent_fn(par = Obj$par, name = "L1_z", 
                                lower = c(-Inf, -0.99), upper = c(Inf, 0.99), bounds = Bounds)
    Bounds = boundsifpresent_fn(par = Obj$par, name = "L2_z", 
                                lower = c(-Inf, -0.99), upper = c(Inf, 0.99), bounds = Bounds)
  }
  if (("OverdispersionConfig" %in% names(TmbData))) {
    if (TmbData[["OverdispersionConfig"]][1] == 0) 
      Bounds = boundsifpresent_fn(par = Obj$par, name = "L1_z", 
                                  lower = c(-Inf, -0.99), upper = c(Inf, 0.99), 
                                  bounds = Bounds)
    if (TmbData[["OverdispersionConfig"]][2] == 0) 
      Bounds = boundsifpresent_fn(par = Obj$par, name = "L2_z", 
                                  lower = c(-Inf, -0.99), upper = c(Inf, 0.99), 
                                  bounds = Bounds)
  }
  Obj$env$inner.control$step.tol <- c(1e-08, 1e-12, 1e-15)[ConvergeTol]
  Obj$env$inner.control$tol10 <- c(1e-06, 1e-08, 1e-12)[ConvergeTol]
  Obj$env$inner.control$grad.tol <- c(1e-08, 1e-12, 1e-15)[ConvergeTol]
  ThorsonUtilities::list_parameters(Obj)
  Return = list(Obj = Obj, Upper = Bounds[, "Upper"], Lower = Bounds[, 
                                                                     "Lower"], Parameters = Parameters, Map = Map, Random = Random)
  class(Return) = "make_model"
  return(Return)
}
