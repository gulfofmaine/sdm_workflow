#####
## Error debugging
#####
library(tidyverse)
library(sf)
library(targets)
#####
## Singularity issue (07/14/2021)
#####

# Summary: After hacking away to find a solution where we can get biomass index abundances for different regions depending on the shapefile, I ran into a memory issue when I tried to do all the regions of interest (strata_use<- data.frame("STRATA" = c("NMFS_and_DFO", "DFO", "NMFS", "Scotian_Shelf", "Georges_Bank", "Gulf_of_Maine", "Southern_New_England", "Mid_Atlantic_Bight"))) so, I ended up reducing the number of regions (strata_use<- data.frame("STRATA" = c("NMFS_and_DFO", "DFO", "NMFS"))). While that solved the memory issue, the final fit returned an error complaining about singularity. I have absolutely no clue what could be causing that error or why it would arise simply from changing the strata. After checking the code quickly, looks like there was one issue with the index_polygon shapefile. Editing that to see how or if things change...No such luck. Trying to remove the joint strata and focus on DFO vs. NMFS. This worked in terms of fitting, but the results plot was really, really bizarre (see Yellowtail flounder2_LogDensity).
setwd(here::here("", "/scratch/aja/TargetsSDM"))
temp_res<- read.csv(here::here("", "/scratch/aja/TargetsSDM/results/tables/Biomass_Index_raw_Yellowtail flounder2.csv"))

plot(temp_res$Time, temp_res$Index_Estimate)

# Just ONE ridiculous outlier...
outlier<- which.max(temp_res$Index_Estimate)
plot(temp_res$Time[-outlier], temp_res$Index_Estimate[-outlier])

temp_res[outlier,]

# 31 time step -- what year was this? Was there something weird in the raw data?
unique(temp_res$Time)

tar_load(tidy_mod_data)

plot_catch<- tidy_mod_data %>%
  filter(., NMFS_SVSPP == 101) %>%
  group_by(., EST_YEAR, SURVEY) %>%
  summarize(., "Total_Bio" = sum(BIOMASS),
            "Total_Abund" = sum(ABUNDANCE))

ggplot() +
  geom_point(data = plot_catch, aes(x = EST_YEAR, y = Total_Bio, color = SURVEY)) +
  theme_bw()

ggplot() +
  geom_point(data = plot_catch, aes(x = EST_YEAR, y = Total_Abund, color = SURVEY)) +
  theme_bw()

# Nothing major there. The time steps, though, seem really bizarre...not sure why they aren't just printed regularly. Oh well..
View(temp_res)

# In the first entry for DFO, we've got time steps every 
temp_res2<- temp_res %>%
  arrange(., Time)

# What is going on with the predictions?
tar_load(vast_predictions)
summary(vast_predictions)

temp<- vast_predictions %>%
  arrange(-Pred)

# 2009-spring seems to be the culprit...
temp2<- tidy_mod_data %>%
  filter(., EST_YEAR == 2009 & SEASON == "SPRING" & NMFS_SVSPP == 105)
summary(temp2)
temp2[which.max(temp2$BIOMASS), ]


# VAST fit
tar_load(vast_fit)
TMBhelper::check_estimability(vast_fit$tmb_list$Obj)
summary(vast_fit$data_list$a_i)


diagnose_hessian <- function(fit,h=NULL, eval.eps=1e-5,evec.eps=1e-2) {
  if(FALSE){
    fit = vast_fit$tmb_list
    h = NULL
    eval.eps=1e-5
    evec.eps=1e-2
  }
  ## pull out the TMB object from the fit
  obj <- fit$Obj
  ee <- environment(obj$fn)
  ## extract parameters
  pp <- ee$last.par[-ee$random]
  ## easiest way to get names corresponding to all of the parameters
  nn <- tryCatch(colnames(vcov(fit, full=TRUE)),
                 ## fall-back position
                 error = function(e) make.unique(names(pp)))
  ## fit$sdr$pdHess
  if ("sdr" %in% names(fit)) {
    cat("bad params according to sdreport:",
        paste(nn[!is.finite(suppressWarnings(sqrt(diag(fit$sdr$cov.fixed))))],
              collapse=", "),"\n")
  }
  ## two ways to compute the Hessian
  ## (1) directly from the objective function, via finite difference+Richardson extrapolation
  ## h1 <- hessian(obj$fn, pp)
  ## (2) use the gradient and compute its Jacobian (faster and probably more stable)
  if (is.null(h)) {
    if (!require(numDeriv)) stop("need numDeriv package installed") 
    h <- jacobian(obj$gr, pp)
  }
  ## double-check we get the same answer (approximately)
  ## all.equal(h1,h,tolerance=1e-5)
  ## now investigate the Hessian
  eigs <- eigen(h)
  ## non-positive definite means some of the eigenvectors are <= 0
  bad <- which(eigs$values/max(eigs$values)<=eval.eps)
  if (length(bad)==0) {
    cat("Hessian seems OK\n")
    return(invisible(h))
  }
  cat(sprintf("max eigenvalue = %1.3g", eigs$values[1]), "\n")
  for (b in bad) {  ## there could be more than one 'bad' direction/eigenvector ..
    cat(sprintf("Hessian eigenvalue %d = %1.3g (relative val = %1.3g)",
                b,eigs$values[b], eigs$values[b]/eigs$values[1]), "\n")
    bad_vec <- eigs$vectors[,b]
    bad_elements <- which(abs(bad_vec) > evec.eps)
    cat("   bad elements:", nn[bad_elements],"\n")
  }
  cat("SDs computed from sqrt(diag(solve(H))):",
      paste(suppressWarnings(sqrt(diag(solve(h)))), collapse=", "),"\n")
  return(invisible(h))
}
diagnose_hessian(m0)
