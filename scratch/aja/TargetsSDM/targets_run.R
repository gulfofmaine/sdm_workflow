##########
##### Executing _targets.R
##########
library(targets)
library(parallel)
library(doFuture)
library(tictoc)
# 
cores_avail<- detectCores()
registerDoFuture()
plan(multisession, workers = cores_avail-2)

# Clean everything?
clean_start<- FALSE
if(clean_start){
  tar_destroy()
}

# First, need to be in the right working directory
getwd()
setwd("~/GitHub/sdm_workflow/scratch/aja/TargetsSDM")

# Checking calls
tar_manifest(fields = "command")

# Graph
# tar_glimpse()
# tar_visnetwork(label = "time", reporter = "forecast", targets_only = TRUE)

# Run it
tic()
tar_make()
toc()

# Check on warnings
warning_ind<- which(!is.na(tar_meta(fields = warnings)$warnings))
tar_meta(fields = warnings)[warning_ind, ]

tar_load(vast_fit)
TMBhelper::check_estimability(vast_fit$tmb_list$Obj)


