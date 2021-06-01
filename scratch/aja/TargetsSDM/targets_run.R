##########
##### Executing _targets.R
##########
library(targets)
library(parallel)
library(doFuture)

cores_avail<- detectCores()
registerDoFuture()
plan(multisession, workers = cores_aval-2)

# Clean everything?
clean_start<- TRUE
if(clean_start){
  tar_destroy()
}

# First, need to be in the right working directory
getwd()
setwd("~/GitHub/sdm_workflow/scratch/aja/TargetsSDM")

# Checking calls
tar_manifest(fields = "command")

# Graph
tar_glimpse()
tar_visnetwork(label = "time", reporter = "forecast", targets_only = TRUE)

# Run it
tar_make()

