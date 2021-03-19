##########
##### Executing _targets.R
##########
library(targets)

# First, need to be in the right working directory
getwd()
setwd("~/GitHub/sdm_workflow/scratch/aja/targets_flow")

# Checking calls
tar_manifest(fields = "command")

# Graph
tar_glimpse()
tar_visnetwork(label = "time", reporter = "forecast", targets_only = )

# Run it
tar_make()

