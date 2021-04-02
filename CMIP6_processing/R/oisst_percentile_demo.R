
####  OISST ClusterR - Getting Quantiles for every year for OISST  
####  2/9/2020
####  Adam Kemberling
####  Goal: Use oisst stack as a tester for processing 5th, mean, 95th percentile arrays


####  Packages  ####
library(raster)
library(gmRi)
library(here)
library(tidyverse)
library(snow)

#  color palette for quick raster displays
temp_pal <- rev(RColorBrewer::brewer.pal(n = 10, name = "RdBu"))


####  Data  ####

# SST from gmRi::
box_paths <- research_access_paths(os.use = "unix")
oisst_path <- box_paths$oisst_mainstays
data_window <- data.frame(lon = c(-100, -60),
                          lat = c(20, 70),
                          time = as.Date(c("2019-01-01", "2020-12-31")))

oisst <- oisst_window_load(oisst_path = oisst_path, data_window = data_window, anomalies = FALSE)


# Starting Point - 2 years
# Objective: Get mean, 5%, 95% raster layer for each month of data, on cell-by-cell basis




####  Functions  ####
source(here("CMIP6_processing/R/sdm_workflow_funs.R"))





####  Workflow  ####
# Map through years to return monthly stats
monthly_percentiles <- map(oisst, timestep_stats)




# Plot one
plot(monthly_percentiles$`2019`$percentile_05$X01, 
     main = "2019 - January\n5th Percentile Temps\nOISSTv2", col = temp_pal)




year_labels <- c("2019", "2020")
all_years_05   <- map2(monthly_percentiles, year_labels, timestep_to_full, stat_group = "percentile_05") %>% stack()
all_years_mean <- imap(monthly_percentiles, timestep_to_full, stat_group = "mean") %>% stack()
all_years_95   <- map2(monthly_percentiles, year_labels, timestep_to_full, stat_group = "percentile_95") %>% stack()



  
  
####  Using clusterR  ####
# # Use clusterR to speed it up - slower for small stacks
# 
# # Get mean + percentile rasters
# beginCluster(4)
# ras.mean <- clusterR(month_subset_ras, calc, args = list(mean, na.rm = T))
# ras.05   <- clusterR(month_subset_ras, calc, args = list(quantile, probs = c(.05), na.rm = T))
# ras.95   <- clusterR(month_subset_ras, calc, args = list(quantile, probs = c(.95), na.rm = T))
# endCluster()


