###  Processing SODA Climatology with XARRAY  ####

####  Libraries  ####
import xarray as xr
import cartopy.crs as crs
import matplotlib.pyplot as plt
import numpy as np
import netCDF4
import os



####  Set up workspace and paths  ####

# Set Workspace:
workspace = "local"

# Root paths for sdm_workflows project - local/docker
root_locations = {
  "local" : "/Users/akemberling/Box/",
  "docker": "/home/jovyan/"}

# Set root based on workspace
box_root = root_locations[workspace]
print(f"Working via {workspace} directory at: {box_root}")

# Set start and end year for climatology
start_year = 1990
end_year   = 2019

# Path to SODA on BOX
soda_path = f"{box_root}RES_Data/SODA/soda3.4.2_mn_ocean_reg_1980.nc"



####  Loading Data  ####
soda = xr.open_dataset(soda_path)

# Coordinates
soda.coords

# Variables
soda.variables

# Dimensions
soda.dims
soda.xt_ocean
soda.yt_ocean
soda.st_ocean

# Attributes
soda.attrs

# Coordinate ordering is time, depth, lat, lon



# #### What does it look like in Rstudio  ####

# First time step, surface layer
test_day = soda.isel(time = 0, st_ocean = 0) 
test_day.temp.plot()
plt.show()

#close plot
plt.clf()

# #first time step, surface layer
test_day.salt.plot()
plt.show()

#close plot
plt.clf()


####  Slicing Surface / Bottom  ####

# Idea: since depth doesnt change, take a single time tep, then get the max depth? or will there be 
# empty cells for all depths







#### Subset Years  ####

# The date format for slice is not right, but the approach is good
#soda = soda.sel(time = slice(f"{start_year}-01-01", f"{end_year}-12-31"))


####  Group on Month to get climatology

# Grouping index for months
months = soda.time.dt.month

# Depth layer indexes:
# surf_depths
# bot depths

# as written, ignores depth
monthly_clim_temp = soda['temp'].groupby(months).mean('time', keep_attrs = True)
monthly_clim_temp = soda['temp'].groupby(months).mean('time', keep_attrs = True)


# What are we looking at?
monthly_clim_temp.sel(month = 2)
#plt.show()


# Notes
# need to sort out the surface/bottom slicing
# surface should be easy, bottom will be tricky because bottom slice depends on topography


