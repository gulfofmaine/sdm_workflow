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



#### What does it look like in Rstudio  ####
test_day = soda.isel(time = 0) #first time step
test_day.temp.plot()
plt.show()

#close plot
plt.clf()



#### Subset Years  ####
soda = soda.sel(time = slice(f"{start_year}-01-01", f"{end_year}-12-31"))


####  Group on Month to get climatology
months = soda.time.dt.month
monthly_clim_temp = soda['temp'].groupby(months).mean('time', keep_attrs = True)

