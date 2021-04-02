#### Comparing Bias Corrections to Observations
# Goal: Side-by-side display of CMIP6 Bias Corrected Data to
# observational data resources for a region/area of interest

####  Load Libraries  
import xarray as xr
import cartopy.crs as crs
import matplotlib.pyplot as plt
import numpy as np
import netCDF4
import regionmask
import geopandas as gpd
import os


####  Load Data
box_root = "/Users/akemberling/Box/RES_Data/"


####  Load Bias Corrected CMIP6 Variables  
cmip_root = f"{box_root}CMIP6/BiasCorrected/"
cmip_surf_temp = xr.open_dataset(f"{cmip_root}surf_temp_OISST_bias_corrected_mean.grd")
cmip_surf_sal  = xr.open_dataset(f"{cmip_root}surf_sal_SODA_bias_corrected_mean.grd")
cmip_bot_temp  = xr.open_dataset(f"{cmip_root}bot_temp_SODA_bias_corrected_mean.grd")
cmip_bot_sal   = xr.open_dataset(f"{cmip_root}bot_sal_SODA_bias_corrected_mean.grd")


####  Load Observation Counterparts 


# 1. SODA
soda_root  = f"{box_root}SODA/"
soda_files = {"bot_sal"  : "SODA_Salt_Red_bottomLayer.nc", 
              "surf_sal" : "SODA_Salt_Red.nc", 
              "bot_temp" : "SODA_Temp_Red_bottomLayer.nc"}
              
bot_sal_data  = xr.open_dataset(f"{soda_root}{soda_files['bot_sal']}")
bot_temp_data = xr.open_dataset(f"{soda_root}{soda_files['bot_temp']}")
surf_sal_data = xr.open_dataset(f"{soda_root}{soda_files['surf_sal']}")

# 2. OISST
oisst_root = f"{box_root}OISST/oisst_mainstays/annual_observations/"

# start and end years for the update
start_yr = 1982
end_yr = 2020

# load the annual files for oisst
fpaths = []
for yr in range(start_yr, end_yr + 1):
    fpaths.append(f'{oisst_location}sst.day.mean.{yr}.v2.nc')
    
# Lazy-load using xr.open_mfdataset
surf_sal_data = xr.open_mfdataset(fpaths, combine = "by_coords", parallel = False)


# 3. Load NE Shelf Shapefile
mask_shape = gpd.read_file(f"{box_root}Shapefiles/NELME_regions/NELME_sf.shp")



observation_vars = [[bot_sal_data], [bot_temp_data], [surf_sal_data], [surf_temp_data]]
