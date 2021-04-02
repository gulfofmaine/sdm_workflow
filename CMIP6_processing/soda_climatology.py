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
start_year = 1985
end_year   = 2014


# Path to SODA data sources on BOX
soda_path = {"surf_sal"  : f"{box_root}RES_Data/SODA/SODA_Salt_Red.nc",
             "bot_sal"   : f"{box_root}RES_Data/SODA/SODA_Salt_Red_bottomLayer.nc",
             "surf_temp" : f"{box_root}RES_Data/SODA/SODA_Temp_Red.nc",
             "bot_temp"  : f"{box_root}RES_Data/SODA/SODA_Temp_Red_bottomLayer.nc"}

# key to match variable name to SODA source
var_key = {"bot_sal"    : "salt",
           "surf_sal"   : "salt",
           "bot_temp"   : "temp",
           "surf_temp"  : "temp"}
           


####  Loading Data  ####

# Instead of loop, just put them all together in one...
surftemp = xr.open_dataset(soda_path["surf_temp"])
surftemp = surftemp.rename_vars({"temp" : "surf_temp"}).drop("st_ocean")

bottemp = xr.open_dataset(soda_path["bot_temp"])
bottemp = bottemp.rename_vars({"temp" : "bot_temp"}).drop("st_ocean")

surfsal = xr.open_dataset(soda_path["surf_sal"])
surfsal = surfsal.rename_vars({"salt" : "surf_sal"}).drop("st_ocean")

botsal = xr.open_dataset(soda_path["bot_sal"])
botsal = botsal.rename_vars({"salt" : "bot_sal"}).drop("st_ocean")





####  Build xr.Dataset from pieces of each  ####


# Drop the depth dimension where it exists, then rebuild
def remove_st_ocean(xr_ds, var):
  """Pull out data as an array, drop st_ocean dimension, rebuild xr.array. 
  Need to pull surface measurement from surface data arrays so depth coordinate
  becomes unnecessary.
  
  Args:
    xr_ds      : xr.ArrayDataset
    var (str)  : String indicating variable to pull and process
  
  """
  
  # Take the data values out as an array
  data = xr_ds[var].values 
  
  # Take all values from each dimension EXCEPT depth
  if var in ["surf_temp", "surf_sal"]:
    data = data[:, 0, :, :] 
  
  # Pull the coordinates
  time     = xr_ds.coords['time']
  xt_ocean = xr_ds.coords['xt_ocean']
  yt_ocean = xr_ds.coords['yt_ocean']

  # Building an xr.array
  xr_manual = xr.DataArray(data, coords = [time, yt_ocean, xt_ocean])
  return xr_manual



# Remove st_ocean from each of them
surftemp_manual = remove_st_ocean(surftemp, "surf_temp")
bottemp_manual  = remove_st_ocean(bottemp, "bot_temp")
surfsal_manual  = remove_st_ocean(surfsal, "surf_sal")
botsal_manual   = remove_st_ocean(botsal, "bot_sal")


# Plot check
botsal_manual.isel(time = 0).plot()
plt.show()
plt.clf()




####  Build a dataset from all four  ####
soda_manual = xr.Dataset({"surf_temp" : surftemp_manual,
                          "surf_sal"  : surfsal_manual,
                          "bot_temp"  : bottemp_manual,
                          "bot_sal"   : botsal_manual})



# Inspect it
soda_full = soda_manual
type(soda_full)

# plot
soda_full.bot_sal.isel(time = 0).plot()
plt.show()
plt.clf()



####  Screening Odd Values  ####

# This step is to catch instances where unreasonably
# high or low values are used to flag NA's

def drop_extremes(xr_ds, soda_var, clamp_low, clamp_high):
  """
  Drop values outside of desired range from xarray dataset.
  
  Args:
    xr_ds              : xr.ArrayDataset to apply data screen to
    soda_var (str)     : variable to screen 
    clamp_low (float)  : lower limit to use for data screen
    clamp_high (float) : upper limit to use for data screen
  
  """
  # Use soda var to pull attribute data
  soda_da = xr_ds[soda_var]
  
  # Set condition flags
  lower_cond = soda_da > clamp_low
  upper_cond = soda_da <  clamp_high
  
  # Reassign to nan where those flags are not true
  soda_da = soda_da.where(lower_cond) 
  soda_da = soda_da.where(upper_cond) 
  
  xr_ds[soda_var] = soda_da
  return xr_ds



# test it
soda_filtered = drop_extremes(xr_ds = soda_full, soda_var = "bot_temp", clamp_low = -50, clamp_high = 100)
soda_filtered = drop_extremes(xr_ds = soda_filtered, soda_var = "bot_sal", clamp_low = -50, clamp_high = 100)
soda_filtered = drop_extremes(xr_ds = soda_filtered, soda_var = "surf_temp", clamp_low = -50, clamp_high = 100)
soda_filtered = drop_extremes(xr_ds = soda_filtered, soda_var = "surf_sal", clamp_low = -50, clamp_high = 100)
soda_filtered







####  Subset Years  ####
soda_slice = soda_filtered.sel(time = slice(f"{start_year}-01-01", f"{end_year}-12-31"))


# Use groupby to get average for each month in reference period
months = soda_slice.time.dt.month # grouping indices
monthly_clim = soda_slice.groupby(months).mean('time', keep_attrs = True)

# does it work without specifying a variable?
monthly_clim.sel(month = 2)["bot_sal"].plot()
plt.title("February Avg. Bottom Salinity")
plt.show()
plt.clf()



# What does the timeline for the whole area look like?:
area_ts = getattr(monthly_clim, "surf_sal").mean(dim = ("yt_ocean", "xt_ocean"))
area_ts.plot()
plt.show()
plt.clf()



####  Save Climatology

# Build out attributes
monthly_clim.attrs = {
  'title'                   : f'MOM5_SODA_3.4.2 - Monthly Climatology {start_year}-{end_year}', 
  'Date Created'            : "2/18/2021",
  'Institution'             : "Gulf  of Maine Research Institute",
  'grid_type'               : 'mosaic', 
  'CDI'                     : 'Climate Data Interface version 1.9.9 (https://mpimet.mpg.de/cdi)', 
  'Conventions'             : 'CF-1.6', 
  'history'                 : 'Fri Feb 12 08:05:52 2021: cdo sellonlatbox,-120,30,0,80 SODA_Temp.nc SODA_Temp_Red.nc\nFri Feb 12 07:14:41 2021: cdo cat sst_soda3.4.2_mn_ocean_reg_1980.nc sst_soda3.4.2_mn_ocean_reg_1981.nc sst_soda3.4.2_mn_ocean_reg_1982.nc sst_soda3.4.2_mn_ocean_reg_1983.nc sst_soda3.4.2_mn_ocean_reg_1984.nc sst_soda3.4.2_mn_ocean_reg_1985.nc sst_soda3.4.2_mn_ocean_reg_1986.nc sst_soda3.4.2_mn_ocean_reg_1987.nc sst_soda3.4.2_mn_ocean_reg_1988.nc sst_soda3.4.2_mn_ocean_reg_1989.nc sst_soda3.4.2_mn_ocean_reg_1990.nc sst_soda3.4.2_mn_ocean_reg_1991.nc sst_soda3.4.2_mn_ocean_reg_1992.nc sst_soda3.4.2_mn_ocean_reg_1993.nc sst_soda3.4.2_mn_ocean_reg_1994.nc sst_soda3.4.2_mn_ocean_reg_1995.nc sst_soda3.4.2_mn_ocean_reg_1996.nc sst_soda3.4.2_mn_ocean_reg_1997.nc sst_soda3.4.2_mn_ocean_reg_1998.nc sst_soda3.4.2_mn_ocean_reg_1999.nc sst_soda3.4.2_mn_ocean_reg_2000.nc sst_soda3.4.2_mn_ocean_reg_2001.nc sst_soda3.4.2_mn_ocean_reg_2002.nc sst_soda3.4.2_mn_ocean_reg_2003.nc sst_soda3.4.2_mn_ocean_reg_2004.nc sst_soda3.4.2_mn_ocean_reg_2005.nc sst_soda3.4.2_mn_ocean_reg_2006.nc sst_soda3.4.2_mn_ocean_reg_2007.nc sst_soda3.4.2_mn_ocean_reg_2008.nc sst_soda3.4.2_mn_ocean_reg_2009.nc sst_soda3.4.2_mn_ocean_reg_2010.nc sst_soda3.4.2_mn_ocean_reg_2011.nc sst_soda3.4.2_mn_ocean_reg_2012.nc sst_soda3.4.2_mn_ocean_reg_2013.nc sst_soda3.4.2_mn_ocean_reg_2014.nc sst_soda3.4.2_mn_ocean_reg_2015.nc sst_soda3.4.2_mn_ocean_reg_2016.nc sst_soda3.4.2_mn_ocean_reg_2017.nc sst_soda3.4.2_mn_ocean_reg_2018.nc sst_soda3.4.2_mn_ocean_reg_2019.nc SODA_Temp.nc\nThu Feb 11 18:55:43 2021: cdo selvar,temp soda3.4.2_mn_ocean_reg_1980.nc sst_soda3.4.2_mn_ocean_reg_1980.nc\nTue Nov  6 14:21:13 2018: ncrcat /aosc/greenland/soda3.4.2/REGRIDED/ocean/soda3.4.2_mn_ocean_reg_1980_01.nc /aosc/greenland/soda3.4.2/REGRIDED/ocean/soda3.4.2_mn_ocean_reg_1980_02.nc /aosc/greenland/soda3.4.2/REGRIDED/ocean/soda3.4.2_mn_ocean_reg_1980_03.nc /aosc/greenland/soda3.4.2/REGRIDED/ocean/soda3.4.2_mn_ocean_reg_1980_04.nc /aosc/greenland/soda3.4.2/REGRIDED/ocean/soda3.4.2_mn_ocean_reg_1980_05.nc /aosc/greenland/soda3.4.2/REGRIDED/ocean/soda3.4.2_mn_ocean_reg_1980_06.nc /aosc/greenland/soda3.4.2/REGRIDED/ocean/soda3.4.2_mn_ocean_reg_1980_07.nc /aosc/greenland/soda3.4.2/REGRIDED/ocean/soda3.4.2_mn_ocean_reg_1980_08.nc /aosc/greenland/soda3.4.2/REGRIDED/ocean/soda3.4.2_mn_ocean_reg_1980_09.nc /aosc/greenland/soda3.4.2/REGRIDED/ocean/soda3.4.2_mn_ocean_reg_1980_10.nc /aosc/greenland/soda3.4.2/REGRIDED/ocean/soda3.4.2_mn_ocean_reg_1980_11.nc /aosc/greenland/soda3.4.2/REGRIDED/ocean/soda3.4.2_mn_ocean_reg_1980_12.nc -o /aosc/greenland/soda3.4.2/REGRIDED/ocean/soda3.4.2_mn_ocean_reg_1980.nc\nTue Nov  6 14:16:16 2018: ncra --rth_flt -v temp,salt,u,v,wt,prho,ssh,mlt,mlp,mls,net_heating,taux,tauy /aosc/greenland/soda3.4.2/REGRIDED/ocean/soda3.4.2_5dy_ocean_reg_1980_01_03.nc /aosc/greenland/soda3.4.2/REGRIDED/ocean/soda3.4.2_5dy_ocean_reg_1980_01_08.nc /aosc/greenland/soda3.4.2/REGRIDED/ocean/soda3.4.2_5dy_ocean_reg_1980_01_13.nc /aosc/greenland/soda3.4.2/REGRIDED/ocean/soda3.4.2_5dy_ocean_reg_1980_01_18.nc /aosc/greenland/soda3.4.2/REGRIDED/ocean/soda3.4.2_5dy_ocean_reg_1980_01_23.nc /aosc/greenland/soda3.4.2/REGRIDED/ocean/soda3.4.2_5dy_ocean_reg_1980_01_28.nc -o /aosc/greenland/soda3.4.2/REGRIDED/ocean/soda3.4.2_mn_ocean_reg_1980_01.nc\nTue Aug 14 14:15:32 2018: cdo -f nc -b F32 copy /aosc/greenland/soda3.4.2/REGRIDED/ocean/soda3.4.2_5dy_ocean_reg_1980_01_03.nc /aosc/greenland/soda3.4.2/REGRIDED/ocean/f32/soda3.4.2_5dy_ocean_reg_1980_01_03.nc', 'filename': './regrid_MOM2SODA.nc                ', 
  'grid_tile'               : '1', 
  'nco_openmp_thread_number': 1, 
  'CDO'                     : 'Climate Data Operators version 1.9.9 (https://mpimet.mpg.de/cdo)'
}

# Review Climatology
monthly_clim

# save to box
monthly_clim.to_netcdf(f"{box_root}RES_Data/SODA/SODA_monthly_climatology{start_year}to{end_year}.nc")
