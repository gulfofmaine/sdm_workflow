# Libraries
import os
import xarray as xr
import numpy as np
import pandas as pd


####  CMIP Model Date Key construction  ####


####  Purpose:
# Due to the unusual calendar type used in these models R hass not been reading in dates 
# that are consistent with their values on xarray
# this script builds lookup tables to use for cross-referencing the dates back into R


####  Goal:
# Load all the CMIP6 runs, put their time dimensions in a table with the filenames
# Use this to ensure R is handling them correctly


####  Set up workspace and paths  ####

# Set Workspace:
workspace = "local"

# Set username
UsrName='mdzaugis'

# Set spp experiment
experiment='ssp1_26'

# Root paths for sdm_workflows project - local/docker
root_locations = {
  "local": f"/Users/{UsrName}/Box/",
  "docker": "/home/jovyan/"}

# Set root based on workspace
box_root = root_locations[workspace]
print(f"Working via {workspace} directory at: {box_root}")


# Path to cmip data sources on BOX
cmip_path = {"surf_sal"  : f"{box_root}RES_Data/CMIP6/{experiment}/SurSalinity/StGrid/",
             "bot_sal"   : f"{box_root}RES_Data/CMIP6/{experiment}/BottomSal/StGrid/",
             "surf_temp" : f"{box_root}RES_Data/CMIP6/{experiment}/SST/StGrid/",
             "bot_temp"  : f"{box_root}RES_Data/CMIP6/{experiment}/BottomT/StGrid/"}


# Pick a variable
cmip_var = "surf_temp"

           
           
           
####  Road Map  ####

####  check all the files in the folders to get the time index for each


# 1. Put full file paths into a list
d = cmip_path[cmip_var]
ssal_files = []
full_paths = []
for file_name in os.listdir(d):
    ssal_files.append(file_name)
    full_path = os.path.join(d, file_name)
    full_paths.append(full_path)

    
# 2. Use zip to loop through the paths and the files
time_indices = []
for full_file, short_name in zip(full_paths, ssal_files):
  
  # load xarray dataset
  try:
    cmip_xr = xr.open_dataset(full_file)
    time_index = cmip_xr.get_index("time").to_numpy()
    time_indices.append(time_index)
    cmip_xr.close()
  except:
    time_indices.append("This file is jank")
    
  # Print the window ranges
  start_date = time_index[0]
  end_date   = time_index[-1]
  print(f"{short_name} time window:")
  print(f"{start_date} - {end_date}")



####  Get the different lengths to determine if historic or projection run  ####


# 1. Get the different lengths
step_lengths = []
for steps in time_indices:
  step_len = len(steps)
  step_lengths.append(step_len)
  
# unique lengths
set(step_lengths)



# 2. Print which ones are strange
for step_len, short_name in zip(step_lengths, ssal_files):
  if step_len not in [780, 1032]:
    print(f"Strange time dimension for {short_name}")
    print(f"For Variable: {cmip_var}")
    print(f"Time dimension length: {step_len}")
  


####  Build Tables for Date Keys  ####


# 1. Put common lengths into tables to use as keys
historic_runs = {}  # should be 780
projection_runs = {}# should be 1032
under_run = {}      # These are ones that don't load at all
over_run = {}       # these are ones that extend beyond what they should
for time_index, short_name in zip(time_indices, ssal_files):
  short_name = short_name.replace(".nc", "")
  if len(time_index) == 780:
    historic_runs[f"{short_name}"] = (time_index)
  if len(time_index) == 1032:
    projection_runs[f"{short_name}"]  = (time_index)
  if len(time_index) > 1032:
    over_run[f"{short_name}"]  = (time_index)
    print(f"{short_name} is jank and has time length of {len(time_index)}")
  if len(time_index) < 780:
    under_run[f"{short_name}"]  = (time_index)
    print(f"{short_name} is jank and has time length of {len(time_index)}")
    


# 2. Convert Dictionaries to pd.dataframe
historic_df = pd.DataFrame.from_dict(historic_runs)  # Historical runs only need to be cataloged once
project_df  = pd.DataFrame.from_dict(projection_runs)
over_df  = pd.DataFrame.from_dict(over_run)

# check if columns are the same across columns, not quite
historic_df.iloc[:, [1]].equals(historic_df.iloc[:, [2]])
historic_df.iloc[:, [2]]


####  Save Keys to Box  ####

# save them to box
folder_location = f"{box_root}RES_Data/CMIP6/{experiment}/DateKeys/"
historic_df.to_csv(f"{folder_location}{cmip_var}/{cmip_var}_historic_runs.csv", index = False)
project_df.to_csv(f"{folder_location}{cmip_var}/{cmip_var}_future_projections.csv", index = False)
over_df.to_csv(f"{folder_location}{cmip_var}/{cmip_var}_over_run.csv", index = False)

