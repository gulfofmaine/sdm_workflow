#!/usr/bin/env python
# “Regrid_script_execute.py”
import cdo                       # for the cdo tools
import glob                      # for file listing
import xarray as xr              # for netcdfs
import numpy as np               # array tools
import fcts
import os
import re
import shutil
import xesmf as xe

cdo = cdo.Cdo()                      # need to initialize a cdo object

print('Enter User Name')
UsrName = input()
print('Enter Dest Group (i.e. RES_Data)')
Group = input()
print('Enter Folder with a / (i.e. CMIP6/)')
Folder = input()

UsrName='mdzaugis'
Group='RES_Data'
Folder='CMIP6/'
ExperimentFolder='CMIP6/SSP1_26/'
Experiment='ssp126'

CMIPpath = fcts.shared_path(user_name=UsrName, group=Group, folder=Folder)
ExperimentPath = fcts.shared_path(user_name=UsrName, group=Group, folder=ExperimentFolder)

dsTOS = glob.glob(f'{ExperimentPath}RawTmpFiles/tos*_{Experiment}*')
dsTOS = glob.glob(f'{ExperimentPath}RawTmpFiles/tos*_historical*')
for file in dsTOS:
    gridfi = f'{CMIPpath}GridFiles/OISST_grid.nc'
    base_filename = os.path.basename(file)
    fileout = f'{ExperimentPath}SST/tmpfiles/stGrid_{base_filename}'
    cdo.remapdis(gridfi,  input=file, output=fileout, options='-f nc')
    ds = xr.open_dataset(fileout)
    savepath = f'{ExperimentPath}SST/StGrid/stGrid_{base_filename}'
    ds_cropped = ds.sel(lon=slice(260, 320), lat=slice(20, 70))
    ds_cropped.to_netcdf(savepath)
    os.remove(fileout)
    print(f'Completed regridding {base_filename} out of {str(len(dsTOS))}')


dsThetao = glob.glob(f'{ExperimentPath}RawTmpFiles/thetao*_{Experiment}*')
dsThetao = glob.glob(f'{ExperimentPath}RawTmpFiles/thetao*_historical*')
file = "/Users/mdzaugis/Box/RES_Data/CMIP6/RawTmpFiles/thetao_CESM2_r4i1p1f1_historical.nc"

## Completed regridding thetao_MRI-ESM2-0_r1i1p1f1_ssp126.nc out of 26 need to find more space to run on box
for file in dsThetao:
    gridfi = f'{CMIPpath}GridFiles/SODA_grid.nc'
    base_filename = os.path.basename(file)
    fileout = f'{ExperimentPath}BottomT/tmpfiles/stGrid_{base_filename}'
    cdo.remapdis(gridfi,  input=file, output=fileout, options='-f nc')
    ds = xr.open_dataset(fileout)
    savepath = f'{ExperimentPath}BottomT/StGrid/stGrid_{base_filename}'
    ds_cropped = ds.sel(longitude=slice(260, 320), latitude=slice(20, 70))
    #ds_cropped = ds_cropped.isel(time=slice(None, 1032))
    ds_cropped.to_netcdf(savepath)
    os.remove(fileout)
    print(f'Completed regridding {base_filename} out of {str(len(dsThetao))}')

dsSO = glob.glob(f'{ExperimentPath}RawTmpFiles/so*_{Experiment}*')
dsSO = glob.glob(f'{ExperimentPath}RawTmpFiles/so*_historical*')
file = "/Users/mdzaugis/Box/RES_Data/CMIP6/RawTmpFiles/so_MRI-ESM2-0_r1i1p1f1_ssp585.nc"
for file in dsSO:
    gridfi = f'{CMIPpath}GridFiles/SODA_grid.nc'
    base_filename = os.path.basename(file)
    fileout = f'{ExperimentPath}BottomSal/tmpfiles/stGrid_{base_filename}'
    cdo.remapdis(gridfi,  input=file, output=fileout, options='-f nc')
    ds = xr.open_dataset(fileout)
    savepath = f'{ExperimentPath}BottomSal/StGrid/stGrid_{base_filename}'
    ds_cropped = ds.sel(longitude=slice(260, 320), latitude=slice(20, 70))
    #ds_cropped = ds_cropped.isel(time=slice(None, 1032))
    ds_cropped.to_netcdf(savepath)
    os.remove(fileout)
    print(f'Completed regridding {base_filename} out of {str(len(dsSO))}')

dsSurSO = glob.glob(f'{ExperimentPath}RawTmpFiles/Surface*_{Experiment}*')
dsSurSO = glob.glob(f'{ExperimentPath}RawTmpFiles/Surface*_historical*')
file = "/Users/mdzaugis/Box/RES_Data/CMIP6/RawTmpFiles/Surface_so_MRI-ESM2-0_r1i1p1f1_ssp585.nc"
for file in dsSurSO:
    gridfi = f'{CMIPpath}GridFiles/SODA_grid.nc'
    base_filename = os.path.basename(file)
    fileout = f'{ExperimentPath}SurSalinity/tmpfiles/stGrid_{base_filename}'
    cdo.remapdis(gridfi,  input=file, output=fileout, options='-f nc')
    ds = xr.open_dataset(fileout, decode_times=False)
    savepath = f'{ExperimentPath}SurSalinity/StGrid/stGrid_{base_filename}'
    ds_cropped = ds.sel(longitude=slice(260, 320), latitude=slice(20, 70))
    ds_cropped.to_netcdf(savepath)
    os.remove(fileout)
    print(f'Completed regridding {base_filename} out of {str(len(dsSurSO))}')
