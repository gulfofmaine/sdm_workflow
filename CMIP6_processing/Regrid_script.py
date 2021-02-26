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

path = fcts.shared_path(user_name=UsrName, group=Group, folder=Folder)

dsTOS = glob.glob(f'{path}RawTmpFiles/tos*_ssp585*')
dsTOS = glob.glob(f'{path}RawTmpFiles/tos*_historical*')
for file in dsTOS:
    gridfi = f'{path}GridFiles/OISST_grid.nc'
    base_filename = os.path.basename(file)
    fileout = f'{path}SST/tmpfiles/stGrid_{base_filename}'
    cdo.remapdis(gridfi,  input=file, output=fileout, options='-f nc')
    ds = xr.open_dataset(fileout)
    savepath = f'{path}SST/StGrid/stGrid_{base_filename}'
    ds_cropped = ds.sel(lon=slice(260, 320), lat=slice(20, 70))
    ds_cropped.to_netcdf(savepath)
    os.remove(fileout)
    print(f'Completed regridding {base_filename} out of {str(len(dsTOS))}')

dsThetao = glob.glob(f'{path}RawTmpFiles/thetao*_ssp585*')
dsThetao = glob.glob(f'{path}RawTmpFiles/thetao*_historical*')
for file in dsThetao:
    gridfi = f'{path}GridFiles/SODA_grid.nc'
    base_filename = os.path.basename(file)
    fileout = f'{path}BottomT/tmpfiles/stGrid_{base_filename}'
    cdo.remapdis(gridfi,  input=file, output=fileout, options='-f nc')
    ds = xr.open_dataset(fileout)
    savepath = f'{path}BottomT/StGrid/stGrid_{base_filename}'
    ds_cropped = ds.sel(longitude=slice(260, 320), latitude=slice(20, 70))
    ds_cropped.to_netcdf(savepath)
    os.remove(fileout)
    print(f'Completed regridding {base_filename} out of {str(len(dsThetao))}')

dsSO = glob.glob(f'{path}RawTmpFiles/so*_ssp585*')
dsSO = glob.glob(f'{path}RawTmpFiles/so*_historical*')
file = "/Users/mdzaugis/Box/RES_Data/CMIP6/RawTmpFiles/so_NorESM2-MM_r1i1p1f1_historical.nc"
for file in dsSO:
    gridfi = f'{path}GridFiles/SODA_grid.nc'
    base_filename = os.path.basename(file)
    fileout = f'{path}BottomSal/tmpfiles/stGrid_{base_filename}'
    cdo.remapdis(gridfi,  input=file, output=fileout, options='-f nc')
    ds = xr.open_dataset(fileout)
    savepath = f'{path}BottomSal/StGrid/stGrid_{base_filename}'
    ds_cropped = ds.sel(longitude=slice(260, 320), latitude=slice(20, 70))
    ds_cropped.to_netcdf(savepath)
    os.remove(fileout)
    print(f'Completed regridding {base_filename} out of {str(len(dsSO))}')

dsSurSO = glob.glob(f'{path}RawTmpFiles/Surface*_ssp585*')
dsSurSO = glob.glob(f'{path}RawTmpFiles/Surface*_historical*')
for file in dsSurSO:
    gridfi = f'{path}GridFiles/SODA_grid.nc'
    base_filename = os.path.basename(file)
    fileout = f'{path}SurSalinity/tmpfiles/stGrid_{base_filename}'
    cdo.remapdis(gridfi,  input=file, output=fileout, options='-f nc')
    ds = xr.open_dataset(fileout, decode_times=True)
    savepath = f'{path}SurSalinity/StGrid/stGrid_{base_filename}'
    ds_cropped = ds.sel(longitude=slice(260, 320), latitude=slice(20, 70))
    ds_cropped.to_netcdf(savepath)
    os.remove(fileout)
    print(f'Completed regridding {base_filename} out of {str(len(dsSurSO))}')

