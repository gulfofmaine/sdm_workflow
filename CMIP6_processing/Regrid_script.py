
import cdo                       # for the cdo tools
import glob                      # for file listing
import xarray as xr              # for netcdfs
import numpy as np               # array tools
import fcts
import os
import re
import shutil

cdo = cdo.Cdo()                      # need to initialize a cdo object

print('Enter User Name')
UsrName = input()
print('Enter Dest Group (i.e. RES_Data)')
Group = input()
print('Enter Folder with a / (i.e. CMIP6/)')
Folder = input()
print('Enter Model Name (i.e. ACCESS-CM2)')
ModelName = input()
print('Enter Model Scenario (i.e. ssp585)')
Scenario = input()

path = fcts.shared_path(user_name=UsrName, group=Group, folder=Folder)
path_BT = os.path.join(path, "BottomT/tmpOriginalGrid/")
path_AllDepths = path + "RawTmpFiles/"+"thetao*"+ModelName+"*.nc"
ncAllBTFiles = glob.glob(path_AllDepths)

for filename in ncAllBTFiles:
    base_filename = os.path.basename(filename)
    temp_ds = xr.open_dataset(fileName)
    bottom_400 = temp_ds.isel(lev=slice(0,18))
    temp_array = bottom_400['thetao']
    depth_indices = fcts.find_deepest_depth_indices(bottom_400)
    ind = xr.DataArray(depth_indices, dims = ['j', 'i'])
    bottomTemps = temp_array.isel(lev=ind)
    ds = bottomTemps.to_dataset()
    ds = ds.rename({'lev': 'bottom'}).reset_coords('bottom')
    save_file = os.path.join(path_BT, base_filename)
    ds.to_netcdf(save_file)
    os.remove(fileName)
    print('Completed finding BT '+ str(f+1) + ' out of '+str(len(ncAllBTFiles)))

ncBTFiles = glob.glob(os.path.join(path, "BottomT/tmpOriginalGrid/thetao*.nc"))

for filename in ncBTFiles:
    gridfi = os.path.join(path, "TestFiles/tempbot400_Omon_CanESM5_ssp585_r1i1p2f1_gn_202001-209912.nc.1x1.nc")
    base_filename = os.path.basename(ilename)
    filein = ilename
    fileout = "".join([path, "BottomT/tmpStd1x1grid/"'1x1_', base_filename])
    cdo.remapdis(gridfi,  input=filein, output=fileout, options='-f nc')
    os.remove(filein)
    print('Completed regridding ' + str(ff+1) + ' out of ' + str(len(ncBTFiles)))

modPath = 'BottomT/tmpStd1x1grid/*thetao*'+ModelName+'*.nc'

ncBTstdGrd = glob.glob(os.path.join(path, modPath))
base_filename = os.path.basename(ncBTstdGrd[0])[0:-17]

base_filename = fcts.nameFile(base_filename, Scenario)

outoutPath = path+'BottomT/Std1x1grid/'+base_filename+'_histoircal_'+Scenario+'.nc'

cdo.cat(input=ncBTstdGrd,
        options='-r',
        output=outoutPath)

for file in ncBTstdGrd:
	os.remove(file)


#SST regrid

path_SST = path + "SST/RawTmpFiles/tos*"+ModelName+"*.nc"
ncTOSFiles = glob.glob(path_SST)

for filename in ncTOSFiles:
    gridfi = os.path.join(path, "TestFiles/tos_Omon_CanESM5_historical_r1i1p2f1_gn_195501-201412.nc.1x1.nc")
    base_filename = os.path.basename(filename)
    filein = filename
    fileout = "".join([path, "SST/tmpStd1x1grid/"'1x1_', base_filename])
    cdo.remapdis(gridfi,  input=filein, output=fileout, options='-f nc')
    os.remove(filein)
    print('Completed regridding ' + str(ff+1) + ' out of ' + str(len(ncTOSFiles)))

grdpath = path + "SST/tmpStd1x1grid/1x1_tos*"+ModelName+"*.nc"
grdfiles = glob.glob(grdpath)
base_filename = os.path.basename(grdfiles[0])[0:-17]
base_filename = fcts.nameFile(base_filename, Scenario)

outoutPath = path+'SST/Std1x1grid/'+base_filename+'_histoircal_'+Scenario+'.nc'

sstfilesin = glob.glob(path+'SST/tmpStd1x1grid/*'+ModelName+'*.nc')

cdo.cat(input=sstfilesin,
        options='-r',
        output=outoutPath)

for file in sstfilesin:
	os.remove(file)
