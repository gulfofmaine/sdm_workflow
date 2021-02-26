from matplotlib import pyplot as plt
import numpy as np
import pandas as pd
import xarray as xr
import zarr
import gcsfs
from dask.diagnostics import ProgressBar
from numpy import unique
import cf_xarray
import operator
import fcts
import glob
import os

# models we care about
source_list = ['IPSL-CM6A-LR',
'CMCC-CM2-SR5',
'MIROC6',
'CanESM5',
'MRI-ESM2-0',
'HadGEM3-GC31-LL',
'GFDL-CM4',
'CESM2',
'CESM2-WACCM',
'CIESM',
'CNRM-CM6-1',
'CNRM-ESM2-1',
'CanESM5-CanOE',
'EC-Earth3',
'EC-Earth3-Veg',
'EC-Earth3-Veg-LR',
'FGOALS-g3',
'FGOALS-f3-L',
'FIO-ESM-2-0',
'GISS-E2-1-G',
'INM-CM4-8',
'INM-CM5-0',
'MIROC-ES2L',
'MRI-ESM2-0',
'NESM3',
'NorESM2-LM',
'NorESM2-MM',
'UKESM1-0-LL']

# enter the var of interest
variable_id = 'thetao'

# enter the table (based on the frequency of measurements)
table_id = 'Omon'

# Enter the experiments of interest
filter_list = ['historical', 'ssp585']


grp1 = 'source_id' # used for grouping normally don't need to change
grp2 = 'member_id' # used for grouping normally don't need to change


# save path
UsrName='mdzaugis'
Group='RES_Data'
Folder='CMIP6/'

path = fcts.shared_path(user_name=UsrName, group=Group, folder=Folder)
# Browse Catalog

# data catalog is stored as a 30MB CSV file

AllModels = pd.read_csv('https://storage.googleapis.com/cmip6/cmip6-zarr-consolidated-stores.csv')

# the columns correspond to the CMIP6 controlled vocab
#Load data

# To access an individual run
df = AllModels.query(f"source_id == 'FIO-ESM-2-0' & variable_id == 'so' & experiment_id == 'historical' & member_id == 'r1i1p1f1' & table_id == 'Omon'")
filteredModels_grid = df.reset_index(drop=True)

df_var = AllModels.query(f"variable_id == '{variable_id}' & table_id == '{table_id}' & experiment_id == @filter_list")
filteredModels = fcts.ExperimentFilter(df_var, grp1, grp2)
filteredModels_grid = filteredModels.query(f"source_id == @source_list").reset_index(drop=True)

# Code below for salinity (surface and bottom) and bottom temperature

# SST is at the bottom of the page

TOP = False  # True if looking for surface, False for bottom

# Only has to be defined once
gcs = gcsfs.GCSFileSystem(token='anon')
#i=0
for i in range(len(filteredModels_grid)):
    source_id = filteredModels_grid.source_id[i]
    member_id = filteredModels_grid.member_id[i]
    experiment_id = filteredModels_grid.experiment_id[i]
    variable_id = filteredModels_grid.variable_id[i]
    if TOP == True:
        savePath = f'{path}RawTmpFiles/Surface_{variable_id}_{source_id}_{member_id}_{experiment_id}.nc'
    else:
        savePath = f'{path}RawTmpFiles/{variable_id}_{source_id}_{member_id}_{experiment_id}.nc'

    # get the path to a specific zarr store 0 index is first on list
    zstore = filteredModels_grid.zstore.values[i]

    # create a mutable-mapping-styly interface to the store
    mapper = gcs.get_mapper(zstore)

    # open it using xarray and zarr
    ds = xr.open_zarr(mapper, consolidated=True)

    lonNames = list(ds.cf[['longitude']].coords)
    latNames = list(ds.cf[['latitude']].coords)

    try:
        vertNames = list(ds.cf[['vertical']].coords)
    except KeyError:
        vertNames = list(ds.cf[['Z']].coords)

    lons = ['lon', 'longitude', 'nav_lon']
    lats = ['lat', 'latitude', 'nav_lat']
    verts = ['lev', 'olevel']

    x_coord = list(set(lonNames).intersection(lons))[0]
    y_coord = list(set(latNames).intersection(lats))[0]
    depth_coord = list(set(vertNames).intersection(verts))[0]

    if len(ds[variable_id][x_coord].dims) == 2:
        multiIndex = True
    else:
        multiIndex = False

    try:
        levUnits = ds[depth_coord].units
    except AttributeError:
        print('No depth units')
        print(ds[depth_coord])
        print('Enter units')
        levUnits = input()

    if levUnits in ['m', 'meters']:
        maxDepth = 400
    elif levUnits in ['cm', 'centimeters']:
        maxDepth = 400 * 100
    else:
        ds[depth_coord]
        print('Check attributes')

   # download atlantic data

    x_coordMin = ds[x_coord].values.min()
    x_coordMax = ds[x_coord].values.max()

    if x_coordMin < 0:
        xmin = -100
        xmax = -40
    else:
        xmin = 260
        xmax = 320

    kwlon = {x_coord: slice(xmin, xmax)}
    kwlat = {y_coord: slice(20, 70)}

    if multiIndex == True:
        # for multi index
        atlantic = ds.where((xmin < ds[x_coord]) & (ds[x_coord] < xmax)
                            & (20 < ds[y_coord]) & (ds[y_coord] < 70),
                            drop=True)
        if experiment_id == 'historical':
            atlantic = atlantic.sel(time=slice('1950-01-01', None))
        elif experiment_id == 'ssp585':
            atlantic = atlantic.isel(time=slice(None, 1032))
        else:
            print("Need to enter date range")
    else:
        # single index
        if experiment_id == 'historical':
            atlantic = ds.sel(**kwlon, **kwlat, time=slice("1950-01-01", None))
        elif experiment_id == 'ssp585':
            atlantic = ds.sel(**kwlon, **kwlat, time=slice(None, '2100-12-31'))
        else:
            print("Need to enter date range")


    if TOP == True:
        kwargs = {depth_coord: 0}
        ds = atlantic.isel(**kwargs)
        ds = ds.rename({depth_coord: 'surface'})

    else:

        kwargs = {depth_coord: slice(0, maxDepth)}
        bottom_400 = atlantic.sel(**kwargs)
        temp_array = bottom_400[variable_id]

        if multiIndex == True:
            dims0 = bottom_400[y_coord].dims[0]
            dims1 = bottom_400[y_coord].dims[1]
        else:
            dims0 = y_coord
            dims1 = x_coord

        depth_indices = fcts.find_deepest_depth_indices_CMIP6(bottom_400, dims0, dims1, variable_id, y_coord, x_coord)
        ind = xr.DataArray(depth_indices, dims=[dims0, dims1])

        kwdepth = {depth_coord: ind}
        dsSel = temp_array.isel(**kwdepth)
        ds = dsSel.to_dataset()
        ds = ds.rename({depth_coord: 'bottom'})

    delayed_obj = ds.to_netcdf(savePath, compute=False)

    with ProgressBar():
        results = delayed_obj.compute()

    print(f'Finished {variable_id}_{source_id}_{member_id}_{experiment_id}.nc')


# Sea Surface Temperature
if variable_id == 'tos':
    source_id = filteredModels_grid.source_id[i]
    member_id = filteredModels_grid.member_id[i]
    experiment_id = filteredModels_grid.experiment_id[i]
    variable_id = filteredModels_grid.variable_id[i]
    savePath = f'{path}RawTmpFiles/{variable_id}_{source_id}_{member_id}_{experiment_id}.nc'

    # get the path to a specific zarr store 0 index is first on list
    zstore = filteredModels_grid.zstore.values[i]

    # create a mutable-mapping-styly interface to the store
    mapper = gcs.get_mapper(zstore)

    # open it using xarray and zarr
    ds = xr.open_zarr(mapper, consolidated=True)

    lons = ['lon', 'longitude', 'nav_lon']
    lats = ['lat', 'latitude', 'nav_lat']

    lonNames = list(ds.cf[['longitude']].coords)
    latNames = list(ds.cf[['latitude']].coords)

    x_coord = list(set(lonNames).intersection(lons))[0]
    y_coord = list(set(latNames).intersection(lats))[0]

    x_coordMin = ds[x_coord].values.min()
    x_coordMax = ds[x_coord].values.max()


    if len(ds[variable_id][x_coord].dims) == 2:
        multiIndex = True
    else:
        multiIndex = False


    if x_coordMin < 0:
        xmin = -100
        xmax = -40
    else:
        xmin = 260
        xmax = 320

    kwlon = {x_coord: slice(xmin, xmax)}
    kwlat = {y_coord: slice(20, 70)}

    if multiIndex == True:
        # for multi index
        atlantic = ds.where((xmin < ds[x_coord]) & (ds[x_coord] < xmax)
                            & (20 < ds[y_coord]) & (ds[y_coord] < 70),
                            drop=True)
        if experiment_id == 'historical':
            atlantic = atlantic.sel(time=slice('1950-01-01', None))
        elif experiment_id == 'ssp585':
            atlantic = atlantic.isel(time=slice(None, 1032))
        else:
            print("Need to enter date range")
    else:
        # single index
        if experiment_id == 'historical':
            atlantic = ds.sel(**kwlon, **kwlat, time=slice("1950-01-01", None))
        elif experiment_id == 'ssp585':
            atlantic = ds.sel(**kwlon, **kwlat, time=slice(None, '2100-12-31'))
        else:
            print("Need to enter date range")

    delayed_obj = atlantic.to_netcdf(savePath, compute=False)

    with ProgressBar():
        results = delayed_obj.compute()




names = {'name': [], 'minDate': [], 'maxData': [], 'length': []}
ncTimes = pd.DataFrame(data=names)

folder = glob.glob(f'{path}SurSalinity/StGrid/*')
for file in folder:
    df = fcts.checkDates(file)
    ncTimes = ncTimes.append(df, ignore_index=True)
