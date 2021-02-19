
import numpy as np
import pandas as pd
import xarray as xr
import dask
import cf_xarray
from dask.diagnostics import ProgressBar
import fcts

UsrName = 'mdzaugis'
Group = 'RES_Data'
Folder = 'SODA/'

file = 'SODA_Salt_Red'

path = fcts.shared_path(user_name=UsrName, group=Group, folder=Folder)

ds = xr.open_dataset(f'{path}{File}.nc', chunks={"time": 10})

print(ds.cf.describe())

ind = fcts.find_deepest_depth_indices(ds=ds, variable_id='salt', x_coord='xt_ocean', y_coord='yt_ocean', depth_coord='st_ocean', maxDepth=400)

kwdepth = {'st_ocean': ind}
var_array = ds['salt']
dsSel = var_array.isel(**kwdepth)
ds = dsSel.to_dataset()

savePath = f'{path}{file}_bottomLayer.nc'

delayed_obj = ds.to_netcdf(savePath, compute=False)

with ProgressBar():
    results = delayed_obj.compute()




