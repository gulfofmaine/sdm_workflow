import numpy as np
import xarray as xr
import pandas as pd
import gcsfs

# functions
def path_switch_unix(x, folder, user_name):
    default = ["/Users/", user_name, "/Box/Mills Lab/Projects/"]
    paths = {
        'RES_Data': "".join(["/Users/", user_name, "/Box/RES_Data/", folder]),
        'Mills Lab': "".join(["/Users/", user_name, "/Box/Mills Lab/", folder]),
        'Cliamte Change Ecology Lab': "".join(["Users/", user_name, "/Box/Climate Change Ecology Lab/", folder]),
        'NSF OKN': "".join(["/Users/", user_name, "/Box/NSF OKN Demo Data/", folder]),
        'root': "".join(["/Users/", user_name, "/Box/", folder])
    }
    return paths.get(x, "".join(default))


def path_switch_windows(x, folder, user_name):
    default = ["C:/Users/", user_name, "/Box/Mills Lab/Projects/", re.sub(".*\\/", "", os.getcwd()), "/"]
    paths = {
        'RES_Data': "".join(["C:/Users/", user_name, "/Box/Res_Data/", folder]),
        'Mills Lab': "".join(["C:/Users/", user_name, "/Box/Mills Lab/", folder]),
        'Cliamte Change Ecology Lab': "".join(["C:/Users/", user_name, "/Box/Climate Change Ecology Lab/", folder]),
        'NSF OKN': "".join(["C:/Users/", user_name, "/Box/NSF OKN Demo Data/", folder]),
        'root': "".join(["C:/Users/", user_name, "/Box/", folder])
    }
    return paths.get(x, "".join(default))


def shared_path(os_use="unix",
                user_name="",
                group="",
                folder=""):
    if os_use == "unix":

        path_out = path_switch_unix(group, folder, user_name)
    elif os_use == 'windows':

        path_out = path_swith_windows(group, folder, user_name)
    else:
        print("OS not recognized")

    return path_out


def nameFile(base_filename, Scenario):
    import re
    if re.search('historical', base_filename):
        base_filename = re.sub('historical_', '', base_filename)

    elif re.search(Scenario, base_filename):

        base_filename = re.sub(Scenario, '', base_filename)
    else:
        print("not recognized")

    return base_filename

# find bottom temp for CMIP6
def find_deepest_depth_indices_CMIP6(ds, dims0, dims1, variable_id, y_coord, x_coord):
    # First get the vertical True/False of valid values
    idx = ds[variable_id].isel(time=0).isnull()
    idx_vals = idx.values

    # Create the initial final array to store indices (integer type)
    depth_indices = np.zeros((len(idx[y_coord][dims0]), len(idx[x_coord][dims1]))).astype(int)

    # Now find the deepest depth where values are True and store in indices array
    for i in range(len(ds[dims1].values)):
        for j in range(len(ds[dims0].values)):
            located = np.where(idx_vals[:, j, i] == False)
            try:
                depth_indices[j, i] = int(located[-1][-1])
            except IndexError:
                depth_indices[j, i] = 1

    return depth_indices


# find bottom temp for any netcdf with depth
def find_deepest_depth_indices(ds, variable_id, y_coord, x_coord, depth_coord, maxDepth=400):

    kwargs = {depth_coord: slice(0, maxDepth)}
    bottom_400 = ds.sel(**kwargs)

    # First get the vertical True/False of valid values
    idx = bottom_400[variable_id].isel(time=0).isnull()
    idx_vals = idx.values


    if len(bottom_400[variable_id][x_coord].dims) == 2:
        multiIndex = True
    else:
        multiIndex = False

    if multiIndex == True:
        dims0 = bottom_400[y_coord].dims[0]
        dims1 = bottom_400[y_coord].dims[1]
    else:
        dims0 = y_coord
        dims1 = x_coord

    # Create the initial final array to store indices (integer type)
    depth_indices = np.zeros((len(idx[y_coord][dims0]), len(idx[x_coord][dims1]))).astype(int)

    # Now find the deepest depth where values are True and store in indices array
    # numba can take a for loop and create a parallel process with @numba.jit
    for i in range(len(bottom_400[dims1].values)):
        for j in range(len(bottom_400[dims0].values)):
            located = np.where(idx_vals[:, j, i] == False)
            try:
                depth_indices[j, i] = int(located[-1][-1])
            except IndexError:
                depth_indices[j, i] = 1

    ind = xr.DataArray(depth_indices, dims=[dims0, dims1])

    return ind

def countExp(x):
    return len(np.unique(x['experiment_id']))

def ExperimentFilter(df, grp1, grp2):

    df3 = df.groupby([grp1, grp2]).apply(lambda x: countExp(x)).reset_index(name='Number_of_exp')  #.query(f'Number_of_exp == {len(filter_list)}')
    df4 = pd.merge(df3, df, how="left", on=['source_id', 'member_id'])
    df5 = df4.groupby([grp1, 'experiment_id']).apply(lambda x: x.iloc[0]).droplevel(0).reset_index(drop=True)
    return df5

def CheckMeta(dfList):
    meta = []
    for i in range(len(filteredModels)):
        # get the path to a specific zarr store 0 index is first on list
        zstore = filteredModels.zstore.values[i]

        # create a mutable-mapping-styly interface to the store
        mapper = gcs.get_mapper(zstore)

        # open it using xarray and zarr
        ds = xr.open_zarr(mapper, consolidated=True)
        attr = ds.attrs
        meta.append(attr)
    res = list(map(operator.itemgetter('nominal_resolution'), meta))
