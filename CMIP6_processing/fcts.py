import numpy as np

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


def find_deepest_depth_indices(ds):
    # First get the vertical True/False of valid values
    idx = ds["thetao"].isel(time = 0).isnull()
    idx_vals = idx.values
    # Create the initial final array to store indices (integer type)
    depth_indices = np.zeros((len(idx.j), len(idx.i))).astype(int)

    # Now find the deepest depth where values are True and store in indices array
    for i in range(len(ds.i.values)):
        for j in range(len(ds.j.values)):
            located = np.where(idx_vals[:, j, i] == False)
            try:
                depth_indices[j, i] = int(located[-1][-1])
            except IndexError:
                depth_indices[j, i] = 1

    return depth_indices

# Eventually try to remove the nested for loop using .map
# dask (with xarray) has implicit parallelism
# numba can take a for loop and create a parallel process with @numba.jit


def nameFile(base_filename, Scenario):
    import re
    if re.search('historical', base_filename):
        base_filename = re.sub('historical_', '', base_filename)

    elif re.search(Scenario, base_filename):

        base_filename = re.sub(Scenario, '', base_filename)
    else:
        print("not recognized")

    return base_filename