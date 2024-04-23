#std library
import glob
import os
import sys
import datetime as dt

#additional libraries
import numpy as np
import xarray as xr

#project libraries
sys.path.insert(0, os.path.abspath('..')) #this would have to be altered to insert ~jkcm/Code/uw-trajectory/ to the path
from uwtrajectory.unified_traj_data import make_trajectory

####################################################################


data_files = glob.glob('/home/disk/eos1/bloss/Runs/EUREC4A/Data/AlongTrajectory/*traj0*nc')
test_file = data_files[0]



ds = xr.open_dataset(test_file)
#parsing in time properly
time = np.array([dt.datetime.strptime(ds.calday.attrs['long_name'], 'Time in days after %HZ on %d %b %Y') 
                 + dt.timedelta(days=i) for i in ds.calday.values])
ds['time'] = time
ds = ds.isel(lat=0, lon=0) # this just removes the useless length-1 dims, which clash with names I need
ds['init_lat'] = ds.lat
ds['init_lon'] = ds.lon
ds = ds.drop(['lat', 'lon'])
ds = ds.rename({'lat_traj': 'lat', 'lon_traj': 'lon'})
ds = make_trajectory(ds, skip=['ERA_ens', 'MODIS_pbl', 'MERRA'], save='/home/disk/eos4/jkcm/testfile.nc')

#this last line will fail because none of the data is downloaded!
