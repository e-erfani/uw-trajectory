import os
import sys
sys.path.insert(0, os.path.abspath('..'))

import glob


from uwtrajectory import unified_traj_data as traj



tdumpfiles = glob.glob('/home/disk/eos4/jkcm/Data/HYSPLIT/tdump/cset_profiles/cset_profile*')
for f in tdumpfiles:
    savefile = f'/home/disk/eos4/jkcm/Data/CSET/profile_backtrajectories/profile_backtraj_{os.path.basename(f)[18:]}.nc'
    print(savefile)
    ds = traj.xarray_from_tdumpfile(f)
    ds = traj.make_trajectory(ds, save=savefile) 
#     ax.plot(ds.lon, ds.lat)
#     ax.plot(ds.lon[-1], ds.lat[-1], '*')