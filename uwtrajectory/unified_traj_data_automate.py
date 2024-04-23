'''
Main Code to compile a trajectory
By Ehsan Erfani @UW, 2022
Original Python code by Johannes Mohrmann @UW, summer 2020
Summary of Modifications by Ehsan:
Changes in parameters and functions to facilitate automation and to generalize the code for projects other than CSET
Addition of optional part 2 (visualization and verification) (This part is in the jupyter notebook)
Note by Ehsan: Some names of functions or variables in this code or other codes of the package might refer to "CSET", but most should work for any other project.
'''

### Importing required libraries
import datetime as dt
import numpy as np
import os
import xarray as xr
import matplotlib as mpl
import matplotlib.pyplot as plt
import glob
import pandas as pd
from itertools import cycle
from geographiclib.geodesic import Geodesic
import time
import matplotlib
from pandas import Series, DataFrame
from netCDF4 import Dataset
from math import log10
from scipy import interpolate
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import warnings
warnings.filterwarnings("ignore")


from MERRA2.add_to_trajectory import add_MERRA_to_trajectory
from ERA5.add_to_trajectory import add_ERA_ens_to_trajectory, \
     add_ERA_sfc_to_trajectory, add_ERA_to_trajectory, add_advection_to_trajectory
from AMSR_Tb.add_to_trajectory import add_AMSR_Tb_to_trajectory
from CERES.add_to_trajectory import add_CERES_to_trajectory
from MODIS_pbl.add_to_trajectory import add_MODIS_pbl_to_trajectory
from SSMI.add_to_trajectory import add_SSMI_to_trajectory
from AMSR.add_to_trajectory import add_AMSR_to_trajectory

import utils
import config
import met_utils
import les_utils


##### Parameters:
all_datasets = ['ERA', 'ERA_sfc', 'ERA_ens', 'MERRA', 'CERES', 'SSMI', 'AMSR', 'AMSR_Tb', 'MODIS_pbl']
skip1        = ['ERA_ens']  ## Which of the above datasets should be skipped?
traj_dirc    = ['forward']  ## trajectory direction, options: 'forward', 'backward'
AUTHORS      = ['Ehsan Erfani (Ehsan@nevada.unr.edu); Rob Wood (robwood2@uw.edu)']

# start_dates   = '201806010900' # '201904231900' # start date of trajectory set in the code "unified_traj_data", format: YYYYMMDDHHMM

## Parameters for iterative runs:
YYYY = ['2018', '2019', '2020', '2021']
MM   = ['06', '07', '08']
DD   = [30  , 31  , 31]

traj_locs     = ['01', '02', '03', '04', '05', '06'] 

traj_loc_longs= ['(37.7 N, 131.4 W): A location towards the northern edge of the NEP Sc deck', \
                 '(35.0 N, 125.0 W): GPCI-S12 location in Lewis (2016)', \
                 '(32.0 N, 129.0 W): GPCI-S11 location in Lewis (2016)', \
                 '(29.0 N, 133.0 W): GPCI-S10 location in Lewis (2016)', \
                 '(26.0 N, 137.0 W): GPCI-S09 location in Lewis (2016)', \
                 '(25.0 N, 125.0 W): NEP location in Sandu et al. (2010)']


'''
Main part: Create the Trajectory:
Required functions
'''

def rounder(values):
    def f(x):
        idx = np.argmin(np.abs(values - x))
        return values[idx]
    return np.frompyfunc(f, 1, 1)    

def xarray_from_cset_flight_trajectory(rfnum, trajnum, trajectory_type='500m_+72'):
    tdump = utils.load_flight_trajectory(traj_loc, start_dates, rfnum, trajnum, trajectory_type=trajectory_type)
    ds = xarray_from_tdump(tdump)
    global_attrs = [{'flight': rfnum},
        {'trajectory': str(trajnum)}]
    for i in global_attrs:  # note: an OrderedDict would be tidier, but does not unpack in order
        ds = ds.assign_attrs(**i)
    return ds
    
def xarray_from_tdumpfile(tdumpfile):
    tdump = utils.read_tdump(tdumpfile).sort_values('dtime')
    ds = xarray_from_tdump(tdump)
    return ds
    

def xarray_from_tdump(tdump):
    ds = xr.Dataset.from_dataframe(tdump).drop(['tnum', 'gnum', 'age'])
    ds = ds.rename({'dtime': 'time'})
    # assigning global attributes
    global_attrs = [
        {'Title': "Unified Trajectory Product"},
        {'institution': "Department of Atmospheric Sciences, University of Washington"},
        {'Authors': AUTHORS},
        {'Creation_Time': str(dt.datetime.utcnow())},
        {'Start_Time_of_Trajectory': start_dates},
        {'Start_Location_of_Trajectory': traj_loc_long},
        {'trajectory_setup': "Trajectories were run isobarically " +
                            "from an initialization height of 500m " +
                            "for 72 hours, using ECMWF ERA5 reanalysis met data"},
        {'Trajectory_model_params': "The trajectories were produced by Ryan Eastman @ UW by using ECMWF ERA5 "+
                   "Lagrangian Trajectory Model using horizontal winds at 950 hPa (within the PBL). "},
        {'Trajectory_model_reference': "Eastman, R., & Wood, R. (2016). Factors controlling low-cloud "+
                      "evolution over the eastern subtropical oceans: A Lagrangian perspective using "+
                      "the A-Train satellites. Journal of Atmospheric Sciences, 73(1), 331-351. "+ 
                      "https://doi.org/10.1175/JAS-D-15-0193.1"}]
    for i in global_attrs:  # note: an OrderedDict would be tidier, but does not unpack in order
        ds = ds.assign_attrs(**i)
    
    # assigning variable attributes
    var_attrs = {
        'lon': {'long_name': 'longitude', 
                'units': 'degrees N'},
        'lat': {'long_name': 'latitude',
                'units': 'degrees E'},
        'fhour': {'long_name': 'forecast_lead_time',
                  'units': 'hours'},
        'pres': {'long_name':'trajectory_pressure',
                 'units': 'hPa'},
        'height': {'long_name': 'trajectory_height_above_ground',
                  'units': 'meters'}}
    for k,v in var_attrs.items():
        ds[k] = ds[k].assign_attrs(**v)
    ds.time.attrs['long_name'] = 'time'
    return ds

def save_trajectory_to_netcdf(ds, location):
    ds.to_netcdf(location)

def add_speeds_to_trajectories(ds):
    """Add speed variables to trajectory. used centered difference of distances traveled
    """
    lats, lons, times = ds.lat.values, ds.lon.values, ds.time.values
    
    heading_starts, heading_ends, seg_speeds = [], [], []
    
    for i in range(len(lats)-1):
        geod = Geodesic.WGS84.Inverse(lats[i], lons[i], lats[i+1], lons[i+1])
        dtime = (times[i+1]-times[i])/np.timedelta64(1, 's')
        heading_starts.append(geod['azi1'])
        heading_ends.append(geod['azi2'])
        seg_speeds.append(geod['s12']/dtime)

    #speeds are centered difference, except at start and end, where they are speeds of 
    #first and last trajectory segments
    #headings are average of end azimuth of previous segment/start azimuth of next geodesic segment,
    #except at start and end, where are just the start/end azimuths of the first/last geodesic
    speeds = np.mean(np.vstack([seg_speeds+[seg_speeds[-1]],[seg_speeds[0]]+seg_speeds]), axis=0)
     #headings = np.mean(np.vstack([[heading_starts[0]]+heading_ends, heading_starts+[heading_ends[-1]]]), axis=0) THIS HAD A BUG
    def radial_mean(h1, h2):
        diff = ((h2-h1)+180)%360-180
        return h1 + diff/2
    headings = radial_mean(np.array([heading_starts[0]]+heading_ends), np.array(heading_starts+[heading_ends[-1]]))
    
    u = speeds*np.cos(np.deg2rad(90-headings))
    v = speeds*np.sin(np.deg2rad(90-headings))
    
    ds['traj_u'] = (('time'), u, {'long_name': 'U component of trajectory velocity', 'units': "m s**-1"})
    ds['traj_v'] = (('time'), v, {'long_name': 'V component of trajectory velocity', 'units': "m s**-1"})
    ds['traj_hdg'] = (('time'), headings, {'long_name': 'Trajectory heading', 'units': 'deg'})
    ds['traj_spd'] = (('time'), speeds, {'long_name': 'Trajectory speed', 'units': "m s**-1"})
    return ds   



def make_trajectory(ds, skip=skip1, save=False):
    ds = add_speeds_to_trajectories(ds)        
    if not 'ERA' in skip:
        print("adding ERA...")
        ds = add_ERA_to_trajectory(ds)
        print('adding advection...')
        ds = add_advection_to_trajectory(ds)
    if not 'ERA_sfc' in skip:
        print('adding ERA sfc data...')
        ds = add_ERA_sfc_to_trajectory(ds)
    if not 'ERA_ens' in skip:
        print('adding ERA ensemble data...')
        ds = add_ERA_ens_to_trajectory(ds)        
    if not 'MERRA' in skip:
        print("adding MERRA...")
        ds = add_MERRA_to_trajectory(ds)
    if not 'CERES' in skip:
        print("adding CERES...")
        ds = add_CERES_to_trajectory(ds)        
    if not 'SSMI' in skip:
        print("adding SSMI...")
        ds = add_SSMI_to_trajectory(ds)
    if not 'AMSR' in skip:
        print("adding AMSR2...")
        ds = add_AMSR_to_trajectory(ds)
    if not 'AMSR_Tb' in skip:
        print("adding AMSR Tb...")
        ds = add_AMSR_Tb_to_trajectory(ds)
    if not 'MODIS_pbl' in skip:
        print("adding MODIS...")
        ds = add_MODIS_pbl_to_trajectory(ds)
    if save:
        save_trajectory_to_netcdf(ds, save)
    return ds

## This function is significantly changed by Ehsan Erfani and is now generalized for future trajectories
def make_CSET_trajectory(rfnum, trajnum, save=False, trajectory_type='500m_+72', skip=skip1):
    ds = xarray_from_cset_flight_trajectory(rfnum, trajnum, trajectory_type)
    ds = make_trajectory(ds, skip=skip, save=save)
    return ds

'''
Run the code iteratively
'''

if __name__ == "__main__":

    for y1 in YYYY:
        for jjj, m1 in enumerate(MM):
            for d1 in range(1,DD[jjj]+1):
                for iii in range(len(traj_locs)):    
                    start_dates   = y1 + m1 + f"{d1:02d}" + '0900' # start date of trajectory set in the code "unified_traj_data", format: YYYYMMDDHHMM
                    traj_loc      = traj_locs[iii]
                    traj_loc_long = traj_loc_longs[iii]

                    force_override = True
                    for case_num, case in config.all_cases.items():
                #         print('working on case {}'.format(case_num))
                        flight = case['TLC_name'].split("_")[0][:4].lower()
                        traj_list = case['TLC_name'].split('_')[1].split('-')

                        for dirn in traj_dirc:
                            nc_dirstring = '48h_backward' if dirn == 'backward' else '72h_forward'
                            for traj in traj_list:
                #                 name = os.path.join(config.trajectory_netcdf_dir, "{}_{}_{}.nc".format(flight, nc_dirstring, traj))
                                name = os.path.join(config.trajectory_netcdf_dir, "st_{}_{}_{}.nc"\
                                                    .format(nc_dirstring, traj_loc, start_dates[:-4]))
                                print("working on {}...".format(os.path.basename(name)))
                #                 if os.path.exists(name):
                #                     print("already exists!")
                #                     if not force_override:
                #                         continue
                #                     else:
                #                         print('overriding')
                #                         os.rename(name, os.path.join(config.trajectory_netcdf_dir, 'old', "sr_{}_{}.nc".format(nc_dirstring, start_dates[:-4])))
                                trajectory_type = '500m_-48' if dirn == 'backward' else '500m_+72'
                #                 print(name)
                                ds = make_CSET_trajectory(rfnum=flight, trajnum=float(traj), save=name, \
                                                          trajectory_type=trajectory_type, skip=skip1);






