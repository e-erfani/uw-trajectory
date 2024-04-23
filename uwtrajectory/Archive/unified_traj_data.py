#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#############################################################
##### Modified by Ehsan Erfani     @UW on 17 Jan. 2022  #####
##### Created by Johannes Mohrmann @UW on 20 Jul. 2020  #####
#############################################################

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

all_datasets = ['ERA', 'ERA_sfc', 'ERA_ens', 'MODIS_pbl', 'MERRA', 'SSMI', 'CERES', 'AMSR_Tb', 'AMSR']

def rounder(values):
    def f(x):
        idx = np.argmin(np.abs(values - x))
        return values[idx]
    return np.frompyfunc(f, 1, 1)

def xarray_from_cset_flight_trajectory(rfnum, trajnum, trajectory_type='500m_+72'):
    tdump = utils.load_flight_trajectory(rfnum, trajnum, trajectory_type=trajectory_type)
    ds = xarray_from_tdump(tdump)
    global_attrs = [{'CSET_flight': rfnum},
        {'flight_trajectory': str(trajnum)}]
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
        {'Title': "CSET Unified Trajectory Product"},
        {'institution': "Department of Atmospheric Sciences, University of Washington"},
        {'contact': "jkcm@uw.edu"},
        {'Creation Time': str(dt.datetime.utcnow())},
        {'trajectory_setup': "Trajectories were run isobarically " +
                            "from an initialization height of 500m " +
                            "for 72 hours, using GDAS analysis met data"},
        {'HYSPLIT_params': "Trajectories run using HYSPLIT (Hybrid Single "+
                   "Particle Lagrangian Integrated Trajectory Model). "+
                   "Acknowledgements to the NOAA Air Resources Laboratory "+
                   "(ARL) for the provision of the HYSPLIT transport and "+
                   "dispersion model used in this publication."},
        {'HYSPLIT_reference': "Stein, A.F., Draxler, R.R, Rolph, G.D., Stunder, "+
                      "B.J.B., Cohen, M.D., and Ngan, F., (2015). NOAA's "+
                      "HYSPLIT atmospheric transport and dispersion modeling "+ 
                      "system, Bull. Amer. Meteor. Soc., 96, 2059-2077, "+
                      "http://dx.doi.org/10.1175/BAMS-D-14-00110.1"}]
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



def make_trajectory(ds, skip=[], save=False):
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
    if not 'MODIS_pbl' in skip:
        print("adding MODIS...")
        ds = add_MODIS_pbl_to_trajectory(ds)
    if not 'MERRA' in skip:
        print("adding MERRA...")
        ds = add_MERRA_to_trajectory(ds)
    if not 'SSMI' in skip:
        print("adding SSMI...")
        ds = add_SSMI_to_trajectory(ds)
    if not 'CERES' in skip:
        print("adding CERES...")
        ds = add_CERES_to_trajectory(ds)
    if not 'AMSR_Tb' in skip:
        print("adding AMSR Tb...")
        ds = add_AMSR_Tb_to_trajectory(ds)
    if not 'AMSR' in skip:
        print("adding AMSR2...")
        ds = add_AMSR_to_trajectory(ds)
    if save:
        save_trajectory_to_netcdf(ds, save)
    return ds

def make_CSET_trajectory(rfnum, trajnum, save=False, trajectory_type='500m_+72', skip=[]):
    ds = xarray_from_cset_flight_trajectory(rfnum, trajnum, trajectory_type)
    ds = make_trajectory(ds, skip=skip, save=save)
    return ds
    
if __name__ == "__main__":

    
    force_override = True
#    for case_num, case in utils.all_cases.items():
    for case_num, case in utils.all_cases[4].items():
        print('working on case {}'.format(case_num))
#         if case_num not in [6, 10]:
#             continue
        flight = case['TLC_name'].split("_")[1][:4].lower()
        traj_list = case['TLC_name'].split('_')[2].split('-')
        for dirn in ['forward', 'backward']:
            nc_dirstring = '48h_backward' if dirn == 'backward' else '72h_forward'
            for traj in traj_list:
#                 if traj not in ['2.3', '6.0']:
#                    continue
                name = os.path.join(config.trajectory_netcdf_dir, "{}_{}_{}.nc".format(flight, nc_dirstring, traj))
                print("working on {}...".format(os.path.basename(name)))
                if os.path.exists(name):
                    print("already exists!")
                    if not force_override:
                        continue
                    else:
                        print('overriding')
                        os.rename(name, os.path.join(config.trajectory_netcdf_dir, 'old', "{}_{}_{}.nc".format(flight, nc_dirstring, traj)))
    #             ds = make_CSET_trajectory(rfnum=flight, trajnum=float(traj), save=name);
                trajectory_type = '500m_-48' if dirn == 'backward' else '500m_+72'
                print(name)
                ds = make_CSET_trajectory(rfnum=flight, trajnum=float(traj), save=name, trajectory_type=trajectory_type, skip=['ERA_ens']);


    #ds = add_ERA_sfc_data(ds)
    #ds = make_CSET_trajectory(rfnum='rf06', trajnum=2.3, save=False)
    #save_trajectory_to_netcdf(ds, r'/home/disk/eos4/jkcm/Data/CSET/model_forcings/rf06_traj_2.3_fullcolumn_withz.nc')

#     all_trajs = {'rf06': [1.6, 2.0, 2.3, 2.6, 3.0],
#                  'rf10': [5.5, 6.0]}


    # for flight, traj_list in all_trajs.items():
    #     for traj in traj_list:
    #         name = os.path.join(config.trajectory_netcdf_dir, "{}_MODIS_traj_{:0.1f}.nc".format(flight, traj))
    #         print("working on {}...".format(os.path.basename(name)))
    #         ds = make_CSET_trajectory(rfnum=flight, trajnum=traj, save=name);


    # ds = make_CSET_trajectory(rfnum='rf06', trajnum=2.3, save=False)
    # save_trajectory_to_netcdf(ds, r'/home/disk/eos4/jkcm/Data/CSET/Lagrangian_project/trajectory_files/rf06_MODIS_traj_2.3.nc')


    # ds = make_CSET_trajectory(rfnum='rf10', trajnum=6.0, save=False)
    # save_trajectory_to_netcdf(ds, r'/home/disk/eos4/jkcm/Data/CSET/Lagrangian_project/trajectory_files/rf10_MODIS_traj_6.0.nc')

