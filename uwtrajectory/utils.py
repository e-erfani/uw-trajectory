# -*- coding: utf-8 -*-
#######################################################
##### Created by Johannes Mohrmann, summer, 2016  #####
##### Modified by Ehsan Erfani (2022):            #####
##### - Now, the start date of a trajectory is    #####
#####   taken from config to help automation      #####
##### - Some changes in names and paths           #####
##### Search my name to see the modifications     #####
#######################################################

import pytz
import os
import re
import pandas as pd
import netCDF4 as nc4
import datetime as dt
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm, rc
from ftplib import FTP
from mpl_toolkits.basemap import Basemap
from time import sleep
from urllib.request import urlopen
from urllib.error import HTTPError
import collections
import matplotlib.path as path
import glob
import xarray as xr
import warnings
import collections
import pickle
import sys
import met_utils as mu
from scipy.interpolate import interp1d
import config
#from . import met_utils as mu
# import xlrd

# %% Parameters


SMALL_SIZE = 16
MEDIUM_SIZE = 20
BIGGER_SIZE = 24

rc('font', size=SMALL_SIZE)          # controls default text sizes
rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
rc('axes', labelsize=SMALL_SIZE)     # fontsize of the x and y labels
rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title
rc('figure', dpi=100)

font = {'family' : 'DejaVu Sans',
        'weight' : 'normal'}
rc('font', **font)


def get_lon_prime(lat, lon, lon0=-140, lat0=30):
        lonp = lon0 + 0.8*(lon-lon0) + 0.4*(lat-lat0)
        return lonp

def nan_helper(y):
    """Helper to handle indices and logical indices of NaNs.

    Input:
        - y, 1d numpy array with possible NaNs
    Output:
        - nans, logical indices of NaNs
        - index, a function, with signature indices= index(logical_indices),
          to convert logical indices of NaNs to 'equivalent' indices
    Example:
        >>> # linear interpolation of NaNs
        >>> nans, x= nan_helper(y)
        >>> y[nans]= np.interp(x(nans), x(~nans), y[~nans])
    """

    return np.isnan(y), lambda z: z.nonzero()[0]
    
    
def gauss2D(shape=(3,3),sigma=0.5):
    """
    2D gaussian mask - should give the same result as MATLAB's
    fspecial('gaussian',[shape],[sigma])
    """
    m,n = [(ss-1.)/2. for ss in shape]
    y,x = np.ogrid[-m:m+1,-n:n+1]
    h = np.exp( -(x*x + y*y) / (2.*sigma*sigma) )
    h[ h < np.finfo(h.dtype).eps*h.max() ] = 0
    sumh = h.sum()
    if sumh != 0:
        h /= sumh
    return h

def closest_index(lat_traj, lon_traj, lat, lon):
    dist = ((lat - lat_traj)**2 + (lon - lon_traj)**2)**(0.5)
    return np.unravel_index(np.nanargmin(dist), dist.shape)


def get_GOES_files_for_dates(date_array):
# if True:
    all_GOES_files = sorted(glob.glob(r'/home/disk/eos4/jkcm/Data/CSET/GOES/VISST_pixel/G15V03.0.NH.*.NC'))
    all_GOES_date_strings = [i[-22:-10] for i in all_GOES_files]
    relevant_dates = [dt.datetime.strftime(i, '%Y%j.%H%M') for i in sorted(as_datetime(date_array))]
    relevant_files = sorted([all_GOES_files[all_GOES_date_strings.index(d)] for d in relevant_dates])
    return relevant_files


def get_ERA_data(var_list, lats, lons, times, pressures, box_degrees=2):
    """Retrieve ERA5 data in a box around a trajectory
    Assumes ERA5 data is 0.3x0.3 degrees
    Returns an xarray Dataset
    """
    space_index = int(np.round(box_degrees/0.3/2)) # go up/down/left/right this many pixels
    assert len(lats) == len(lons) == len(times) == len(pressures)
    unique_days = set([as_datetime(i).date() for i in times])
    files = [os.path.join(ERA_source, "ERA5.pres.NEP.{:%Y-%m-%d}.nc".format(i))
             for i in unique_days]
    return_ds = xr.Dataset(coords={'time': times})
    with xr.open_mfdataset(sorted(files)) as data:
        for var in var_list:
            vals = []
            for (lat, lon, time, pres) in zip(lats, lons%360, times, pressures):
                x = data[var].sel(longitude=slice(lon - box_degrees/2, lon + box_degrees/2),
                                  latitude=slice(lat + box_degrees/2, lat - box_degrees/2))
                y = x.sel(method='nearest', tolerance=np.timedelta64(minutes=59), time=time)
                z = y.sel(method='nearest', tolerance=50, level=pres)
                #this applies a 2D gaussian the width of z, i.e. sigma=box_degrees
                gauss = gauss2D(shape=z.shape, sigma=z.shape[0])
                filtered = z.values * gauss
                vals.append(np.sum(filtered))
            da = xr.DataArray(data=vals, coords={'time': times}, dims=['time'])
            return_ds[var] = da
    return return_ds

lev_map = {'1': 0.0100, '2': 0.0200, '3': 0.0327, '4': 0.0476,
           '5': 0.0660, '6': 0.0893, '7': 0.1197, '8': 0.1595,
           '9': 0.2113, '10': 0.2785, '11': 0.3650, '12': 0.4758,
           '13': 0.6168, '14': 0.7951, '15': 1.0194, '16': 1.3005,
           '17': 1.6508, '18': 2.0850, '19': 2.6202, '20': 3.2764,
           '21': 4.0766, '22': 5.0468, '23': 6.2168, '24': 7.6198,
           '25': 9.2929, '26': 11.2769, '27': 13.6434, '28': 16.4571,
           '29': 19.7916, '30': 23.7304, '31': 28.3678, '32': 33.8100,
           '33': 40.1754, '34': 47.6439, '35': 56.3879, '36': 66.6034,
           '37': 78.5123, '38': 92.3657, '39': 108.6630, '40': 127.8370,
           '41': 150.3930, '42': 176.9300, '43': 208.1520, '44': 244.8750,
           '45': 288.0830, '46': 337.5000, '47': 375.0000, '48': 412.5000,
           '49': 450.0000, '50': 487.5000, '51': 525.0000, '52': 562.5000,
           '53': 600.0000, '54': 637.5000, '55': 675.0000, '56': 700.0000,
           '57': 725.0000, '58': 750.0000, '59': 775.0000, '60': 800.0000,
           '61': 820.0000, '62': 835.0000, '63': 850.0000, '64': 865.0000,
           '65': 880.0000, '66': 895.0000, '67': 910.0000, '68': 925.0000,
           '69': 940.0000, '70': 955.0000, '71': 970.0000, '72': 985.0000}

pres_map = {}
for k, v in lev_map.items():
    pres_map[v] = int(k)

def get_MERRA_level(pressure):
    a, b = zip(*[(float(k), v) for k, v in lev_map.items()])
    levels = sorted(a)
    pressures = sorted(b)

    return(interp1d(pressures, levels)(pressure))

def MERRA_lev(lev, invert=False, lev_map=lev_map):
    if invert:
        pres_map = {}
        for k, v in lev_map.items():
            pres_map[str(v)] = int(k)
        lev_map = pres_map
    if isinstance(lev, collections.Iterable):
        pres = [lev_map[str(int(i))] for i in lev]
    else:
        pres = lev_map[int(float(str(lev)))]
    return pres


def get_MERRA_data(var_list, lats, lons, times, pressures, box_degrees=2):
    """Retrieve ERA5 data in a box around a trajectory
    Assumes ERA5 data is 0.3x0.3 degrees
    Returns an xarray Dataset
    """
    # Merra lat spacing is 0.5 deg (n-s), lon-spacing is 0.625 (e-w)
    #lat_space_index = int(np.round(box_degrees/0.5/2)) # go up/down this many pixels
    #lon_space_index = int(np.round(box_degrees / 0.625 / 2))  # go left-right this many pixels
    assert len(lats) == len(lons) == len(times) == len(pressures)
    unique_days = set([as_datetime(i).date() for i in times])
    files = [os.path.join(MERRA_source, "svc_MERRA2_400.inst3_3d_aer_Nv.{:%Y%m%d}.nc4".format(i))
             for i in unique_days]
    return_ds = xr.Dataset(coords={'time': times})
    with xr.open_mfdataset(sorted(files)) as data:
        for var in var_list:
            vals = []
            for (lat, lon, time, pres) in zip(lats, (lons+180)%360-180, times, pressures):
                x = data[var].sel(lon=slice(lon - box_degrees/2, lon + box_degrees/2),
                                  lat=slice(lat - box_degrees/2, lat + box_degrees/2))
                y = x.sel(method='nearest', tolerance=dt.timedelta(minutes=179), time=time)
                z = y.sel(method='nearest', tolerance=1, lev=get_MERRA_level(pres))
                #this applies a 2D gaussian the width of z, i.e. sigma=box_degrees
                gauss = gauss2D(shape=z.shape, sigma=max(z.shape))
                filtered = z.values * gauss
                vals.append(np.sum(filtered))
            da = xr.DataArray(data=vals, coords={'time': times}, dims=['time'])
            return_ds[var] = da
    return return_ds

def dep_get_MERRA_var(varname, nc, make_vol=True):

    if varname == 'SALT':
        salt_names = ['SS{:03}'.format(i) for i in range(1, 3)]
        var = np.sum([nc[name][:] for name in salt_names], axis=0)
        var_name_str = 'Sea Salt Mixing Ratio (all bins)'
        units = nc['SS001'].units
    elif varname == 'DUST':
        dust_names = ['DU{:03}'.format(i) for i in range(1, 2)]
        var = np.sum([nc[name][:] for name in dust_names], axis=0)
        var_name_str = 'Dust Mixing Ratio (all bins)'
        units = nc['DU001'].units
    elif varname == 'BC':
        BC_names = ['BCPHILIC', 'BCPHOBIC']
        var = np.sum([nc[name][:] for name in BC_names], axis=0)
        var_name_str = 'Black Carbon Mixing Ratio (total)'
        units = nc['BCPHILIC'].units
    elif varname == 'OC':
        OC_names = ['OCPHILIC', 'OCPHOBIC']
        var = np.sum([nc[name][:] for name in OC_names], axis=0)
        var_name_str = 'Organic Carbon Mixing Ratio (total)'
        units = nc['OCPHILIC'].units
    elif varname == 'SG':
        SG_names = ['DMS', 'MSA', 'SO2']
        var = np.sum([nc[name][:] for name in SG_names], axis=0)
        var_name_str = 'Sulfur Compounds Mixing Ratio (total)'
        units = nc['OCPHILIC'].units
    elif varname == 'AEROSOL':
        aa_names = ['SS001', 'SS002', 'DU001', 'BCPHILIC', 'BCPHOBIC',
                    'OCPHILIC', 'OCPHOBIC', 'SO4']
        var = np.sum([nc[name][:] for name in aa_names], axis=0)
        var_name_str = 'Particulate Aerosol Mixing Ratio (total)'
        units = nc['BCPHILIC'].units
    else:
        var = nc[varname][:].squeeze()  # (time, pres, lats, lons)
        var_name_str = (' ').join(nc[varname].long_name.split('_'))
        units = nc[varname].units
    # Sorting out units and variable names, converting to volumetric
    if units == "kg kg-1" and np.mean(var) < 0.01:  # dealing with aerosol
        if make_vol:
            dens = nc['AIRDENS'][:].squeeze()  # (time, pres, lats, lons)
            var = var * dens * 10**9
            units = r' (${\mu}g  m^{-3}$)'
        else:
            var = var*10**9
            units = r' (${\mu}g  kg^{-1}$)'
    else:
        units = (r'\_').join(units.split('_'))
    var_str = var_name_str + ' (' + units + ')'
    return var, var_str

def CSET_date_from_table(date, time):
    """return datetime object from CSET Lookup Table-formatted date and time
    """
    d = as_datetime(dt.datetime.strptime(str(int(date)), '%m%d%y') + dt.timedelta(seconds=time))
    return d


def add_leg_sequence_labels(df, start_times, end_times, legs, sequences):
        """add leg labels to insitu data."""
#         data = self.flight_data
        sequence_array = np.empty(df.time.values.shape, dtype='U1')
        leg_array = np.empty(df.time.values.shape, dtype='U1')
        df['leg'] = (('time'), leg_array)
        df['sequence'] = (('time'), sequence_array)
        for s, e, leg, seq in zip(start_times, end_times, legs, sequences):
            which_times = np.logical_and(as_datetime(df['time'].values) >= s,
                                         as_datetime(df['time'].values) <= e)
            df['leg'][which_times] = leg
            df['sequence'][which_times] = seq
        df = df.set_coords(['leg', 'sequence'])#, inplace=True)
        return df, sequences
#         self.sequences = sorted(list(set(sequences)))

def flightpair_from_flight(flight):
    if isinstance(flight, str):
        if len(flight) == 4:
            flight = int(flight[2:])
        else:
            flight = int(flight)
    if not flight in range(2, 16):
        raise ValueError('invalid flight number')
    if flight % 2 == 0:
        return ('rf{:02d}_rf{:02d}'.format(flight, flight + 1))
    elif flight % 2 == 1:
        return ('rf{:02d}_rf{:02d}'.format(flight - 1, flight))


def get_waypoint_data(flight, waypoint_type='a'):
    # selecting wp file
    flightpair = flightpair_from_flight(flight)
    # This line modified by Ehsan Erfani
    floc = config.traj_txt_dir / '{}_waypoints'.format(waypoint_type)    
    wpfile = os.path.join(floc, flightpair.upper() + '_{}_waypoints.txt'.format(waypoint_type))

    # parsing
    def parseFunc(y, m, d, H, M):
        return dt.datetime(int(y), int(m), int(d), int(H), int(M))

    columns = ['lab', 'outlat', 'outlon', 'out_Y', 'out_M', 'out_D', 'out_HH', 'out_MM',
               'retlat', 'retlon', 'ret_Y', 'ret_M', 'ret_D', 'ret_HH', 'ret_MM']
    if waypoint_type == 'b':
        columns.append('dist')
    data = pd.read_table(wpfile, names=columns, skiprows=3, engine='python',
                         date_parser=parseFunc, index_col='lab',
                         sep='\s+', # delim_whitespace=True, 
                         parse_dates={'out_time': ['out_Y', 'out_M', 'out_D', 'out_HH', 'out_MM'],
                                      'ret_time': ['ret_Y', 'ret_M', 'ret_D', 'ret_HH', 'ret_MM']})
    return (data)


def qv_from_p_T_RH(p, T, RH):
    """p in Pa, T in K, Rh in pct. return is in kg/kg
    """
    es = 611.2*np.exp(17.67*(T-273.15)/(T-29.65))
    qvs = 0.622*es/(p-es)
    qv = qvs*RH/100
    return qv
    
#     qvs = 0.622*es/(p-0.378*es)
#     rvs = qvs/(1-qvs)
#     rv = RH/100. * rvs
#     qv = rv/(1+rv)
#     return qv


def load_flight_trajectory(traj_loc, start_dates, flight, number, trajectory_type='500m_+72'):
    # These three lines modified by Ehsan Erfani
    out_date = start_dates
    trajectory_loc = config.traj_txt_dir # '{}'.format(flightpair)    
    trajectory_name = r'{}.{}.{}.airmass_trajectories_{}.txt'.format(config.traj_prefix, traj_loc, out_date, trajectory_type)
    # date format: {:%Y%m%d%H%M}
    t_data = read_tdump(os.path.join(trajectory_loc, trajectory_name))
    return (t_data.sort_values('dtime'))


def load_flightplan(infile):
    with open(infile, 'rb') as readfile:
        flightplan = pickle.load(readfile)
    return flightplan

def load_flight_file(infile):
    """
    loads a flight file from disk.
    Opposite of make_flight_file
    """
    with open(infile, 'rb') as readfile:
        flightplan = pickle.load(readfile)
    return flightplan

def read_CSET_Lookup_Table(path=None, rf_num=None, sequences=None, legs=None, variables=None):
    """Read in data from the CSET Lookup Table.

    Arguments
    ----------
    path : str
        string representing the location of the lookup table
    rf_num : str or list, int
        list of integers of research flights, or 'all' for all flights
    legs : str or list, str
        list of strings representing the LEG NAMES for which variables should
        be retrieved, or 'all' for all variables
            b: below cloud
            c: in cloud
            a: above cloud
            p: porpoise
            m: Mather sounding
            k: Kona sounding
            f: ferry
            u: up sounding
            d: down sounding
    sequences : str or list, str
        list of strings representing the SEQUENCE NAMES for which variables should
        be retrieved, or 'all' for all defined sequences. first sequence of each flight
        is 'A', last is C-E depending on how many sequences were performed.
        NOTE: 'all' is NOT the same as leaving this as None (default). 'all' will
        explicitly look at all the sequences, so any upper-level data would be excluded.
        None will look only at rf_num and specified legs, ignoring sequences entirely.
    variables: list, str
        list of strings representing variables you want as list. leave blank to
        get error message with all options. Useful ones are 'Date', 'ST', 'ET'
        for date, start time, and end time
    Returns
    ----------
    ret_dict : dict
        dictionary with m+2 entries, where m is the number of requested vars:
            'rf': an array of length n the research flights
            'sequence': an array of length n of the sequences
            for each variable, a dictionary with units and a length n array
            of variable values

    """
    # warnings.warn("NOTE: usage change Feb 2018: sequences now refers to flight sequence (A,B,...) "
    #       "and legs refers to portion of flight ('b', 'p'), etc. see docstring")
    if path is None:
        path = r'/home/disk/eos4/jkcm/Data/CSET/LookupTable_all_flights.xls'
    sheet = xlrd.open_workbook(path).sheet_by_index(0)
    leg_colnum = np.argwhere(np.array(sheet.row_values(11)) == 'Leg Name').flatten()[0]
    all_legs = [str(i) for i in sheet.col_values(leg_colnum)[18:]]
    flight_colnum = np.argwhere(np.array(sheet.row_values(11)) == 'Research Flight Number').flatten()[0]
    all_flights = [int(i) for i in sheet.col_values(flight_colnum)[18:]]
    seq_colnum = np.argwhere(np.array(sheet.row_values(11)) == 'Sequence Name').flatten()[0]
    all_sequences = [str(i) for i in sheet.col_values(seq_colnum)[18:]]
    abb_cell = [str(i.value) for i in sheet.col_slice(0, 0, 10)]
    val_cell = [str(i.value) for i in sheet.col_slice(1, 0, 10)]
    varab = [str(i.value) for i in sheet.row_slice(12, 3, 39)]
    vname = [str(i.value).ljust(28) for i in sheet.row_slice(11, 3, 39)]
    vunit = [str(i.value).ljust(6) for i in sheet.row_slice(16, 3, 39)]

    if legs == 'all':
        legs = [str(i) for i in set(all_legs)]
    elif isinstance(legs, str):
        legs = [legs]
    if rf_num == 'all':
        rf_num = [i for i in set(all_flights)]
    elif isinstance(rf_num, int):
        rf_num = [rf_num]
    if sequences == 'all':
        sequences = [str(i) for i in set(all_sequences)]
    elif isinstance(sequences, str):
        sequences = [sequences]

    # if there is missing input, print some helpful information
    mess = "Missing or incorrect input, printing help"
    if rf_num is None or not set(rf_num) <= set(all_flights):
        mess += ("\n\nspecify the RESEARCH FLIGHTS (rf_num) you want as list."
                 "\noptions are {}".format(str([i for i in set(all_flights)])))
        mess += "\nor select 'all'"
    if legs is None or not set(legs) <= set(all_legs):
        abbs = ['{}: {}'.format(a, b) for (a, b) in zip(abb_cell, val_cell)]
        mess += ("\n\nspecify the LEG NAMES (legs) you want as list.\n"
                 "options are: \n{}".format('\n'.join(abbs)))
        mess += "\nor select 'all'"
    if sequences is not None and not set(sequences) <= set(all_sequences):
        mess += ("\n\neither leave SEQUENCE NAMES (seqs) blank to \n"
                 "ignore sequences, or else specify as list, or select 'all'")
    if variables is None or not set(variables) <= set(varab):
        var = ['{}{}{}'.format(a.ljust(14), b, c) for (a, b, c) in
               zip(varab, vname, vunit)]
        mess += ("\n\nspecify the VARIABLES (variables) you want as list.\n"
                 "options are: \n{}".format('\n'.join(var)))
    if len(mess) > 41:
        raise ValueError(mess)

    # otherwise return the requested values as a dict of dicts\
    rows, = np.where(
        [False]*18 + [True if (l in legs and f in rf_num) else False
                      for l, f in zip(all_legs, all_flights)])
    if sequences is not None:
        seqrows, = np.where(
            [False]*18 + [True if s in sequences else False
                          for s in all_sequences])
        rows = np.intersect1d(rows, seqrows)
    cols, = np.where(
        [False]*3 + [True if v in variables else False for v in varab])
    rf = np.array([int(sheet.cell(r, flight_colnum).value) for r in rows])
    leg = np.array([str(sheet.cell(r, leg_colnum).value) for r in rows])
    seq = np.array([str(sheet.cell(r, seq_colnum).value) for r in rows])

    ret_dict = {'rf': rf, 'leg': leg, 'seq': seq}
    for c in cols:
        varname = str(sheet.cell(12, c).value)
        units = str(sheet.cell(16, c).value)
        values = np.array([sheet.cell(r, c).value for r in rows])
        ret_dict[varname] = {'units': units, 'values': values}
    return ret_dict


def get_leg_times_by_sequence(flightnum, sequence, leg):
    path = r'/home/disk/eos4/jkcm/Data/CSET/LookupTable_all_flights.xls'
    flight = read_CSET_Lookup_Table(path, rf_num=flightnum, sequences=[sequence],
                                       legs=[leg], variables=['Date', 'ST', 'ET'])
    start_times = as_datetime([CSET_date_from_table(d, t) for d, t in
                   zip(flight['Date']['values'], flight['ST']['values'])])
    end_times = as_datetime([CSET_date_from_table(d, t) for d, t in
                 zip(flight['Date']['values'], flight['ET']['values'])])
    sounding_times = list(zip(flight['rf'], start_times, end_times))
    return(sounding_times[0][1], sounding_times[0][2])


def read_CSET_data(fname, var_list=None,
                   start_date=None, end_date=None):
    """read in CSET UHSAS .nc file and returns requested variables
    """
    with nc4.Dataset(fname, 'r') as nc:
        timevar = nc.variables['Time']
        date = nc4.num2date(timevar[:], units=timevar.units)

        if start_date is None:
            start_date = date[0]
        if end_date is None:
            end_date = date[-1]
        indx = np.logical_and(date >= start_date, date <= end_date)
        date = date[indx]
        ret_dict = {'Date': date}
        for var_name in var_list:
            ret_dict[var_name] = nc.variables[var_name][:].squeeze()[indx]
    return ret_dict

def get_GOES_data(variable_list, lat, lon, time, degree, dlat=12, dlon=21):
    def GOES_file_from_date(time, location, max_distance=3):
        """Return the goes filename corresponding to the time in the location folder
        max_distance is the max number of hours away we are allowed to validly look
        dlat is number of lat indexes per degree, same for dlon
        """
#        offs = 0 if time.minute < 30 else 1
        f, b = np.arange(max_distance) + 1, -np.arange(max_distance)
        offs = np.hstack(zip(f,b)) if time.minute > 30 else np.hstack(zip(b,f))
        for off in offs:
            file = "G15V03.0.NH.{:%Y%j.%H}00.PX.08K.NC".format(time + dt.timedelta(hours=int(off)))   
            if os.path.exists(os.path.join(location, file)):
                return os.path.join(location, file)
        raise IOError("no GOES file found!")
        
    file_name = GOES_file_from_date(time, GOES_source)
    with xr.open_dataset(file_name) as data:
#    if True:
#        data = xr.open_dataset(file_name)
        [k for (k,v) in data.coords.items()]
        ret_dict = {}
        lats = data.coords['latitude'].values
        lons = data.coords['longitude'].values
        ilat, ilon = closest_index(lat, lon, lats, lons)
#        dlat = lats[ilat-1,ilon] - lats[ilat,ilon]
#        dlon = lons[ilat,ilon+1] - lons[ilat,ilon]
#        delta_lat = degree/2/dlat
#        delta_lon = degree/2/dlon
#        lat_mask = np.logical_and(lats > lat - degree/2., lats < lat + degree/2.)
#        lon_mask = np.logical_and(lats > lat - degree/2., lats < lat + degree/2.)
#        crd_mask = np.logical_and(lat_mask, lon_mask)
        
        delta_lat = int((degree/2)*dlat)
        delta_lon = int((degree/2)*dlon)
#        print(delta_lat)
#        print(delta_lon)
        latslice = slice(ilat-delta_lat,ilat+delta_lat)
        lonslice = slice(ilon-delta_lon,ilon+delta_lon)
        ret_dict['lat'] = lats[latslice,lonslice]
        ret_dict['lon'] = lons[latslice,lonslice]
        for variable in variable_list:
#            variable = 'visible_count'
            if variable not in data.data_vars.keys():
                raise ValueError("{} variable not in dataset!")
            vardata = data.data_vars[variable].values[latslice,lonslice]
            ret_dict[variable] = vardata
#            ret_dict[variable] = data.data_vars[variable].loc[dict(image_y=latslice,
#                                                                image_x=lonslice)]
        return ret_dict


def get_flight_start_end_times(rf_num, lookup_table_path):
    if rf_num == 16:
        start_time = dt.datetime(2015, 8, 12, 15, 25)
        end_time = dt.datetime(2015, 8, 12, 22, 5)
        return (start_time, end_time)
    x = read_CSET_Lookup_Table(lookup_table_path, rf_num=[rf_num],
                               legs=['m', 'k'], variables=['Date', 'ST', 'ET'])
    if rf_num % 2 == 0:  # westward leg, m is start, k is end
        start_time = CSET_date_from_table(x['Date']['values'][0], x['ST']['values'][0])
        end_time = CSET_date_from_table(x['Date']['values'][1], x['ET']['values'][1])
    else:  # eastward leg
        start_time = CSET_date_from_table(x['Date']['values'][0], x['ST']['values'][0])
        end_time = CSET_date_from_table(x['Date']['values'][1], x['ET']['values'][1])
    return (start_time, end_time)
    

def make_landmask_dep(lats, lons):
    
    def points_in_polys(points, polys):
        result = []
#        mask = np.empty_like(points)*False
        for poly in polys:
#            mask = path.contains_points(points, poly)
            polypath = path.Path(poly)
            mask = polypath.contains_points(points)
#            result.extend(points[mask])
            points = points[~mask]
        
        return np.array(result)

    m = Basemap(projection='moll',lon_0=0,resolution='c')
    m.drawcoastlines()
    m.fillcontinents(color='coral',lake_color='aqua')
    x, y = m(lons, lats)
#    loc = np.c_[x, y]
#    loc = np.array([(a, b) for a, b in zip(x, y)])
    loc = np.array([(a, b) for a, b in zip(x.ravel(), y.ravel())])
    polys = [p.boundary for p in m.landpolygons]
#    path = path.Path 
    land_loc = points_in_polys(loc, polys)
    mask = np.array([True if a in land_loc else False for a in loc]).reshape(x.shape)
    return mask
    
def make_landmask(lats, lons):
    
    m = Basemap(projection='cyl', resolution='c')
    x, y = m(lons.ravel(), lats.ravel())
    locations = np.c_[x, y]
    polygons = [path.Path(p.boundary) for p in m.landpolygons]
    result = np.zeros(len(locations), dtype=bool) 
    for polygon in polygons:
        result += np.array(polygon.contains_points(locations))
    return result.reshape(lats.shape)


def varcheck(fname, attr):
    with nc4.Dataset(fname) as dataset:
        if attr in list(dataset.variables.keys()):
#            print 'okay'
            return True
        else:
            print(fname)
            return False

def get_hysplit_files(run_date, run_hours):
    """Get HYSPLIT files required to run trajectories, return as list of files
    run_date: date of trajectory initialization
    run_hours: hours of trajectory. negative number means backward trajectory

    """
    today = dt.datetime.today()
    start_date = min(run_date, run_date + dt.timedelta(hours=run_hours))
    end_date = max(run_date, run_date + dt.timedelta(hours=run_hours))

    days_since_start = (today.date() - start_date.date()).days
    days_since_end = (today.date() - end_date.date()).days

    file_list = []

    while days_since_start > 0:  # add all analysis files from previous days
        date_to_add = today - dt.timedelta(days=days_since_start)
        if date_to_add > end_date:
            break
        try:
            f, d = get_hysplit_analysis(date_to_add)
            file_list.append(f)
        except ValueError:
            print(('could not find analysis for {}'.format(date_to_add)))
        days_since_start -= 1

    if days_since_end < 1:  # trajectory either ends today or in future
        f, d = get_hysplit_appended_files(today)
        file_list.append(f)
        f, d = get_hysplit_forecast_files(today)
        file_list.append(f)

    return file_list


def get_hysplit_analysis(date):
    """
    gets hysplit analysis file for day in date.
    if the file is already acquired, will not download it again.
    if the file does not exist yet raises error.
    """
    ftp = FTP('arlftp.arlhq.noaa.gov')
    ftp.login()
    ftp.cwd('/archives/gdas0p5')
    rx = re.compile('{:%Y%m%d}_gdas0p5\Z'.format(date))
    files = sorted(filter(rx.match, ftp.nlst()))
    if len(files) == 0:
        raise ValueError("ARL: No analysis available for {:%Y%m%d} yet...".format(date))
    newest = files[-1]
    savedir = os.path.join(HYSPLIT_source, 'analysis')
    if not os.path.exists(savedir):
        os.makedirs(savedir)
    print(("ARL: Attempting to find analysis file {} locally...".format(newest)))
    if os.path.isfile(os.path.join(savedir, newest)):
        print("ARL: File already acquired, not downloading it again.")
    else:
        print("ARL: File not found, will grab it from archives.")
        try:
            ftp.retrbinary("RETR " + newest,
                           open(os.path.join(savedir, newest), 'wb').write)
        except:
            print("ARL: Error in ftp transfer.")
            raise
        print('ARL: Analysis file successfully downloaded')

    savedfile = os.path.join(savedir, newest)
    print(('ARL: {}'.format(savedfile)))
    return savedfile, date


def get_hysplit_appended_files(date=None):
    """
    Gets most recent HYSPLIT appended files on date.
    Returns file location and initialization time (in the appended
    case that means the end of the file, so gfsa for 18Z on the 12th
    is relevant from 18Z on the 10th through the 12th, for instance)
    """
    f, d = get_hysplit_forecast_files(date, model='gfsa')
    return f, d


def get_hysplit_forecast_files(date=None, model='gfsf'):
    """
    Gets most recent HYSPLIT forecast files on date.
    Finds most recent file on ARL server. If it already exists on disk,
    does nothing and returns location on disk and initialization date.
    If it does not exist on disk, downloads and then returns the same.
    """
    def try_FTP_connect(ftpname):
        counter = 0
        while True:
            try:
                ftp = FTP(ftpname)
                return ftp
            except Exception as e:
                counter += 1
                sleep(1)
                if counter > 20:
                    raise e

    if date is None:
        date = dt.datetime.utcnow()

    ftp = try_FTP_connect('arlftp.arlhq.noaa.gov')
    ftp.login()
    ftp.cwd('/forecast/{:%Y%m%d/}'.format(date))
    rx = re.compile('hysplit.*.{}\Z'.format(model))
    files = list(filter(rx.match, ftp.nlst()))
    if len(files) == 0:  # too early in the day
        print(('ARL: no recent {} matches, looking at yesterday instead'.format(model)))
        date = date - dt.timedelta(days=1)
        ftp.cwd('/forecast/{:%Y%m%d/}'.format(date))
        files = list(filter(rx.match, ftp.nlst()))
    newest = files[-1]

    savedir = os.path.join(HYSPLIT_source, 'forecast',
                           '{:%Y%m%d}'.format(date))
    if not os.path.exists(savedir):
        os.makedirs(savedir)

    print(("ARL: Attempting to find {} for {:%Y-%m-%d}...".format(newest, date)))
    if os.path.isfile(os.path.join(savedir, newest)):
        print("ARL: File already acquired, not downloading it again.")
    else:
        print("ARL: File not found, will grab it from server.")
        try:
            ftp.retrbinary("RETR " + newest,
                           open(os.path.join(savedir, newest), 'wb').write)
        except:
            print("AR:L Error in ftp transfer.")
            raise
        print('ARL: File successfully downloaded')

    inittime = int(newest.split('.')[-2][1:3])
    initdate = date.replace(hour=inittime, minute=0, second=0,
                            microsecond=0)
    savedfile = os.path.join(savedir, newest)
    print(("ARL: file saves as {}".format(savedfile)))
    return(savedfile, initdate)


def write_control_file(start_time, coords, hyfile_list, hours, vertical_type, init_height,
                       tdumpdir):
    """
    This file generates the CONTROL files used for running the trajectories.
    start_time - the datetime object of when the trajectory should start
    coords - list of decimal [lat, lon] pairs. N and E are positive.
    hyfile_list - list of HYSPLIT source files on which to run model
    hours- negative hours means backwards run
    vertical_type:
        0 'data' ie vertical velocity fields
        1 isobaric
        2 isentropic
        3 constant density
        4 constant internal sigma coord
        5 from velocity divergence
        6 something wacky to convert from msl to HYSPLIT's above ground level
        7 spatially averaged vertical velocity
    """

    fl = os.path.join(HYSPLIT_workdir, 'CONTROL')
    f = open(fl, 'w')

    f.write(start_time.strftime('%y %m %d %H\n'))
    f.writelines([str(len(coords)), '\n'])
    for j in coords:
        f.write('{} {} {}\n'.format(str(j[0]), str(j[1]), init_height))
    f.writelines([str(hours), '\n'])

    f.writelines([str(vertical_type), '\n', '10000.0\n'])

    f.write('{}\n'.format(len(hyfile_list)))
    for hyfile in hyfile_list:
        f.writelines([
            os.path.dirname(hyfile), os.sep, '\n',
            os.path.basename(hyfile), '\n'])

    f.writelines([tdumpdir, os.sep, '\n', 'tdump',
                  start_time.strftime('%Y%m%dH%H%M'), '\n'])
    f.close()
    return os.path.join(tdumpdir, 'tdump'+start_time.strftime('%Y%m%dH%H%M'))


def read_tdump(tdump):
    """
    Read a tdump file as output by the HYSPLIT Trajectory Model
        Returns a pandas DataFrame object.
    """
    def parseFunc(y, m, d, H, M):
        ## Modified by Ehsan Erfani to adapt ERA5 trajectory. The original code is for HYSPLIT.
#         return dt.datetime(int('20'+y), int(m), int(d), int(H), int(M))
        return dt.datetime(int(y), int(m), int(d), int(H), int(M))
    columns = ['tnum', 'gnum', 'y', 'm', 'd', 'H', 'M', 'fhour', 'age', 'lat',
               'lon', 'height', 'pres']

    tmp = pd.read_table(tdump, nrows=100, header=None)
    l = [len(i[0]) for i in tmp.values]
    skiprows = l.index(max(l))
    D = pd.read_table(tdump, names=columns,
                      skiprows=skiprows,
                      engine='python',
                      sep='\s+', # delim_whitespace=True, 
                      parse_dates={'dtime': ['y', 'm', 'd', 'H', 'M']},
                      date_parser=parseFunc,
                      index_col='dtime')
    return D



def bmap(ax=None, drawlines=True, llr=None, par_labs=[1, 1, 0, 0], mer_labs=[0, 0, 1, 1], 
         merspace=15, parspace=15, **kwargs):

    if ax is None:
        fig, ax = plt.subplots()

    if llr is None:
        lat_range = latlon_range['lat']
        lon_range = latlon_range['lon']
    else:
        lat_range = llr['lat']
        lon_range = llr['lon']
    if 'projection' not in kwargs.keys():
        kwargs['projection'] = 'cyl'
        kwargs['rsphere'] =(6378137.00, 6356752.3142)
    m = Basemap(llcrnrlon=lon_range[0], llcrnrlat=lat_range[0],
                urcrnrlon=lon_range[1],  urcrnrlat=lat_range[1],
                ax=ax, resolution='l', **kwargs)
    if drawlines:
        m.drawparallels(np.arange(-90., 90., parspace), labels=par_labs, fontsize=14)
        m.drawmeridians(np.arange(-180., 180., merspace), labels=mer_labs, fontsize=14)
    m.drawcoastlines()
    m.fillcontinents(color="white", lake_color="white")

    return m

    
def read_flightpath(flightfile):
    """read in flight file netcdf and return as dict.
    """
    with nc4.Dataset(flightfile, 'r') as flt_nc:
        lats = flt_nc.variables['LATC'][:].copy()
        lons = flt_nc.variables['LONC'][:].copy()
        alt = flt_nc.variables['ALT'][:].copy()
        timevar = flt_nc.variables['Time']
        date = nc4.num2date(timevar[:], units=timevar.units)
    if isinstance(lats, np.ma.core.MaskedArray):
        m = np.logical_or(lats.mask, lons.mask)
        lats = lats.data[~m]
        lons = lons.data[~m]
        alt = alt.data[~m]
        date = date[~m]
    fp = {'lats': lats, 'lons': lons, 'date': date,
          'alt': alt}
    return fp
    
    

def gridder(SW, NW, NE, SE, numlats=6, numlons=6):
    """each point is a [lat lon] corner of the desired area"""
    lat_starts = np.linspace(SW[0], NW[0], numlats)
    lon_starts = np.linspace(SW[1], SE[1], numlons)
    lat_ends = np.linspace(SE[0], NE[0], numlats)
    lon_ends = np.linspace(NW[1], NE[1], numlons)
    lat_weight = np.linspace(0., 1., numlats)
    lon_weight = np.linspace(0., 1., numlons)
    lat = (1. - lon_weight[:, None])*lat_starts[None, :] +\
        lon_weight[:, None]*lat_ends[None, :]
    lon = (1. - lat_weight[:, None])*lon_starts[None, :] +\
        lat_weight[:, None]*lon_ends[None, :]
    l = []
    for i in range(numlats):
        for j in range(numlons):
            l.append((lat[j, i], lon[i, j]))
    return(l)


def plot_gridpoints(coords, outfile=None):

    fig = plt.figure()
    ax = fig.add_axes([0, 0, 1, 1])
    m = bmap(ax=ax, proj='cyl', drawlines=True)

    m.drawgreatcircle(-121.3, 38.6, -156, 19.8, linestyle='--', c='black')
    colors = cm.rainbow(np.linspace(0, 1, len(coords)))
    for i, crd in enumerate(coords):
        m.plot(crd[1], crd[0], '*', c=colors[i], latlon=True, ms=12, label=i)
        x, y = m(crd[1]+.5, crd[0]+.5)
        ax.annotate(str(i), xy=(x, y), xytext=(x, y), xycoords='data',
                    textcoords='data', fontsize=6)

    if outfile is not None:
        ax.patch.set_visible(False)
        fig.savefig(outfile, dpi=300, transparent=True, bbox_inches='tight',
                    pad_inches=0)


def plot_trajectory(date=None, filename=None):
    if date is None and filename is None:
        print('give me a date (YYYY-MM-DD) or a file, dummy')
        return
    elif date:
        datet = dt.datetime.strptime(date, '%Y-%m-%d')
        filename = os.path.join(trajectory_dir, 'tdump'+datet.strftime('%Y%m%dH%H%M'))

    fig, ax, m_ax = make_map_plot()
    add_tdump_to_plot(m_ax, filename)
    return


def make_map_plot(ax=None, llr=None, **kwargs):
    if ax is None:
        fig, ax = plt.subplots(figsize=(7, 8))
    else:
        fig = ax.get_figure()
    m_ax = bmap(ax=ax, llr=llr, **kwargs)
#    m_ax.drawgreatcircle(-121.3, 38.6, -156, 19.8, linestyle='--', c='black')
    m_ax.plot(-121.3, 38.6, 's', ms=8, c='black', latlon=True)
    m_ax.plot(-156, 19.8, '*', ms=12, c='black', latlon=True)
#    m_ax.plot(-118.2, 33.77, 's', ms=8, c='red', latlon=True)
    return fig, ax, m_ax


def nan_correlate(x, y):
    x, y = np.array(x), np.array(y)
    index = np.logical_and(~np.isnan(x), ~np.isnan(y))
    return np.corrcoef(x[index], y[index])[0][1]


def plot_single(t, m=None, c=None, i=None):
    m.plot(t.lon.values, t.lat.values, c=c, latlon=True, label=t.tnum[0])
    m.plot(t.lon.values[::6], t.lat.values[::6], '.', c=c, latlon=True)
    m.plot(t.lon.values[0], t.lat.values[0], '*', c=c, latlon=True, ms=12)
    m.plot(t.lon.values[-1], t.lat.values[-1], 's', c=c, latlon=True, ms=8)
    if i is not None:
        plt.annotate(str(i), xy=(t.lon.values[0]+.5, t.lat.values[0]+.5))

    return m


def add_tdump_to_plot(m_ax, tdump):

    T = read_tdump(tdump)
    t = T.groupby('tnum')
    colors = cm.rainbow(np.linspace(0, 1, len(list(t.groups.keys()))))
    for i, k in enumerate(t.groups.keys()):
        m_ax = plot_single(t.get_group(k), m=m_ax, c=colors[i], i=i)

    return


def get_pesky_GOES_files():
    badfiles = []
    with open(r'/home/disk/p/jkcm/Code/Lagrangian_CSET/GOES_Extractor.log', 'r') as f:
            for line in f:
                if r'/home/disk/eos4/mcgibbon/nobackup/GOES' in line:
                    if line not in badfiles:
                        badfiles.append(line)

    with open(r'/home/disk/p/jkcm/Code/Lagrangian_CSET/flawed_GOES.log', 'w') as g:
        for line in sorted(badfiles):
            if os.path.exists(line[:-1]):
                size = '{:3.0f}'.format(os.path.getsize(line[:-1])/1024)
#                print size
            else:
                size = 'NA '
            replace_GOES_file(line[:-1])
            g.writelines(size + '    ' + line)


def replace_GOES_file(filename, savedir=None):
    oldfilename = os.path.basename(filename)
    year = int(oldfilename[12:16])
    date = dt.datetime(year, 1, 1) + dt.timedelta(days=int(oldfilename[16:19]) - 1)
    newfilename = 'prod.goes-west.visst-pixel-netcdf.{:%Y%m%d}.{}'.format(
        date, oldfilename)
    floc = 'prod/goes-west/visst-pixel-netcdf/{:%Y/%m/%d}/'.format(date)
    server = r'http://cloudsgate2.larc.nasa.gov/'
    url = server + floc + newfilename

    try:
        response = urlopen(url)
    except HTTPError:
        print('could not find file!')
        return
    print('file found, downloading')
    if savedir is None:
        savedir = GOES_source
    print(('old size is {}KB'.format(os.path.getsize(filename)/1024.)))
    if os.path.dirname(filename) == savedir:
        print('moving old file')
        if not os.path.exists(os.path.join(savedir, 'old')):
            os.makedirs(os.path.join(savedir, 'old'))
        os.rename(filename, os.path.join(savedir, 'old', oldfilename))

    save_file = os.path.join(savedir, oldfilename)
    with open(save_file, 'wb') as fp:
        while True:
            chunk = response.read(16384)
            if not chunk:
                break
            fp.write(chunk)

    print(('new size = {}KB'.format(os.path.getsize(save_file)/1024.)))

    
def as_datetime(date, timezone=pytz.UTC):
    "Converts all datetimes types to datetime.datetime with TZ = UTC"
    def to_dt(d, timezone):
        """does all the heavy lifting
        """
        supported_types = (np.datetime64, dt.datetime)
        if not isinstance(d, supported_types):
            raise TypeError('type not supported: {}'.format(type(d)))
        if isinstance(d, np.datetime64):
            # TODO: add timezoneawareness here
            ts = (d - np.datetime64('1970-01-01T00:00:00Z')) / np.timedelta64(1, 's')
            d = dt.datetime.utcfromtimestamp(ts)
        if isinstance(d, pd.Timestamp):
            d = d.to_datetime()
        if isinstance(d, dt.datetime):
            if d.tzinfo is None:
                return(d.replace(tzinfo=timezone))
            else:
                return(d.astimezone(timezone))

    if isinstance(date, (collections.Sequence, np.ndarray)):
        return np.array([to_dt(x, timezone) for x in date])
    return to_dt(date, timezone)
    
    
datemap = {'20150701': 'RF01',
             '20150707': 'RF02',
             '20150709': 'RF03',
             '20150712': 'RF04',
             '20150714': 'RF05',
             '20150717': 'RF06',
             '20150719': 'RF07',
             '20150722': 'RF08',
             '20150724': 'RF09',
             '20150727': 'RF10',
             '20150729': 'RF11',
             '20150801': 'RF12',
             '20150803': 'RF13',
             '20150807': 'RF14',
             '20150809': 'RF15',
             '20150812': 'RF16'}


def get_data_from_dropsonde(file):
#    file = os.path.join(dropsonde_dir, 'D20150712_201424_PQC.nc')

    data = xr.open_dataset(file)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        index = data.GPSAlt.values < 4000

    ret = {}
    ret['TIME']=as_datetime(data.time_offset.values[index])
    ret['GGLAT']=data.Lat.values[index]
    ret['GGLON']=data.Lon.values[index]
    ret['GGALT']=data.GPSAlt.values[index]
    ret['RHUM']=data.RH.values[index]
    ret['ATX']=data.Temp.values[index]+273.15
    ret['PSX']=data.Press.values[index]
    ret['DPXC']= data.Dewpt.values[index]+273.15
    ret['QV'] = mu.qv_from_p_T_RH(ret['PSX']*100, ret['ATX'], ret['RHUM'])*1000
    ret['MR'] = ret['QV']/(1-ret['QV']/1000)
    ret['TVIR'] = mu.tvir_from_T_w(ret['ATX'], ret['MR']/1000)
    ret['DENS'] = mu.density_from_p_Tv(ret['PSX']*100, ret['TVIR'])  

    ret['THETA']= mu.theta_from_p_T(ret['PSX'], ret['ATX'])
    ret['THETAE']= mu.thetae_from_t_tdew_mr_p(ret['ATX'], ret['DPXC'], ret['MR']/1000, ret['PSX']*100) #equiv pot temp, K we can get this if we really want

    ret['QL'] = np.full_like(ret['PSX'], fill_value=np.nan)
    ret['THETAL'] = np.full_like(ret['PSX'], fill_value=np.nan)
    ret['PLWCC']=  np.full_like(ret['PSX'], fill_value=np.nan)
    return ret


def date_interp(dates_new, dates_old, vals_old, bounds_error=False):
    if not isinstance(dates_new, (collections.Sequence, np.ndarray)):
        dates_new = np.array([dates_new])
    dates_new = as_datetime(dates_new)
    dates_old = as_datetime(dates_old)
    ref = min(min(dates_old), min(dates_new))
    d_new = [(i-ref).total_seconds() for i in dates_new]
    d_old = [(i-ref).total_seconds() for i in dates_old]
    vals_new = interp1d(d_old, vals_old, bounds_error=bounds_error)(d_new).squeeze()
    if vals_new.shape == ():
        return vals_new.item()
    return vals_new

def get_cloud_only_vals(dataset, flip_cloud_mask=False):

    # cloud if ql_cdp > 0.01 g/kg and RH > 95%
    lwc_cdp = dataset['PLWCD_LWOI']
    rhodt = dataset['RHODT']
    mr = dataset['MR']
    cheat_airdens = rhodt/mr
    lwmr_cdp = lwc_cdp/cheat_airdens
    lw_index = lwmr_cdp > 0.01
    RH_index = dataset['RHUM'] > 95
    cloud_index = np.logical_and(RH_index, lw_index)
    if flip_cloud_mask:
        cloud_index = np.logical_not(cloud_index)
    return dataset.isel(time=cloud_index)
