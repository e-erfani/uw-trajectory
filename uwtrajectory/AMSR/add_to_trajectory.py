###############################################
## Modifications by Ehsan Erfani (2022):
## - Now, the path and file name are specified
##   in config file
##############################################

import xarray as xr
import numpy as np
import pandas as pd
import utils, met_utils, config
from config import AMSR_fmt

def add_AMSR_to_trajectory(ds, box_degrees=2):
    lats, lons, times = ds.lat.values, ds.lon.values, pd.DatetimeIndex(ds.time.values)
#    data = xr.open_mfdataset('/home/disk/eos9/jkcm/Data/amsr/rss/all/amsr_unified_*-*.nc', combine='by_coords')
    data = xr.open_mfdataset(AMSR_fmt, combine='by_coords')

    dummy = np.full_like(lats, fill_value=np.nan)
    vals_dict = {'AMSR_CWV': {'name': 'vapor', 'means': dummy.copy(), 'stds': dummy.copy(), 'n_samples': dummy.copy()},
                'AMSR_LWP': {'name': 'cloud', 'means': dummy.copy(), 'stds': dummy.copy(), 'n_samples': dummy.copy()},
                'AMSR_SST': {'name': 'sst', 'means': dummy.copy(), 'stds': dummy.copy(), 'n_samples': dummy.copy()},
                }
    for i, (lat, lon, time) in enumerate(zip(lats, lons%360, times)):
        lst = time.hour + 24*((lon+180)%360-180)/360
        for orbit in [0,1]:
            data_subs = data.sel(time=time.replace(hour=0), orbit_segment=orbit,
                                 longitude=slice(lon - box_degrees/2, lon + box_degrees/2),
                                 latitude=slice(lat - box_degrees/2, lat + box_degrees/2))
            x = np.nanmean(data_subs.UTCtime)
            if np.isnan(x):
#                 print('out of swath, skipping')
                continue
            slice_time = data_subs.time.values + np.timedelta64(int(x), 'h') + np.timedelta64(int((x-int(x))*60), 'm')
            time_miss = (slice_time-time)/np.timedelta64(1, 'h')
            if np.abs(time_miss)<1: #within an hour of data sample
                for var in vals_dict.values():
                    var['means'][i] = np.nanmean(data_subs[var['name']].values)
                    var['stds'][i] = np.nanstd(data_subs[var['name']].values)
                    var['n_samples'][i] = np.sum(~np.isnan(data_subs[var['name']].values))
                                                                  
    for var,vals in vals_dict.items():
        attrs = data[vals['name']].attrs
        ds[var] = (('time'), vals['means'], attrs)
        attrs.update(long_name=attrs['long_name']+', standard deviation over box')
        ds[var+'_std'] = (('time'), vals['stds'], attrs)
        ds[var+'_nsamples'] = (('time'), vals['n_samples'], {'long_name': ds[var].attrs['long_name']+', number of samples'})
        

        ds.attrs['AMSR2_params'] = f'AMSR2 statistics computed over a {box_degrees}-deg average centered on trajectory'
    ds.attrs['AMSR2_reference'] = f"AMSR2 data are produced by Remote Sensing Systems. Data are available at www.remss.com/missions/amsr."
    
    return ds
