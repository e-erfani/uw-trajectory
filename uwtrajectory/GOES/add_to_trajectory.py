from .. import config
import xarray as  xr
import numpy as np
import pandas as pd
import glob

def add_GOES_to_trajectory(ds, box_degrees=2):
    lats, lons, times = ds.lat.values, ((ds.lon.values+180)%360)-180, ds.time.values


    cths, lwps, cfs = [], [], []
    for i, (lat, lon, time) in enumerate(zip(lats, lons, pd.to_datetime(times))):
        try:
            goes_file = glob.glob(config.GOES_file_fmt.format(time.year, time.dayofyear, time.hour))[0]
        except IndexError:
            print("no goes found")
            cths.append(np.nan)
            lwps.append(np.nan)
            cfs.append(np.nan)
            continue
        goes_data = xr.open_dataset(goes_file, chunks={})
        mask = np.logical_and.reduce(((goes_data.latitude>lat-box_degrees/2).values,
                                      (goes_data.latitude<lat+box_degrees/2).values,
                                      (goes_data.longitude<lon+box_degrees/2).values,
                                      (goes_data.longitude>lon-box_degrees/2).values))
#         print(np.sum(mask))
        goes_cth = np.nanmean(goes_data.cloud_top_height.values[mask])
        goes_lwp = np.nanmean(goes_data.cloud_lwp_iwp.values[mask])
        goes_phase = goes_data.cloud_phase.values[mask]
        goes_cf = np.nansum(goes_phase==1)/np.nansum(np.logical_or(goes_phase==1, goes_phase==4))

        cths.append(goes_cth)
        lwps.append(goes_lwp)
        cfs.append(goes_cf)
    
    ds['GOES_CTH'] = (('time'), np.array(cths), {'long_name': 'GOES_cloud_top_height'})
    ds['GOES_LWP'] = (('time'), np.array(lwps), {'long_name': 'GOES_liquid_water_path'})
    ds['GOES_CF'] = (('time'), np.array(cfs), {'long_name': 'GOES_cloud_fraction'})
    return ds