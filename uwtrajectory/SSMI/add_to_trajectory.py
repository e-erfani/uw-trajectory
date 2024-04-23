###############################################
## Modifications by Ehsan Erfani (2022):
## - Now, the path and file name are specified
##   in config file
##############################################

import config
import numpy as np
import xarray as xr
import glob

def add_SSMI_to_trajectory(ds, box_degrees=2, hour_tolerance=0.5):
    lats, lons, times = ds.lat.values, ds.lon.values, ds.time.values
    cloud_vals = np.full_like(lats, fill_value=np.nan)
    vapor_vals = np.full_like(lats, fill_value=np.nan)
    wspd_vals = np.full_like(lats, fill_value=np.nan)
    cloud_vals_std = np.full_like(lats, fill_value=np.nan)
    vapor_vals_std = np.full_like(lats, fill_value=np.nan)
    wspd_vals_std = np.full_like(lats, fill_value=np.nan)
    count_vals = np.full_like(lats, fill_value=np.nan)
    total_vals = np.full_like(lats, fill_value=np.nan)
#    sats = ['f15', 'f16' ,'f17', 'f18']
    for sat in config.sats:
        ssmi_data = xr.open_mfdataset(glob.glob(config.SSMI_file_fmt.format(sat)), concat_dim='time', combine='by_coords')
        for i, (lat, lon, time) in enumerate(zip(lats, lons%360, times)):
            for orbit_segment in [0,1]:
                ds_sub = ssmi_data.sel(time=time-np.timedelta64(1-orbit_segment, 'D'), method='nearest', tolerance=np.timedelta64(24, 'h')).sel(orbit_segment=orbit_segment)
                #note to future users: because of the longitude of CSET, the 1-day offset is included. This is a hardcoded value to deal with the fact that the UTC
                # day of SSMI doesn't line up when you're near the antimeridian. Future use should deal with this better, by deciding whether or not to add a day
                # based on the longitude being considered. 
                ds_sub2 = ds_sub.sel(longitude=slice(lon - box_degrees/2, lon + box_degrees/2),
                                          latitude=slice(lat - box_degrees/2, lat + box_degrees/2))
                sel_date = ds_sub2.UTCtime
                nsample = 1-(np.sum(ds_sub2.nodata)/np.size(ds_sub2.nodata)).values
                if nsample < 0.5:
#                     print('no samples')
        #             print(f'skipping {time}')
                    continue
                else:

                    meantime = np.nanmean(ds_sub2.UTCtime.values)
                    sampletime = ds_sub2.time.values + np.timedelta64(int(meantime), 'h') + np.timedelta64(int(60*(meantime - int(meantime))), 'm')
                    miss = (time-sampletime)/np.timedelta64(1, 'h')
                    if np.abs(miss)<hour_tolerance:
                         #print(f'{sat}: found data at {time}, off by {miss} hours. sample fill is {nsample:.0%}')
                        # print(np.sum(~np.isnan(ds_sub2.cloud).values), np.sum(~np.isnan(ds_sub2.vapor).values), np.sum(~np.isnan(ds_sub2.wspd_mf)).values)
                        # print(np.size(ds_sub2.nodata.values)-np.sum(ds_sub2.nodata.values), np.sum(~ds_sub2.nodata).values)
                        cloud_vals[i] = np.nanmean(ds_sub2.cloud)
                        vapor_vals[i] = np.nanmean(ds_sub2.vapor)
                        wspd_vals[i] = np.nanmean(ds_sub2.wspd_mf)
                        count_vals[i] = np.sum(~np.isnan(ds_sub2.cloud))
                        total_vals[i] = np.size(ds_sub2.cloud.values)
                        cloud_vals_std[i] = np.nanstd(ds_sub2.cloud)
                        vapor_vals_std[i] = np.nanstd(ds_sub2.vapor)
                        wspd_vals_std[i] = np.nanstd(ds_sub2.wspd_mf)
#                     else:
#                         print('outside hour tolerance')
#                     break
    ds['SSMI_LWP'] = (('time'), np.array(cloud_vals), ds_sub2.cloud.attrs)
    ds['SSMI_LWP_std'] = (('time'), np.array(cloud_vals_std), ds_sub2.cloud.attrs)
    ds['SSMI_WVP'] = (('time'), np.array(vapor_vals), ds_sub2.vapor.attrs)
    ds['SSMI_VWP_std'] = (('time'), np.array(vapor_vals_std), ds_sub2.vapor.attrs)
    ds['SSMI_WSPD'] = (('time'), np.array(wspd_vals), ds_sub2.wspd_mf.attrs)
    ds['SSMI_WSPD_std'] = (('time'), np.array(wspd_vals_std), ds_sub2.wspd_mf.attrs)
    ds['SSMI_n_samples'] = (('time'), np.array(count_vals), {'long_name': 'SSMI number of data samples'})
    ds['SSMI_n_total'] = (('time'), np.array(total_vals), {'long_name': 'SSMI total number of pixels'})
#     print(np.size(ds_sub2.values))
    ds['SSMI_WSPD_std'] = (('time'), np.array(wspd_vals_std), ds_sub2.wspd_mf.attrs)
    for i in ['SSMI_LWP_std', 'SSMI_VWP_std', 'SSMI_WSPD_std']:
        ds[i].attrs['long_name'] = ds[i].attrs['long_name']+' standard deviation over box'
        
    ds.attrs['SSMI_params'] = f'SSM/I data added from satellites {", ".join(config.sats).upper()}; statistics computed over a {box_degrees}-deg average centered on trajectory'
    ds.attrs['SSMI_reference'] = f"SSM/I and SSMIS data are produced by Remote Sensing Systems. Data are available at www.remss.com/missions/ssmi."

    return ds
