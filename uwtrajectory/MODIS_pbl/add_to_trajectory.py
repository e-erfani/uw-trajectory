###############################################
## Modifications by Ehsan Erfani (2022):
## - Now, the path and file name are specified
##   in config file
##############################################

import utils, config
import numpy as np
import xarray as xr
import glob

def add_MODIS_pbl_to_trajectory(ds, box_degrees=3):
    lats, lons,times = ds.lat.values, ((ds.lon.values+180)%360)-180, ds.time.values
    MODIS_day_idx = np.argwhere([i.hour == 23 for i in utils.as_datetime(times)]).squeeze()
    MODIS_night_idx = np.argwhere([i.hour == 11 for i in utils.as_datetime(times)]).squeeze()
    dayfile = config.MODIS_pbl_dayfile
    nightfile = config.MODIS_pbl_nightfile
    vals = []   
    stds = []
    nanfrac = []
    medians = []
    counts = []
    mins = []
    maxs = []

    for i in range(len(times)):
        if i in MODIS_day_idx:
            f = dayfile
        elif i in MODIS_night_idx:
            f = nightfile
        else:
            vals.append(np.nan)
            stds.append(np.nan)
            nanfrac.append(np.nan)
            medians.append(np.nan)
            counts.append(np.nan)
            mins.append(np.nan)
            maxs.append(np.nan)
            continue
        with xr.open_dataset(f) as data:
            lat, lon, time = lats[i], lons[i], utils.as_datetime(times[i])
            t_idx = np.argwhere(np.logical_and(data['days'].values == time.timetuple().tm_yday, 
                                   data['years'].values == time.year))[0][0]
            x = data['cth'].sel(longitude=slice(lon - box_degrees/2, lon + box_degrees/2),
                  latitude=slice(lat + box_degrees/2, lat - box_degrees/2))
            z = x.isel(time=t_idx).values
            vals.append(np.nanmean(z))
            medians.append(np.nanmedian(z))
            counts.append(np.sum(~np.isnan(z)))
            stds.append(np.nanstd(z))
            nanfrac.append(np.sum(np.isnan(z))/z.size)
            try:
                mins.append(np.nanmin(z))
                maxs.append(np.nanmax(z))
            except ValueError:
                mins.append(np.nan)
                maxs.append(np.nan)
    ds['MODIS_CTH'] = (('time'), np.array(vals), {'long_name': 'MODIS cloud top height, box mean', 'units': 'km'})
    ds['MODIS_CTH_median'] = (('time'), np.array(medians), {'long_name': 'MODIS cloud top height, box median', 'units': 'km'})
    ds['MODIS_CTH_std'] = (('time'), np.array(stds), {'long_name': 'MODIS cloud top height, box standard deviation', 'units': 'km'})
    ds['MODIS_CTH_min'] = (('time'), np.array(mins), {'long_name': 'MODIS cloud top height, box min', 'units': 'km'})
    ds['MODIS_CTH_max'] = (('time'), np.array(maxs), {'long_name': 'MODIS cloud top height, box max', 'units': 'km'})
    ds['MODIS_CTH_nanfrac'] = (('time'), np.array(nanfrac), {'long_name': 'MODIS missing pixel fraction', 'units': '0-1'})
    ds['MODIS_CTH_n_samples'] = (('time'), np.array(counts), {'long_name': 'MODIS CTH number of samples'})
    
    
    
    ds.attrs['MODIS_params'] = f'MODIS Aqua cloud-top height from Eastman et al. (2017) computed over a {box_degrees}-deg average centered on trajectory'
    ds.attrs['MODIS_reference'] = f"Eastman, R., Wood, R., & O, K. T. (2017). The Subtropical Stratocumulus-Topped Planetary Boundary Layer: A Climatology and the Lagrangian Evolution. Journal of the Atmospheric Sciences, 74(8), 2633–2656. https://doi.org/10.1175/JAS-D-16-0336.1"
    return ds