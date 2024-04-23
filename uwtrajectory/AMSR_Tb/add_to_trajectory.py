import numpy as np
import xarray as xr
import pandas as pd

def add_AMSR_Tb_to_trajectory(ds, box_degrees=2):
    lats, lons, times = ds.lat.values, ds.lon.values%360, ds.time.values
    # daynight = {}
    amsr_prcp_mean = np.full(len(lats), fill_value=np.nan)
    amsr_prcp_std = np.full(len(lats), fill_value=np.nan)
    amsr_prcp_count = np.full(len(lats), fill_value=np.nan)
    amsr_prcp_dist = np.full((len(lats), 9), fill_value=np.nan)
    for dn in ['day', 'night']:
        # print(f'processing amsr {dn}time')
        precip_files = np.unique([f"/home/disk/eos5/rmeast/rain_rates_89/{i.year}_v82/AMSR2_89GHz_pcp_est_{i.year}_{i.dayofyear:03}_{dn}.nc" 
                                                    for i in pd.DatetimeIndex(ds.time.values)])
        precip_data = xr.open_mfdataset(precip_files, combine='nested', concat_dim='time')
        year = precip_data.time_vars.isel(yr_day_utc=0).values
        day = precip_data.time_vars.isel(yr_day_utc=1).values
        utc = precip_data.time_vars.isel(yr_day_utc=2).values
        total_secs = (utc*3600)
        secs = total_secs//1
        msecs = 1000*total_secs%1
        dtime = np.datetime64(f'{np.median(year):0.0f}-01-01')+np.timedelta64(1, 'D')*(day-1)+\
                np.timedelta64(1, 's')*(secs)+np.timedelta64(1, 'ms')*(msecs)
        precip_data['time'] = (('time'), dtime)
        precip_data = precip_data.drop(labels=['time_vars'])
        precip_data['longitude'] = precip_data['longitude']%360


        # filtering out the misaligned times, couple scans at day start/end
        tdiff = np.diff(precip_data.time)/np.timedelta64(1, 'h')
        ups = np.argwhere(tdiff>20).flatten()
        downs = np.argwhere(tdiff<-20).flatten()
        if len(ups)<len(downs):
            if len(ups) == 0: # starting with a misalignment and only one error:
                ups = np.insert(ups, 0, 0) 
            elif ups[0]>downs[0] : # also starting with a misalignment; this is done in a kludgy way to avoid an index error
                ups = np.insert(ups, 0, 0) # chop off from 0 to first down
        good = np.ones(precip_data.time.shape).astype(bool)
        for (u,d) in zip(ups, downs):
            good[u:d+1] = np.zeros_like(good[u:d+1])
        precip_data = precip_data.isel(time=good)

        #actual statistics come here

        for i, (lat, lon, time) in enumerate(zip(lats, lons, times)):
            try:
                ds_sub = precip_data.sel(time=slice(time-np.timedelta64(1, 'h'), time+np.timedelta64(1, 'h')))
            except (ValueError, KeyError) as e:
                raise e
            good = np.logical_and(np.abs(ds_sub.latitude-lat)<1, np.abs(ds_sub.longitude-lon)<1)
            n_good = np.sum(good).values
            n_agood = 0
            if n_good > 100: # at least 100 actual pixels in the swath
                actual_good = ds_sub.rain_stats.isel(prob_rate_rwr_max=0).where(good)>0
                n_agood = np.sum(actual_good).values
                # print(f'{dn}: {i}, {time} inside swath')
            else:
                pass #outside of swath
            if n_agood > 100: # at least 100 pixels with a valid retrieval
                g_rate = ds_sub.rain_stats.isel(prob_rate_rwr_max=1).where(good).values.flatten()[good.values.flatten()]
                g_rate = g_rate[~np.isnan(g_rate)] # this should be redundant

                for var in [g_rate]:
                    var = var[~np.isnan(var)]
                    sorted_var = sorted(var)
                    cumsum = np.cumsum(sorted_var)
                    pctiles = np.linspace(0.1*np.sum(var), 0.9*np.sum(var), 9)
                    pct_vals = np.array([sorted_var[j] for j in [np.argmax(cumsum>=i) for i in pctiles]])

                    
                amsr_prcp_mean[i] = np.nanmean(g_rate)
                amsr_prcp_std[i] = np.nanstd(g_rate)
                amsr_prcp_count[i] = n_agood
                amsr_prcp_dist[i] = pct_vals
            else:
                pass # leave it as nans
            
    ds['AMSR_prcp_mean'] = (('time'), amsr_prcp_mean, {'long_name': 'AMSR Tb mean rain rate', 'units': 'mm hr^-1'})
    ds['AMSR_prcp_std'] = (('time'), amsr_prcp_std, {'long_name': 'AMSR Tb rain rate standard deviation', 'units': 'mm hr^-1'})
    ds['AMSR_prcp_n_samples'] = (('time'), amsr_prcp_count, {'long_name': 'AMSR Tb rain rate sample count', 'units': ''})
    ds['AMSR_prcp_dist'] = (('time', 'AMSR_prcp_pctile'), amsr_prcp_dist, 
                            {'long_name': 'AMSR Tb rain rate cumulative distribution', 'units': 'mm hr^-1'})
    ds['AMSR_prcp_pctile'] = (('AMSR_prcp_pctile'), np.linspace(10, 90, 9), 
                              {'long_name': 'rain below value in amsr_prcp_dist explains this much of all precipitation'})
    

    ds.attrs['AMSR_params'] = f'AMSR data statistics computed over a {box_degrees}-deg average centered on trajectory, only when > 100 samples are present'
    ds.attrs['AMSR_reference'] = f'AMSR precipitation is derived from 89 GHz brightness temperature, from Eastman, R., Lebsock, M., & Wood, R. (2019). Warm Rain Rates from AMSR-E 89-GHz Brightness Temperatures Trained Using CloudSat Rain-Rate Observations. Journal of Atmospheric and Oceanic Technology, 36(6), 1033–1051. https://doi.org/10.1175/JTECH-D-18-0185.1'

    return ds    
    
    
    
    
def dec_old_add_amsr_to_trajectory(ds, box_degrees=2):     
    lats, lons, times = ds.lat.values, ds.lon.values, ds.time.values

    try:
        precip_data = xr.open_mfdataset(np.unique([f'/home/disk/eos9/jkcm/Data/rain/{i.year}/AMSR2_89GHz_pcp_est_{i.year}_{i.dayofyear:03}_day_gridded.nc' 
                                                   for i in pd.DatetimeIndex(ds.time.values)]) , combine='by_coords')
        
        means_dict, stds_dict = dict(), dict()
        for i, (lat, lon, time) in enumerate(zip(lats, lons%360, times)):
            ds_sub = precip_data.sel(date=time, method='nearest', tolerance=np.timedelta64(24, 'h'))  
            ds_sub2 = ds_sub.sel(longitude=slice(lon - box_degrees/2, lon + box_degrees/2),
                                 latitude=slice(lat - box_degrees/2, lat + box_degrees/2))
            # print(ds_sub2.time)       

            print((ds_sub2.time.values - time)/np.timedelta64(1, 'h'))

        

            for var in  ['rain_rate_mean', 'rain_rwr_mean', 'rain_prob_mean']:
                    
                    means_dict.setdefault(var, []).append(np.nanmean(ds_sub2[var]))
                    stds_dict.setdefault(var, []).append(np.nanstd(ds_sub2[var]))

            for var in means_dict.keys():
                attrs = attr_dict.setdefault(var, ds_sub2[var].attrs)
                ds['AMSR_TB_'+var] = (('time'), np.array(means_dict[var]), attrs)
                attrs.update(long_name=attrs['long_name']+', standard deviation over box')
                ds['AMSR_TB_'+var+'_std'] = (('time'), np.array(stds_dict[var]), attrs)

    except (FileNotFoundError, AttributeError) as e:
        raise e
    #     if isinstance(e, FileNotFoundError):
    #         print('could not find precip file for date, continuing:' + str(date))
    #         print(e)
    #     elif isinstance(e, AttributeError):
    #         print('attribute error, likely no data found' + str(date))
    #         print(e)
            
    #     df['amsr_tb_rate'] = np.nan
    #     df['amsr_tb_rwr'] = np.nan
    #     df['amsr_tb_prob'] = np.nan
        
    # except KeyError as e:
    #     print(precip_data.date)
    #     print(date)
    #     print(e)
    
    return ds