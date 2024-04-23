###############################################################
## Modifications by Ehsan Erfani (2022):
## - Addition of a few variables: 'iwp_low_1h', 'iwp_total_1h',
##                          'lwp_total_1h', 'cldarea_total_1h',
## - Calculation of CRE
###############################################################

import xarray as  xr
import numpy as np
import utils, met_utils, config
from config import CERES_fmt

def add_CERES_to_trajectory(ds, box_degrees=2):
    lats, lons, times = ds.lat.values, ds.lon.values, ds.time.values

    ceres_file = xr.open_mfdataset(CERES_fmt, combine='by_coords')

    attr_dict ={'net_cre': {'long_name': 'Net Cloud Radiative Effect', 'units': 'W m-2'},
            'sw_cre': {'long_name': 'Shortwave Cloud Radiative Effect', 'units': 'W m-2'},
            'lw_cre': {'long_name': 'Longwave Cloud Radiative Effect', 'units': 'W m-2'}}
    
    means_dict, stds_dict = dict(), dict()
    
    lw_cre_time     = []; sw_cre_time     = []; net_cre_time     = []   ## This line added by Ehsan Erfani
    lw_cre_std_time = []; sw_cre_std_time = []; net_cre_std_time = []   ## This line added by Ehsan Erfani

    for i, (lat, lon, time) in enumerate(zip(lats, lons%360, times)):
        ds_sub  = ceres_file.sel(time=time, method='nearest', tolerance=np.timedelta64(24, 'h'))
        ds_sub2 = ds_sub.sel(lon=slice(lon - box_degrees/2, lon + box_degrees/2),
                             lat=slice(lat - box_degrees/2, lat + box_degrees/2))
        
        #####################################
        ####### Modified by Ehsan Erfani:       
        LW_CRE  = - ds_sub2['toa_lw_all_1h'] + ds_sub2['toa_lw_clr_1h']
        SW_CRE  =   ds_sub2['toa_sw_all_1h'] - ds_sub2['toa_sw_clr_1h']
        net_CRE =   SW_CRE + LW_CRE
        
        lw_cre_time.append (np.nanmean(LW_CRE)); sw_cre_time.append (np.nanmean(SW_CRE)); net_cre_time.append(np.nanmean(net_CRE))
        lw_cre_std_time.append (np.nanstd(LW_CRE)); sw_cre_std_time.append (np.nanstd(SW_CRE)); net_cre_std_time.append(np.nanstd(net_CRE))
        
        # CERES extra variables:   # 
        for var in ['toa_alb_all_1h', 'toa_alb_clr_1h', 'toa_sw_all_1h', 'toa_sw_clr_1h', 'toa_lw_all_1h', 'toa_lw_clr_1h',\
                    'toa_solar_all_1h', 'cldarea_low_1h' , 'cldarea_total_1h', 'cldhght_top_low_1h', 'cldtau_low_1h', 'lwp_low_1h',\
                    'lwp_total_1h', 'iwp_low_1h', 'iwp_total_1h', 'solar_zen_angle_1h', 'cldwatrad_low_1h',\
                    'adj_atmos_sw_down_all_surface_1h', 'adj_atmos_sw_up_all_surface_1h', 'adj_atmos_lw_down_all_surface_1h',\
                    'adj_atmos_lw_up_all_surface_1h']:     
        #####################################            
            
#             if var == 'lw_cre':
#                 means_dict.setdefault(var, []).append(np.nanmean(LW_CRE))
#                 stds_dict.setdefault (var, []).append(np.nanstd (LW_CRE))             
#             elif var == 'sw_cre':
#                 means_dict.setdefault(var, []).append(np.nanmean(SW_CRE))
#                 stds_dict.setdefault (var, []).append(np.nanstd (SW_CRE))            
#             elif var == 'net_cre':
#                 means_dict.setdefault(var, []).append(np.nanmean(net_CRE))
#                 stds_dict.setdefault (var, []).append(np.nanstd (net_CRE))             
            means_dict.setdefault(var, []).append(np.nanmean(ds_sub2[var]))
            stds_dict.setdefault (var, []).append(np.nanstd (ds_sub2[var]))
                     
    for var in means_dict.keys():
        attrs = attr_dict.setdefault(var, ds_sub2[var].attrs)
        ds['CERES_'+var] = (('time'), np.array(means_dict[var]), attrs)
        attrs.update(long_name=attrs['long_name']+', standard deviation over box')
        ds['CERES_'+var+'_std'] = (('time'), np.array(stds_dict[var]), attrs)
    
    #adding ND
    nd = 1.4067 * 10**4 * (ds.CERES_cldtau_low_1h.values ** 0.5) / (ds.CERES_cldwatrad_low_1h**2.5)
    ds['CERES_Nd'] = (('time'), nd, {'long_name': 'cloud droplet number concentration',
     'units': 'cm**-3'})
 
    #####################################
    ####### Modified by Ehsan Erfani:
    ds['CERES_lw_cre' ] = (('time'), np.array(lw_cre_time ), {'long_name': 'Longwave Cloud Radiative Effect' , 'units': 'W m-2'})
    ds['CERES_sw_cre' ] = (('time'), np.array(sw_cre_time ), {'long_name': 'Shortwave Cloud Radiative Effect', 'units': 'W m-2'})
    ds['CERES_net_cre'] = (('time'), np.array(net_cre_time), {'long_name': 'Net Cloud Radiative Effect'      , 'units': 'W m-2'})

    ds['CERES_lw_cre_std' ] = (('time'), np.array(lw_cre_std_time ), {'long_name': 'Std of LW Cloud Radiative Effect'  , 'units': 'W m-2'})
    ds['CERES_sw_cre_std' ] = (('time'), np.array(sw_cre_std_time ), {'long_name': 'Std of SW Cloud Radiative Effect'  , 'units': 'W m-2'})
    ds['CERES_net_cre_std'] = (('time'), np.array(net_cre_std_time ), {'long_name': 'Std of Net Cloud Radiative Effect', 'units': 'W m-2'})
    #####################################

    ds.attrs['CERES_params'] = f'CERES data is from SYN1deg hourly product; statistics computed over a {box_degrees}-deg average centered on trajectory'
    ds.attrs['CERES_reference'] = f'CERES data available from NASA LARC at ceres.larc.nasa.gov/data.  doi: 10.1175/JTECH-D-12-00136.1, doi: 10.1175/JTECH-D-15-0147.1 \n CERES Nd derived following Painemal and Zuidema (2011), eqn 7 (doi:10.1029/2011JD016155).'
    return ds
    