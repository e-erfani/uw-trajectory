###############################################
## Modifications by Ehsan Erfani (2022):
## - Now, the path and file name are specified
##   in config file
## - ice cloud water content and ice snow water
##   content are now added to trajectory files.
## - Most of print commands are commented for
##   better run during numeroucs cases.
##############################################

import numpy as np
import xarray as xr
import os

import utils, met_utils, config
from LoopTimer import LoopTimer

def add_ERA_ens_to_trajectory(ds, box_degrees=2):
    
    lats, lons, times = ds.lat.values, ds.lon.values, ds.time.values
    space_index = int(np.round(box_degrees/0.25/2)) # go up/down/left/right this many pixels
    
    
    unique_days = set([utils.as_datetime(i).date() for i in times])
    ens_files = [os.path.join(config.ERA_ens_source, config.ERA_ens_fmt.format(i)) for i in unique_days]
    

    with xr.open_mfdataset(sorted(ens_files), combine='by_coords') as data:
#         data = data.rename({'level': 'ens_level'})
        ds.coords['number'] = data.coords['number']
#         ds.coords['ens_level'] = data.coords['ens_level']
        
                    
        if 'w' in data.data_vars.keys() and 'sp' in data.data_vars.keys():
            data['dspdt'] = (data.sp.dims, np.gradient(data.sp, np.median(np.gradient(data.time.values)/np.timedelta64(1, 's')), axis=0),
                    {'units': "Pa s**-1", 'long_name': "Surface pressure tendency", 'standard_name': 'tendency_of_surface_air_pressure'})
            data['w_corr'] = (data.w.dims, data.w - data.dspdt, {'units': data.w.units, 'long_name': 'Vertical velocity (sp-corrected)'})
        
        
        for var in data.data_vars.keys():
            var_shape = data[var].isel(time=0, latitude=0, longitude=0).shape
            vals = []
            for (lat, lon, time) in zip(lats, lons, times):
                if lat > np.max(data.coords['latitude']) or lat < np.min(data.coords['latitude']) or \
                    lon > np.max(data.coords['longitude']) or lon < np.min(data.coords['longitude']):
                    print(f'out of range of data" {lat}, {lon}, {time}')
                    vals.append(np.full(var_shape, float('nan'), dtype='float'))
                    continue
                x = data[var].sel(longitude=slice(lon - box_degrees/2, lon + box_degrees/2),
                                  latitude=slice(lat + box_degrees/2, lat - box_degrees/2))
                z = x.sel(method='nearest', time=time, tolerance=np.timedelta64(1, 'h'))
                #this applies a 2D gaussian the width of z, i.e. sigma=box_degrees
                gauss_shape = tuple([v for v,i in zip(z.shape,z.dims) if i in ['latitude', 'longitude'] ])
                gauss = utils.gauss2D(shape=gauss_shape, sigma=gauss_shape[-1])
                filtered = z * gauss
                vals.append(filtered.sum(dim=('latitude', 'longitude')).values)
            ds['ERA_ens_'+var] = (tuple(x for x in data[var].dims if x not in ['latitude', 'longitude']), np.array(vals), data[var].attrs)

#     print('adding ensemble temperatures...')
#     ens_temp_files = [os.path.join(utils.ERA_ens_temp_source, i) for i in sorted(os.listdir(utils.ERA_ens_temp_source))]
    
#     with xr.open_mfdataset(sorted(ens_temp_files), combine='by_coords') as data:
#         data = data.rename({'level': 'ens_level'})
#         #ds.coords['number'] = data.coords['number']
#         #ds.coords['ens_level'] = data.coords['ens_level']
        
#         for var in data.data_vars.keys():
#             var_shape = data[var].isel(time=0, latitude=0, longitude=0).shape
#             vals = []
#             for (lat, lon, time) in zip(lats, lons, times):
#                 if lat > np.max(data.coords['latitude']) or lat < np.min(data.coords['latitude']) or \
#                     lon > np.max(data.coords['longitude']) or lon < np.min(data.coords['longitude']):
#                     print('out of range of data')
#                     print(lat, lon, time)
#                     print(np.max(data.coords['latitude'].values), np.min(data.coords['latitude'].values))
#                     print(np.max(data.coords['longitude'].values), np.min(data.coords['longitude'].values))
#                     vals.append(np.full(var_shape, float('nan'), dtype='float'))
#                     continue
#                 x = data[var].sel(longitude=slice(lon - box_degrees/2, lon + box_degrees/2),
#                                   latitude=slice(lat + box_degrees/2, lat - box_degrees/2))
#                 z = x.sel(method='nearest', time=time, tolerance=np.timedelta64(1, 'h'))
#                 #this applies a 2D gaussian the width of z, i.e. sigma=box_degrees
#                 gauss_shape = tuple([v for v,i in zip(z.shape,z.dims) if i in ['latitude', 'longitude'] ])
#                 gauss = utils.gauss2D(shape=gauss_shape, sigma=gauss_shape[-1])
#                 filtered = z * gauss
#                 vals.append(filtered.sum(dim=('latitude', 'longitude')).values)
#             ds['ERA_ens_'+var] = (tuple(x for x in data[var].dims if x not in ['latitude', 'longitude']), np.array(vals), data[var].attrs)
    return ds


def add_ERA_to_trajectory(ds, box_degrees=2):
    """Retrieve ERA5 data in a box around a trajectory
    Assumes ERA5 data is 0.3x0.3 degrees
    Returns an xarray Dataset
    """
    lats, lons, times = ds.lat.values, ds.lon.values, ds.time.values
    space_index = int(np.round(box_degrees/0.3/2)) # go up/down/left/right this many pixels
    unique_days = set([utils.as_datetime(i).date() for i in times])
    files = [os.path.join(config.ERA_source, config.ERA_fmt.format(i))
             for i in unique_days]
#     flux_files = [os.path.join(utils.ERA_source, "ERA5.flux.NEP.{:%Y-%m-%d}.nc".format(i))
#          for i in unique_days]
    with xr.open_mfdataset(sorted(files), combine='by_coords') as full_dataset:
        data = full_dataset.sel(latitude=slice(np.max(lats)+box_degrees, np.min(lats)-box_degrees), longitude=slice(np.min(lons)-box_degrees, np.max(lons)+box_degrees))
#         print('prepping data...') 
        ds.coords['level'] = data.coords['level']
        
        #adding in q:
        T = data['t'].values 
        RH = data['r'].values
        p = np.broadcast_to(data.coords['level'].values[None, :, None, None], T.shape)*100
        q = utils.qv_from_p_T_RH(p, T, RH)
        data['q'] = (('time', 'level', 'latitude', 'longitude'), q)
        data['q'] = data['q'].assign_attrs({'units': "kg kg**-1", 
                                'long_name': "specific_humidity",
                                'dependencies': 'ERA_t, ERA_p, ERA_r'})
        MR = q/(1-q)
        data['MR'] = (('time', 'level', 'latitude', 'longitude'), MR)
        data['MR'] = data['MR'].assign_attrs({'units': "kg kg**-1", 
                                'long_name': "mixing_ratio",
                                'dependencies': 'ERA_t, ERA_p, ERA_r'})
        
        # adding gradients in for z, t, and q. Assuming constant grid spacing.
        for var in ['t', 'q', 'z', 'u', 'v', 'MR']:
#             print(f'creating gradient in {var}...', end="\r")
            [_,_,dvardj, dvardi] = np.gradient(data[var].values)
            dlatdy = 360/4.000786e7  # degrees lat per meter y
            def get_dlondx(lat) : return(360/(np.cos(np.deg2rad(lat))*4.0075017e7))

            lat_spaces = np.diff(data.coords['latitude'].values)
            lon_spaces = np.diff(data.coords['longitude'].values)
            try:
                assert(np.allclose(lat_spaces, -0.25, atol=0.01) and np.allclose(lon_spaces, 0.25, atol=0.05))
            except AssertionError as e:
                print(np.unique(lat_spaces))
                print(np.unique(lon_spaces))
                raise e
            dlondi = np.mean(lon_spaces)
            dlatdj = np.mean(lat_spaces)
            dlondx = get_dlondx(data.coords['latitude'].values)
            dvardx = dvardi/dlondi*dlondx[None,None,:,None]
            dvardy = dvardj/dlatdj*dlatdy
            data['d{}dx'.format(var)] = (('time', 'level', 'latitude', 'longitude'), dvardx)
            data['d{}dy'.format(var)] = (('time', 'level', 'latitude', 'longitude'), dvardy)

        grad_attrs = {'q': {'units': "kg kg**-1 m**-1",
                            'long_name': "{}_gradient_of_specific_humidity",
                            'dependencies': "ERA_t, ERA_p, ERA_r"},
                      't':  {'units': "K m**-1",
                            'long_name': "{}_gradient_of_temperature",
                            'dependencies': "ERA_t"},
                      'z':  {'units': "m**2 s**-2 m**-1",
                            'long_name': "{}_gradient_of_geopotential",
                            'dependencies': "ERA_z"},
                      'u': {'units': "m s**-1 m**-1",
                            'long_name': "{}_gradient_of_zonal_wind",
                            'dependencies': "ERA_u"},
                      'v': {'units': "m s**-1 m**-1",
                            'long_name': "{}_gradient_of_meridional_wind",
                            'dependencies': "ERA_v"},
                      'MR': {'units': "kg kg**-1 m**-1",
                            'long_name': "{}_gradient_of_mixing_ratio",
                            'dependencies': "ERA_t, ERA_p, ERA_r"}}

        for key, val in grad_attrs.items():
            for (n, drn) in [('x', 'eastward'), ('y', 'northward')]:
                attrs = val.copy()
                var = 'd{}d{}'.format(key, n)
                attrs['long_name'] = attrs['long_name'].format(drn)
                data[var] = data[var].assign_attrs(attrs)
#         print('\ndone')
        for var in data.data_vars.keys():
            vals = []
#             print(f'working on {var}...')
            lt = LoopTimer(len(lats))
            for i, (lat, lon, time) in enumerate(zip(lats, lons, times)):
                lt.update()
                if lat > np.max(data.coords['latitude']) or lat < np.min(data.coords['latitude']) or \
                    lon > np.max(data.coords['longitude']) or lon < np.min(data.coords['longitude']):
                    print(f'out of range of data" {lat}, {lon}, {time}')
                    print(np.max(data.coords['latitude']), np.min(data.coords['latitude']), np.max(data.coords['longitude']) , np.min(data.coords['longitude']))
                    vals.append(np.full_like(data.coords['level'], float('nan'), dtype='float'))
                    continue
                x = data[var].sel(longitude=slice(lon - box_degrees/2, lon + box_degrees/2),
                                  latitude=slice(lat + box_degrees/2, lat - box_degrees/2))

                z = x.sel(method='nearest', tolerance=np.timedelta64(1, 'h'), time=time)
                #this applies a 2D gaussian the width of z, i.e. sigma=box_degrees
                gauss = utils.gauss2D(shape=z.shape[1:], sigma=z.shape[0])
                filtered = z.values * gauss
                vals.append(np.sum(filtered, axis=(1,2)))
#             print('')
            ds['ERA_'+var] = (('time', 'level'), np.array(vals))
            ds['ERA_'+var] = ds['ERA_'+var].assign_attrs(data[var].attrs)
             
        t_1000 = ds.ERA_t.sel(level=1000).values
        theta_700 = met_utils.theta_from_p_T(p=700, T=ds.ERA_t.sel(level=700).values)
        LTS = theta_700-t_1000
        ds['ERA_LTS'] = (('time'), np.array(LTS))
        ds['ERA_LTS'] = ds['ERA_LTS'].assign_attrs(
                {"long_name": "Lower tropospheric stability",
                 "units": "K",
                 "_FillValue": "NaN"})
        t_dew = t_1000-(100-ds.ERA_r.sel(level=1000).values)/5
        lcl = met_utils.calculate_LCL(t=t_1000, t_dew=t_dew, z=ds.ERA_z.sel(level=1000).values/9.81)
        z_700 = ds.ERA_z.sel(level=700).values/9.81
        gamma_850 = met_utils.get_moist_adiabatic_lapse_rate(ds.ERA_t.sel(level=850).values, 850)
        eis = LTS - gamma_850*(z_700-lcl)
        ds['ERA_EIS'] = (('time'), np.array(eis))
        ds['ERA_EIS'] = ds['ERA_EIS'].assign_attrs(
                {"long_name": "Estimated inversion strength",
                 "units": "K",
                 "_FillValue": "NaN"})
        
#         with xr.open_mfdataset(sorted(flux_files), combine='by_coords') as flux_data:
#             for var in flux_data.data_vars.keys():
#                 vals = []
#                 for (lat, lon, time) in zip(lats, lons%360, times):
#                     if lat > np.max(flux_data.coords['latitude']) or lat < np.min(flux_data.coords['latitude']) or \
#                         lon > np.max(flux_data.coords['longitude']) or lon < np.min(flux_data.coords['longitude']):
#                         print(f'out of range of data" {lat}, {lon}, {time}')
#                         vals.append(float('nan'))
#                         continue
#                     x = flux_data[var].sel(longitude=slice(lon - box_degrees/2, lon + box_degrees/2),
#                                           latitude=slice(lat + box_degrees/2, lat - box_degrees/2))
#                     z = x.sel(method='nearest', time=time, tolerance=np.timedelta64(1, 'h'))
#                     gauss = utils.gauss2D(shape=z.shape, sigma=z.shape[0])
#                     filtered = z.values * gauss
#                     vals.append(np.sum(filtered))
#                 ds['ERA_'+var] = (('time'), np.array(vals))
#                 ds['ERA_'+var] = ds['ERA_'+var].assign_attrs(flux_data[var].attrs)

    ds.attrs['ERA_params'] = f'ERA5 data acquired from ECWMF Copernicus at cds.climate.copernicus.eu/. statistics computed over a {box_degrees}-deg average centered on trajectory. EIS and LTS computed according to Wood and Bretherton (2006) and Klein and Hartmann (1993) respectively.'
    ds.attrs['ERA_reference'] = 'Copernicus Climate Change Service (C3S) (2017): ERA5: Fifth generation of ECMWF atmospheric reanalyses of the global climate . Copernicus Climate Change Service Climate Data Store (CDS), date of access. https://cds.climate.copernicus.eu/cdsapp#!/home'
         
    return ds



    
def add_ERA_sfc_to_trajectory(ds, box_degrees=2):
    
    lats, lons, times = ds.lat.values, ds.lon.values, ds.time.values
    space_index = int(np.round(box_degrees/0.3/2)) # go up/down/left/right this many pixels
    unique_days = set([utils.as_datetime(i).date() for i in times])
    sfc_files = [os.path.join(config.ERA_source, config.ERA_sfc_fmt.format(i))
             for i in unique_days]
    with xr.open_mfdataset(sorted(sfc_files), combine='by_coords') as data:
         for var in data.data_vars.keys():
            vals = []
            for (lat, lon, time) in zip(lats, lons, times):
                if lat > np.max(data.coords['latitude']) or lat < np.min(data.coords['latitude']) or \
                    lon > np.max(data.coords['longitude']) or lon < np.min(data.coords['longitude']):
                    print('out of range of data')
                    print(lat, lon, time)
                    vals.append(float('nan'))
                    continue
                x = data[var].sel(longitude=slice(lon - box_degrees/2, lon + box_degrees/2),
                                  latitude=slice(lat + box_degrees/2, lat - box_degrees/2))
                z = x.sel(method='nearest', tolerance=np.timedelta64(59, 'm'), time=time)
                gauss = utils.gauss2D(shape=z.shape, sigma=z.shape[0])
                filtered = z.values * gauss
                vals.append(np.sum(filtered))
            ds['ERA_'+var] = (('time'), np.array(vals))
            ds['ERA_'+var] = ds['ERA_'+var].assign_attrs(data[var].attrs)
    return ds
    
    
    
def add_advection_to_trajectory(ds):
    
    """Add advection to trajectory after adding ERA data
    TODO make data dependencies explicit
    """
    names = dict(u='ERA_u', v='ERA_v', u_t='traj_u', v_t='traj_v',
                   dtdx='ERA_dtdx', dtdy='ERA_dtdy', dqdx='ERA_dqdx', dqdy='ERA_dqdy', dMRdx='ERA_dMRdx', dMRdy='ERA_dMRdy')
    assert np.all([i in ds.data_vars.keys() for i in names.values()])
    rel_adv_of_T = -((ds[names['u']].values-ds[names['u_t']].values[:, None])*ds[names['dtdx']].values + \
                   (ds[names['v']].values-ds[names['v_t']].values[:, None])*ds[names['dtdy']].values)
    rel_adv_of_q = -((ds[names['u']].values-ds[names['u_t']].values[:, None])*ds[names['dqdx']].values + \
                   (ds[names['v']].values-ds[names['v_t']].values[:, None])*ds[names['dqdy']].values)
    rel_adv_of_MR = -((ds[names['u']].values-ds[names['u_t']].values[:, None])*ds[names['dMRdx']].values + \
                   (ds[names['v']].values-ds[names['v_t']].values[:, None])*ds[names['dMRdy']].values)
    T_adv_attr = {'units': "K s**-1", 
                  'long_name': "trajectory_relative_advection_of_temperature",
                  'dependencies': 'ERA_t, traj_u, traj_v, ERA_u, ERA_v'}
    q_adv_attr = {'units': "kg kg**-1 s**-1", 
                  'long_name': "trajectory_relative_advection_of_specific_humidity",
                  'dependencies': 'ERA_q, traj_u, traj_v, ERA_u, ERA_v'}
    MR_adv_attr = {'units': "kg kg**-1 s**-1", 
                  'long_name': "trajectory_relative_advection_of_mixing ratio",
                  'dependencies': 'ERA_q, traj_u, traj_v, ERA_u, ERA_v'}
        
    ds['ERA_T_adv'] = (('time', 'level'), rel_adv_of_T)
    ds['ERA_T_adv'] = ds['ERA_T_adv'].assign_attrs(**T_adv_attr)
    
    ds['ERA_q_adv'] = (('time', 'level'), rel_adv_of_q)
    ds['ERA_q_adv'] = ds['ERA_q_adv'].assign_attrs(**q_adv_attr)
    
    ds['ERA_MR_adv'] = (('time', 'level'), rel_adv_of_MR)
    ds['ERA_MR_adv'] = ds['ERA_MR_adv'].assign_attrs(**MR_adv_attr)
    return ds
