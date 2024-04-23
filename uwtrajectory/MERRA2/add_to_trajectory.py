###############################################
## Modifications by Ehsan Erfani (2022):
## - Now, the path and file name are specified
##   in config file
## - New reference is added for calculation of
##   MERRA2 mass from number concentration  
##############################################

import os
import conda

conda_file_dir = conda.__file__
conda_dir = conda_file_dir.split('lib')[0]
proj_lib = os.path.join(os.path.join(conda_dir, 'share'), 'proj')
os.environ["PROJ_LIB"] = proj_lib

from mpl_toolkits.basemap import Basemap #```

#os.chdir(r"../")
import utils, les_utils, config
import xarray as  xr
import numpy as np
#os.chdir(r"./MERRA2")


def add_MERRA_to_trajectory(ds, box_degrees=2):
    """Add MERRA-inferred aerosol number concentrations to trajectory.
    """ 
    lats, lons, times = ds.lat.values, ds.lon.values, utils.as_datetime(ds.time.values)
    unique_days = set([utils.as_datetime(i).date() for i in times])
    files = [os.path.join(config.MERRA_dir, config.MERRA_fmt.format(i))     for i in unique_days]
    ######################
    ## Modified by Ehsan Erfani to include file format for years after 2020:
    if times[0].year > 2020:
        files = [os.path.join(config.MERRA_dir, config.MERRA_fmt_20s.format(i)) for i in unique_days]
    #####################
    
#     if location=='nep':
#         files = [os.path.join("/home/disk/eos4/jkcm/Data/MERRA/3h/", "more_vertical", "MERRA2_400.inst3_3d_aer_Nv.{:%Y%m%d}.nc4.nc4".format(i))
#                 for i in unique_days]
    
#     elif location=='sea':
#         files = [os.path.join("/home/disk/eos4/jkcm/Data/MERRA/sea/new/", "MERRA2_400.inst3_3d_aer_Nv.{:%Y%m%d}.SUB.nc".format(i))
#                 for i in unique_days]

    with xr.open_mfdataset(sorted(files), combine='by_coords') as merra_data:
        
        if np.abs(np.mean(ds.lon.values))>90: # i.e. our trajectory is closer to 180 than it is to 0 lon. force to 0-360
            merra_data.coords['lon'] = merra_data['lon']%360
            lons = lons%360
        else:
            merra_data.coords['lon'] = (merra_data['lon']+180)%360-180
            lons = (lons+180)%360-180

        merra_data = merra_data.sel(lat=slice(np.min(lats)-2, np.max(lats)+2), lon=slice(np.min(lons)-2, np.max(lons)+2))

        
        # calculating MERRA pressure levels and heights
        dz = merra_data.DELP/(9.81*merra_data.AIRDENS)
        assert(merra_data.DELP.dims[1]=='lev')
        assert(dz.dims[1]=='lev')
        assert(dz.lev[-1]==72)     
        z = dz[:,::-1,:,:].cumsum(axis=1)[:,::-1,:,:]
        z.load()
        z[:, :-1, :, :] = (z.values[:, 1:, :, :]+z.values[:, :-1, :, :])/2
        z[:, -1,:,:] = z[:, -1,:,:]/2
        p = -merra_data.DELP[:,::-1,:,:].cumsum(axis=1)[:,::-1,:,:]
        p = p + merra_data.PS
        p.load()
        p[:, :-1, :, :] = (p.values[:, 1:, :, :]+p.values[:, :-1, :, :])/2
        p[:, -1,:,:] = (p[:, -1,:,:]+merra_data.PS)/2

        merra_data['H'] = z
        merra_data.H.attrs = {'long_name': 'mid_layer_heights', 'units': 'm'}
        merra_data['PL'] = p
        merra_data.PL.attrs = {'long_name': 'mid_level_pressure', 'units': 'Pa'}

        
        vals_to_add = ['Na_tot', 'Na_tot_corr', 'H', 'PL', 'RH', 'AIRDENS']
        na_tot = np.zeros_like(merra_data.SS001.values)
        
        new_vals = []
        for varname,params in les_utils.merra_species_dict_colarco.items():
            vals_to_add.append(varname)
            var = merra_data[varname]

            num=les_utils.mass_to_number(mass=var, air_density=merra_data.AIRDENS.values, shape_params=params)

            na_tot = na_tot+num
            merra_data[varname+'_Na'] = (('time', 'lev', 'lat', 'lon'), num)
            new_vals.append(varname+'_Na')
        merra_data['Na_tot'] = (('time', 'lev', 'lat', 'lon'), na_tot, {'long_name': 'total aerosol number concentration, >100 um', 'units': 'cm**-3'})

        merra_data['Na_tot_corr'] = (('time', 'lev', 'lat', 'lon'), np.exp(1.24*np.log(na_tot) + 0.18), {'long_name': 'total aerosol number concentration, >100 um, corrected to aircraft', 'units': 'cm**-3'})  
        merra_data['Na_tot_corr_BL_logfit'] = (('time', 'lev', 'lat', 'lon'), np.exp(0.63*np.log(na_tot) + 2.42), {'long_name': 'total aerosol number concentration, >100 um, corrected to aircraft (boundary layer obs only)', 'units': 'cm**-3'})  
                   
                   
        ds = ds.assign_coords(lev = les_utils.MERRA_lev(merra_data.lev))
        
        merra_data = merra_data.assign_coords(lev = les_utils.MERRA_lev(merra_data.lev))                
        
        for var in vals_to_add+new_vals:
            var_shape = merra_data[var].isel(time=0, lat=0, lon=0).shape
            vals = []
            for (lat, lon, time) in zip(lats, lons, times):
                if lat > np.max(merra_data.coords['lat']) or lat < np.min(merra_data.coords['lat']) or \
                    lon > np.max(merra_data.coords['lon']) or lon < np.min(merra_data.coords['lon']):
                    print(f'out of range of data" {lat}, {lon}, {time}')
                    print(merra_data.coords['lat'])
                    print(merra_data.coords['lon'])
                    raise ValueError()
                    continue
                try:
                    x = merra_data[var].sel(lon=slice(lon - box_degrees/2, lon + box_degrees/2),
                                            lat=slice(lat - box_degrees/2, lat + box_degrees/2))
                except KeyError as e:
                    print(var)
                    print(lon, lat)
                    print(merra_data.lon)
                    print(merra_data.lat)

                    raise e
                z = x.sel(method='nearest', time=time, tolerance=np.timedelta64(2, 'h'))
                #this applies a 2D gaussian the width of z, i.e. sigma=box_degrees
                gauss_shape = tuple([v for v,i in zip(z.shape,z.dims) if i in ['lat', 'lon'] ])
                gauss = utils.gauss2D(shape=gauss_shape, sigma=gauss_shape[-1])
                filtered = z * gauss
                vals.append(filtered.sum(dim=('lat', 'lon')).values)
            if var in vals_to_add:
                attrs = merra_data[var].attrs
            elif var in new_vals:
                attrs = {'long_name': merra_data[var[:-3]].long_name + ', inferred aerosol number concentration', 'units':'cm**-3'}
            ds['MERRA_'+var] = (tuple(x for x in merra_data[var].dims if x not in ['lat', 'lon']), np.array(vals), attrs)
    
    
    
    
    ds['MERRA_Na_tot_mass'] = ds.MERRA_OCPHILIC + ds.MERRA_OCPHOBIC + ds.MERRA_BCPHILIC + \
                              ds.MERRA_BCPHOBIC + ds.MERRA_SO4 + ds.MERRA_DU001 + ds.MERRA_DU002 +\
                              ds.MERRA_DU003 +ds.MERRA_DU004 + ds.MERRA_DU005 + ds.MERRA_SS001 + \
                              ds.MERRA_SS002 + ds.MERRA_SS003 + ds.MERRA_SS004 + ds.MERRA_SS005
   
    # #akn=aitken = everything below 80nm
    # #acc = accumulution = everything between 80 and 1000
    # #crs=coarse = everything above 1000
    mass_acc_dict = {}
    mass_akn_dict = {}
    mass_crs_dict = {}
    num_acc_dict = {}
    num_akn_dict = {}
    num_crs_dict = {}
    
    
    
    for x in ['MERRA_OCPHILIC', 'MERRA_OCPHOBIC', 'MERRA_BCPHILIC', 'MERRA_BCPHOBIC', 'MERRA_SO4']:
        params = les_utils.merra_species_dict_colarco[x.split('_')[1]]
        data = ds[x]

        rho = ds.MERRA_AIRDENS.values 
        n0 = les_utils.get_n0(mass=data.values, density=params['density'], r_max=50, r_min=0.001, 
                              std_dev=params['geometric_std_dev'], mode_radius=params['mode_radius'])
        mass_acc_dict[x] = les_utils.get_m_subset(density=params['density'], n0=n0, r_min=0.08, r_max=1, 
                                                  std_dev=params['geometric_std_dev'], mode_radius=params['mode_radius'])
        mass_akn_dict[x] = les_utils.get_m_subset(density=params['density'], n0=n0, r_min=0.01, r_max=0.08, 
                                                  std_dev=params['geometric_std_dev'], mode_radius=params['mode_radius'])
        mass_crs_dict[x] = les_utils.get_m_subset(density=params['density'], n0=n0, r_min=1, r_max=params['upper'], 
                                                  std_dev=params['geometric_std_dev'], mode_radius=params['mode_radius'])
        num_acc_dict[x] = les_utils.get_n_subset(n0, r_min=0.08, r_max=1, 
                                                 std_dev=params['geometric_std_dev'], mode_radius=params['mode_radius'])
        num_akn_dict[x] = les_utils.get_n_subset(n0, r_min=0.01, r_max=0.08, 
                                                 std_dev=params['geometric_std_dev'], mode_radius=params['mode_radius'])
        num_crs_dict[x] = les_utils.get_n_subset(n0, r_min=1, r_max=params['upper'], 
                                                 std_dev=params['geometric_std_dev'], mode_radius=params['mode_radius'])
    
        ds[x+'_n0'] = (('time', 'lev'), n0*rho)


    mass_acc_attrs = {'long_name': 'accumulation mode aerosol mass',
                        'units': 'kg kg**-1'}
    mass_akn_attrs = {'long_name': 'aikten mode aerosol mass',
                        'units': 'kg kg**-1'}
    mass_crs_attrs = {'long_name': 'coarse mode aerosol mass',
                        'units': 'kg kg**-1'}

    ds['MERRA_acc_mass'] = (('time', 'lev'), np.sum(list(mass_acc_dict.values()), axis=0) + \
                                ds.MERRA_DU001.values + ds.MERRA_SS002.values + ds.MERRA_SS003.values,
                                mass_acc_attrs)
    ds['MERRA_akn_mass'] = (('time', 'lev'), np.sum(list(mass_akn_dict.values()), axis=0) + \
                                ds.MERRA_SS001.values,
                                mass_akn_attrs)
    ds['MERRA_crs_mass'] = (('time', 'lev'), np.sum(list(mass_crs_dict.values()), axis=0) + \
                                ds.MERRA_DU002.values + ds.MERRA_DU003.values + ds.MERRA_DU004.values + ds.MERRA_DU005.values + \
                                ds.MERRA_SS004.values + ds.MERRA_SS005.values,
                                mass_crs_attrs)
    
    num_acc_attrs = {'long_name': 'accumulation mode aerosol number',
                        'units': 'kg kg**-1'}
    num_akn_attrs = {'long_name': 'aikten mode aerosol number',
                        'units': 'kg kg**-1'}
    num_crs_attrs = {'long_name': 'coarse mode aerosol number',
                        'units': 'kg kg**-1'}
    
    ds['MERRA_acc_num'] = (('time', 'lev'), np.sum(list(num_acc_dict.values()), axis=0) + \
                                ds.MERRA_DU001_Na.values + ds.MERRA_SS002_Na.values + ds.MERRA_SS003_Na.values,
                                mass_acc_attrs)
    
    
    ds.lev.attrs['long_name'] = 'model level pressure'
    ds.lev.attrs['units'] = 'millibars'
    
    ds.attrs['MERRA_params'] = f'MERRA-2 data primarily downloaded from NASA GMAO, and statistics computed over a {box_degrees}-deg average centered on trajectory. For aerosol estimates (Na), equivalent aerosol number is computed based on aerosol mass consistent with the MERRA2-assumed aerosol optical properties. For more details, see: Erfani et al. (2022) Simulating aerosol lifecycle impacts on the subtropical stratocumulus-to-cumulus transition using large eddy simulations, J. Geophys. Res. Atmos. https://doi.org/10.1029/2022JD037258'
    ds.attrs['MERRA_reference'] = 'MERRA-2 data available at https://gmao.gsfc.nasa.gov/reanalysis/MERRA-2.'
    
    return ds


def new_add_MERRA_to_trajectory(ds, box_degrees=2):
    
    ds['MERRA_Na_tot_mass'] = ds.MERRA_OCPHILIC + ds.MERRA_OCPHOBIC + ds.MERRA_BCPHILIC + ds.MERRA_BCPHOBIC + ds.MERRA_SO4 + ds.MERRA_DU001 + ds.MERRA_DU002 + ds.MERRA_DU003 +ds.MERRA_DU004 + \
                       ds.MERRA_DU005 + ds.MERRA_SS001 + ds.MERRA_SS002 + ds.MERRA_SS003 + ds.MERRA_SS004 + ds.MERRA_SS005
    
    #akn=aitken = everything below 80nm
    #acc = accumulution = everything between 80 and 1000
    #crs=coarse = everything above 1000
    
    mass_acc_dict = {}
    mass_aik_dict = {}
    mass_crs_dict = {}
    num_acc_dict = {}
    num_aik_dict = {}
    num_crs_dict = {}
    
    # for x in ['MERRA_OCPHILIC', 'MERRA_OCPHOBIC', 'MERRA_BCPHILIC', 'MERRA_BCPHOBIC', 'MERRA_SO4']:
        
    
    # ds['MERRA_Na_acc_mass'] = 
    # ds['MERRA_Na_akn_mass'] = 
    # ds['MERRA_Na_crs_mass'] = 
    
    # ds['MERRA_akn_num'] = 
    # ds['MERRA_acc_num'] = 
    # ds['MERRA_crs_num'] = 
    
    
    
    
    return ds

    #aitken = low-100nm
    #add aikten NUMBER
    #add aikten mass
    