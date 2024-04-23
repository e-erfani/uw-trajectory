#############################################################################
## Modified by Ehsan Erfani (2022):
## - year, month, and day are provided as arguments to facilitate automation.
## Make sure to change the parameters based on your own need.
## Parameters: data type, variable names, domain, temporal resolution, 
##             horizontal resolution, storage path, ...
#############################################################################

import cdsapi
import datetime as dt
from ecmwfapi import ECMWFDataServer
import os
import sys
#####
def get_surface_level_ERA5_Data(date, levels, area, saveloc, id_string=''):
    datestr = dt.datetime.strftime(date, '%Y-%m-%d')
    c = cdsapi.Client()

    c.retrieve(
        'reanalysis-era5-single-levels',
        {
            'product_type': 'reanalysis',
            'format': 'netcdf',

           'variable': [
               '10m_u_component_of_wind', '10m_v_component_of_wind', '2m_dewpoint_temperature',
               '2m_temperature', 'boundary_layer_height', 'cloud_base_height',
               'convective_precipitation', 'convective_rain_rate', 'high_cloud_cover',
               'instantaneous_large_scale_surface_precipitation_fraction', 'large_scale_precipitation', 'large_scale_precipitation_fraction',
               'large_scale_rain_rate', 'low_cloud_cover', 'mean_convective_precipitation_rate',
               'mean_large_scale_precipitation_rate', 'medium_cloud_cover', 'sea_surface_temperature',
               'surface_latent_heat_flux', 'surface_pressure', 'surface_sensible_heat_flux',
               'toa_incident_solar_radiation', 'total_cloud_cover', 'total_column_cloud_liquid_water',
               'total_column_rain_water', 'total_column_water_vapour', 'total_precipitation',
              ],

            'year': f'{date:%Y}',
            'month': f'{date:%m}',
            'day': f'{date:%d}',
            'time': [
                '00:00', '01:00', '02:00',
                '03:00', '04:00', '05:00',
                '06:00', '07:00', '08:00',
                '09:00', '10:00', '11:00',
                '12:00', '13:00', '14:00',
                '15:00', '16:00', '17:00',
                '18:00', '19:00', '20:00',
                '21:00', '22:00', '23:00',
            ],
            'area': area,
        },
        os.path.join(saveloc, f'ERA5.sfc.{id_string}.{datestr}.nc'))


######
def test_api():
    date = dt.datetime(int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]))
            
    all_levels = ['1', '2', '3', '5', '7', '10', '20', '30', '50', '70', '100', '125', '150', '175', '200', '225', '250', '300', '350', '400', '450', '500', '550', '600', '650', '700', '750', '775', '800', '825', '850', '875', '900', '925', '950', '975', '1000',]
    sea_area = [60, -180, 0, -100]
    saveloc = '/home/disk/eos3/erfani/Data/NEP/ERA5/NEW/'
    id_string='NEP'

    get_surface_level_ERA5_Data(date, all_levels, sea_area, saveloc, id_string)


######
test_api()
