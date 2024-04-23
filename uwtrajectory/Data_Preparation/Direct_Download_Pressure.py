#############################################################################
## Modified by Ehsan Erfani (2022):
## - year, month, and day are provided as arguments to facilitate automation.
## - Ice variables are added.
## Make sure to change the parameters based on your own need.
## Parameters: data type, variable names, domain, temporal resolution, 
##             horizontal resolution, levels, storage path, ...
#############################################################################

import cdsapi
import datetime as dt
from ecmwfapi import ECMWFDataServer
import os
import sys
#####
def get_pressure_level_ERA5_Data(date, levels, area, saveloc, id_string=''):
    datestr = dt.datetime.strftime(date, '%Y-%m-%d')
    c = cdsapi.Client()

    c.retrieve(
        'reanalysis-era5-pressure-levels',
        {
            'product_type': 'reanalysis',
            'format': 'netcdf',
            'variable': ['fraction_of_cloud_cover', 
                         'geopotential', 
                         'ozone_mass_mixing_ratio', 
                         'relative_humidity',
                         'specific_cloud_liquid_water_content',
                         'specific_rain_water_content',
                         'temperature',
                         'u_component_of_wind',
                         'v_component_of_wind',
                         'vertical_velocity',
                         'specific_cloud_ice_water_content',
                         'specific_snow_water_content',
                        ],
            'pressure_level': levels,
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
        os.path.join(saveloc, f'ERA5.pres.{id_string}.{datestr}.nc'))


######
def test_api():
    date = dt.datetime(int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]))
            
    all_levels = ['1', '2', '3', '5', '7', '10', '20', '30', '50', '70', '100', '125', '150', '175', '200', '225', '250', '300', '350', '400', '450', '500', '550', '600', '650', '700', '750', '775', '800', '825', '850', '875', '900', '925', '950', '975', '1000',]
    sea_area = [60, -180, 0, -100]
    saveloc = '/home/disk/eos3/erfani/Data/NEP/ERA5/NEW/'
    id_string='NEP'

    get_pressure_level_ERA5_Data(date, all_levels, sea_area, saveloc, id_string)


######
test_api()
