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
            'stream'  : 'enda',
            'format': 'netcdf',
            'variable': ['fraction_of_cloud_cover', 
                         'relative_humidity',
                         'temperature',
                         'u_component_of_wind',
                         'v_component_of_wind',
                         'vertical_velocity',
                        ],
            'pressure_level': levels,
            'year': f'{date:%Y}',
            'month': f'{date:%m}',
            'day': f'{date:%d}',
            'time': [
                '00:00', '03:00', '06:00',
                '09:00', '12:00', '15:00',
                '18:00', '21:00',
            ],
            'area': area,
            'grid': [.5, .5]
        },
        os.path.join(saveloc, f'ERA5.enda.pres.{id_string}.{datestr}.nc'))


######
def test_api():
    date = dt.datetime(int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]))
            
    all_levels = ['1', '2', '3', '5', '7', '10', '20', '30', '50', '70', '100', '125', '150', '175', '200', '225', '250', '300', '350', '400', '450', '500', '550', '600', '650', '700', '750', '775', '800', '825', '850', '875', '900', '925', '950', '975', '1000',]
    sea_area = [60, -160, 0, -110]
    saveloc = '/home/disk/eos3/erfani/Data/Ship_tracks/ERA5/ensemble/'
    id_string='NEP'

    get_pressure_level_ERA5_Data(date, all_levels, sea_area, saveloc, id_string)


######
test_api()



