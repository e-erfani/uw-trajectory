



#need to acquire surface data, pres data, and ens data
import cdsapi
import datetime as dt
from ecmwfapi import ECMWFDataServer
import os


def test_api():
    date = dt.datetime(2015, 4, 1)
    bl_levels = "700/750/775/800/825/850/875/900/925/950/975/1000"
#     get_pressure_level_ERA5_Data(date, bl_levels)
    
    
    
    all_levels = ['1', '2', '3', '5', '7', '10', '20', '30', '50', '70', '100', '125', '150', '175', '200', '225', '250', '300', '350', '400', '450', '500', '550', '600', '650', '700', '750', '775', '800', '825', '850', '875', '900', '925', '950', '975', '1000',]
    sea_area = [60, -160, 0, -110]
    saveloc = '/home/disk/eos4/jkcm/Data/ORACLES/ERA5/'
    id_string='SEA'

    get_pressure_level_ERA5_Data(date, all_levels, sea_area, saveloc, id_string)
    
    
def old_get_pressure_level_ERA5_Data(date, levels):
    datestr = dt.datetime.strftime(date, '%Y-%m-%d')

    server = ECMWFDataServer()
    server.retrieve({
        "class": "ea",  # Do not change
        "dataset": "era5",  # Do not change
        "expver": "1",  # Do not change
        "stream": "oper",
    # can be "oper", "wave", etcetera; see ERA5 catalogue (http://apps.ecmwf.int/data-catalogues/era5 ) and ERA5 documentation (https://software.ecmwf.int/wiki/display/CKB/ERA5+data+documentation )
        "type": "an",  # can be an (Analysis) or fc (forecast) or 4v (4D variational analysis)
        "levtype": "pl",  # can be "sfc", "pl", "ml", etcetera; see ERA5 documentation
        "param": "u/v/w/r/z/t/o3",
    # Parameters you want to retrieve. For available parameters see the ERA5 documentation. Specify here using shortName or paramID, and separated by '/'.
        "levelist": levels,
        "date": datestr,  # Set a single date as "YYYY-MM-DD" or a range as "YYYY-MM-DD/to/YYYY-MM-DD".
        "time": "00:00:00/01:00:00/02:00:00/03:00:00/04:00:00/05:00:00/06:00:00/07:00:00/08:00:00/09:00:00/10:00:00/11:00:00/12:00:00/13:00:00/14:00:00/15:00:00/16:00:00/17:00:00/18:00:00/19:00:00/20:00:00/21:00:00/22:00:00/23:00:00",
    # If above you set "type":"an", "time" is the time of analysis. If above you set "type":"fc", "time" is the initialisation time of the forecast.
        "step": "0",
    # The forecast step. If above you set "type":"an", set "step":"0". If above you set "type":"fc", set "step" > 0.
        "grid": "0.3/0.3",
    # Optional. The horizontal resolution in decimal degrees. If not set, the archived grid as specified in the data documentation is used.
        "area": "45/-160/15/-115",
    # Optional. Subset (clip) to an area. Specify as N/W/S/E in Geographic lat/long degrees. Southern latitudes and western longitudes must be
        # given as negative numbers. Requires "grid" to be set to a regular grid, e.g. "0.3/0.3".
        "format": "netcdf",
    # Optional. Output in NetCDF format. Requires that you also specify 'grid'. If not set, data is delivered in GRIB format, as archived.
        "target": "/home/disk/eos4/jkcm/Data/CSET/ERA5/ERA5.pres.NEP.{}.nc".format(datestr),
    # Change this to the desired output path and file name, e.g. "data1.nc" or "./data/data1.grib". The default path is the current working directory.
    })


levels = ['1', '2', '3', '5', '7', '10', '20', '30', '50', '70', '100', '125', '150', '175', '200', '225', '250', '300', '350', '400', '450', '500', '550', '600', '650', '700', '750', '775', '800', '825', '850', '875', '900', '925', '950', '975', '1000',]
sea_area = [60, -160, 0, -110]

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






# import cdsapi

# c = cdsapi.Client()

# c.retrieve(
#     'reanalysis-era5-single-levels',
#     {
#         'product_type': 'reanalysis',
#         'variable': [
#             '10m_u_component_of_wind', '10m_v_component_of_wind', '2m_dewpoint_temperature',
#             '2m_temperature', 'boundary_layer_height', 'cloud_base_height',
#             'convective_precipitation', 'convective_rain_rate', 'high_cloud_cover',
#             'instantaneous_large_scale_surface_precipitation_fraction', 'large_scale_precipitation', 'large_scale_precipitation_fraction',
#             'large_scale_rain_rate', 'low_cloud_cover', 'mean_convective_precipitation_rate',
#             'mean_large_scale_precipitation_rate', 'medium_cloud_cover', 'sea_surface_temperature',
#             'surface_latent_heat_flux', 'surface_pressure', 'surface_sensible_heat_flux',
#             'toa_incident_solar_radiation', 'total_cloud_cover', 'total_column_cloud_liquid_water',
#             'total_column_rain_water', 'total_column_water_vapour', 'total_precipitation',
#         ],
#         'year': '2015',
#         'month': '07',
#         'day': [
#             '01', '02', '03',
#             '04', '05', '06',
#             '07', '08', '09',
#             '10', '11', '12',
#             '13', '14', '15',
#             '16', '17', '18',
#             '19', '20', '21',
#             '22', '23', '24',
#             '25', '26', '27',
#             '28', '29', '30',
#             '31',
#         ],
#         'time': [
#             '00:00', '01:00', '02:00',
#             '03:00', '04:00', '05:00',
#             '06:00', '07:00', '08:00',
#             '09:00', '10:00', '11:00',
#             '12:00', '13:00', '14:00',
#             '15:00', '16:00', '17:00',
#             '18:00', '19:00', '20:00',
#             '21:00', '22:00', '23:00',
#         ],
#         'area': [
#             60, -160, 0,
#             -110,
#         ],
#         'format': 'netcdf',
#     },
#     'download.nc')

test_api()