"""
Created Dec 2 2020
author: Hans Mohrmann (jkcm@uw.edu)

some utility functions for getting ERA5 data from the Copernicus Data Store
"""
import cdsapi
import xarray



# import cdsapi

# c = cdsapi.Client()

# c.retrieve(
#     'cams-global-reanalysis-eac4',
#     {
#         'date': '2017-07-28/2017-09-06',
#         'format': 'netcdf',
#         'variable': [
#             'carbon_monoxide', 'dust_aerosol_0.03-0.55um_mixing_ratio', 'dust_aerosol_0.55-0.9um_mixing_ratio',
#             'dust_aerosol_0.9-20um_mixing_ratio', 'fraction_of_cloud_cover', 'geopotential',
#             'hydrophilic_black_carbon_aerosol_mixing_ratio', 'hydrophilic_organic_matter_aerosol_mixing_ratio', 'hydrophobic_black_carbon_aerosol_mixing_ratio',
#             'hydrophobic_organic_matter_aerosol_mixing_ratio', 'ozone', 'relative_humidity',
#             'sea_salt_aerosol_0.03-0.5um_mixing_ratio', 'sea_salt_aerosol_0.5-5um_mixing_ratio', 'sea_salt_aerosol_5-20um_mixing_ratio',
#             'specific_cloud_liquid_water_content', 'specific_humidity', 'sulphate_aerosol_mixing_ratio',
#             'temperature', 'u_component_of_wind', 'v_component_of_wind',
#             'vertical_velocity',
#         ],
#         'pressure_level': [
#             '500', '600', '700',
#             '800', '850', '900',
#             '925', '950', '1000',
#         ],
#         'time': [
#             '00:00', '03:00', '06:00',
#             '09:00', '12:00', '15:00',
#             '18:00', '21:00',
#         ],
#         'area': [
#             10, -35, -60,
#             15,
#         ],
#     },
#     'download.nc')
