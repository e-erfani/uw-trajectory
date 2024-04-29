# University of Washington Trajectory Colocated Extraction

### Authors:
This package was developed by Johannes Mohrmann and was set up for CSET 2015 field campaign. 
Ehsan Erfani modified the package for 2019 ship track cases with some other updates, outlined below.

### Original notes:
Maintainer: Johannes Mohrmann

Code repository for the running of wind-driven trajectories, based on model or reanalysis gridded winds, and the colocation and extraction of geospatial datasets along trajectories. 

Todo: all of it!

###### Notes on geographic regions:
CSET (NEP) should 0N-60N, and 160W-110W (200E-250E)
ORACLES (SEA) should be 40S-5N, 20W-10E (340E-10E)

###### Notes on data acquisition:
for ERA5, while there exists a python API, it is not instantaneous. Either the code can be made asynchronous, or (faster) a user can download the data themselves.
SSMI: downloaded from RSS - can be automated (just wget), but then transformed with classified-code. need to add code to this repo
CERES: currently downloading from larc, using data order tool

#
# 
## Notes added by Ehsan Erfani
### Main codes:

- unified_traj_data.ipynb: This is the final code that needs to be run. It uses all other subroutines and functions to read various data and create the data files along a single trajectory.
- unified_traj_data_automate.py: Similar to previous code, but it is designed to automate the process of creating numerous trajectories. 

- config.py: All parameters and file paths/names are prescribed in this code, so make sure they are correct for your datasets.

- add_to_trajectory.py: There is a code with this name in each data directory (e.g. MERRA2, ERA5, ...) that contains the details of extracting and analyzing that dataset.

#
### Note regarding the versions and libraries:

Recommend using Python3 installed by Anaconda. 

You can check the libraries in the file "conda_env.yml", but the versions of the libraries are outdated.

#
### Data acquisition and Preparation

##### Summary:

![datasets](https://github.com/e-erfani/uw-trajectory/assets/28571068/ec10c38d-b56c-40bd-b754-218bbea604ed)

A Python code is created to download all these datasets for multiple years and months [The final version to be uploaded].


#### More details are provided regarding data acquisition and preparation:

Any example code needed for this purpose is provided in this directory: [uwtrajectory/Data_Preparation](https://github.com/e-erfani/uw-trajectory/tree/main/uwtrajectory/Data_Preparation)

#### ERA5 surface and pressure:
- To create an account, log in, and create .cdsapirc file in your $HOME Linux, check this [website](https://cds.climate.copernicus.eu/api-how-to)
- To download the data, check this [website](https://confluence.ecmwf.int/display/CKB/How+to+download+ERA5)
- Install the required packages:

conda install -c conda-forge cdsapi

conda config --add channels conda-forge

conda install ecmwf-api-client

- **Note:** recommend downloading **a separate ERA5 file for each DAY** since this is the format compatible with the trajectory package.
- Example Python scripts (Direct_Download_Pressure.py, Direct_Download_Surface.py, and Direct_Download_Ensemble.py) are created that take the yyyy mm dd as arguments, submit the request to cdc, and download ERA5 data directly to the Linux server. Make sure to modify the script for your own path, data type, variables, levels, and region. Also, invoke the command with & to run in the background (it might take hours or more to get one daily file): 

nohup ipython [Direct_Download_Pressure.py](https://github.com/e-erfani/uw-trajectory/blob/main/uwtrajectory/Data_Preparation/Direct_Download_Pressure.py) 2019 6 10 &

nohup ipython [Direct_Download_Surface.py](https://github.com/e-erfani/uw-trajectory/blob/main/uwtrajectory/Data_Preparation/Direct_Download_Surface.py) 2019 6 10 &


#### MERRA2 aerosols:
- **Note:** recommend downloading **a separate MERRA2 file for each DAY** since this is the format compatible with the trajectory package.
- The aerosol data we need is called: M2I3NVAER and is [available here](https://disc.gsfc.nasa.gov/datasets/M2I3NVAER_5.12.4/summary)
- Click on Subset, then select "Get File Subsets using OPeNDAP", then select the date range and region (0 to 65; -180 to -110), and click "Get Data". 
- After the data is generated, download the list of files, transfer it to Linux and run it using a command similar to this:

wget --load-cookies ~/.urs_cookies --save-cookies ~/.urs_cookies --auth-no-challenge=on --keep-session-cookies --content-disposition -i subset_M2I3NVAER_5.12.4_20220130_024125.txt &
- First-time setup: registration, linking the account, and a few more steps for wget configuration are needed (all instructions are provided during the data generation)


#### CERES:
- **Note:** recommend downloading **a separate CERES file for each MONTH** since this is the format compatible with the trajectory package.
- Go to NASA LARC website to download [CERES SYN 1-deg data](https://ceres-tool.larc.nasa.gov/ord-tool/jsp/SYN1degEd41Selection.jsp)
- Select variables (below), levels (surface, TOA, ...), all-sky and/or clear-sky, temporal resolution (hourly), region (-65 to 65; 180 to 320), and date range (preferably a full month).
- Note: the current code is setup for these variables: 

Observed TOA Fluxes:

toa_sw_clr_1h, toa_sw_all_1h, toa_lw_clr_1h, toa_lw_all_1h, toa_wn_clr_1h, toa_wn_all_1h, toa_net_clr_1h, toa_net_all_1h, toa_alb_clr_1h, toa_alb_all_1h, toa_solar_all_1h, 

Observed cloud variables:

cldarea_low_1h, cldarea_total_1h, cldtau_low_1h, cldtau_total_1h, cldtau_lin_low_1h, cldtau_lin_total_1h, lwp_low_1h, lwp_total_1h, iwp_low_1h, iwp_total_1h, cldhght_top_low_1h, cldhght_top_total_1h, cldwatrad_low_1h, cldwatrad_total_1h, 

Number of observations:

toa_sw_num_obs_all_1h, toa_lw_num_obs_all_1h, 

Adjusted All sky fluxes:

adj_atmos_sw_up_all_toa_1h, adj_atmos_sw_up_all_surface_1h, adj_atmos_sw_down_all_toa_1h, adj_atmos_sw_down_all_surface_1h, adj_atmos_lw_up_all_toa_1h, adj_atmos_lw_up_all_surface_1h, adj_atmos_lw_down_all_toa_1h, adj_atmos_lw_down_all_surface_1h,

Adjusted Clear sky fluxes:

adj_atmos_sw_up_clr_toa_1h, adj_atmos_sw_up_clr_surface_1h, adj_atmos_sw_down_clr_toa_1h, adj_atmos_sw_down_clr_surface_1h, adj_atmos_lw_up_clr_toa_1h, adj_atmos_lw_up_clr_surface_1h, adj_atmos_lw_down_clr_toa_1h, adj_atmos_lw_down_clr_surface_1h, 

Surface Data:

solar_zen_angle_1h, aux_ocean_1h

- Order the data (registration and/or login required). It will take hours for the order to get ready.
- The download link will be sent via email. Use wget to get the data. An example command: 

wget https://ceres-tool.larc.nasa.gov/ord-tool/data1//CERES_2022-02-02:21717/dir1/CERES_SYN1deg-1H_Terra-Aqua-MODIS_Ed4.1_Subset_20190401-20190425.nc

- Possibly it is needed to concatenate multiple files to create one file for each month. A command similar to this should work:

cdo mergetime [[INPUT1]].nc [[INPUT2]].nc [[OUTPUT]].nc


#### AMSR:
- This example of the Linux Bash command downloads all the daily files for April 2019 (Website valid as of 1/26/2022):

wget -r --no-parent --reject "index.html*" https://data.remss.com/amsr2/bmaps_v08/y2019/m04/

or 

wget https://data.remss.com/amsr2/ocean/L3/v08.2/daily/2019/RSS_AMSR2_ocean_L3_daily_2019-04-{01..30}_v08.2.nc  (not used in this study)

- The original data format is binary. Use this script [RSS_bytemaps_to_netcdf.ipynb](https://github.com/e-erfani/uw-trajectory/blob/main/uwtrajectory/Data_Preparation/RSS_bytemaps_to_netcdf.ipynb) to unify the data and save in NetCDF format.
- **Note:** This script is prepared to convert data for **a full month**, so make sure data is available for a specific month, and use the script iteratively if you need to convert more than one month of data. The output of the script will be **a separate file for each month**.


#### SSMI:
- This example of the Linux Bash command downloads all the F18 daily binary files for April 2019 (Website valid as of 1/26/2022) 
- Note that there are three SSMI instruments: F16, F17, and F18, and different commands are needed for each instrument. 
- Also, use bmaps_v07 for F16 and F17):

wget -r --no-parent --reject "index.html*" https://data.remss.com/ssmi/f18/bmaps_v08/y2019/m04/

- The original data format is binary. Use this script [RSS_bytemaps_to_netcdf.ipynb](https://github.com/e-erfani/uw-trajectory/blob/main/uwtrajectory/Data_Preparation/RSS_bytemaps_to_netcdf.ipynb) to unify the data and save in NetCDF format.
- **Note:** This script is prepared to convert data for **a full month**, so make sure data is available for a specific month, and use the script iteratively if you need to convert more than one month of data. The output of the script will be **a separate file for each month**.


#### MODIS_pbl and AMSR_TB:
- These datasets were generated by Ryan Eastman @UW.


#### MODIS/Aqua Aerosol Cloud:
- This is not included in the current trajectory package, but the guide is provided for future references.
- Go to [LAADS DAAC website](https://ladsweb.modaps.eosdis.nasa.gov/missions-and-measurements/products/MYD08_D3#overview) and select ["Data Archive"](https://ladsweb.modaps.eosdis.nasa.gov/archive/allData/61/MYD08_D3/) to get the data for the whole world or select ["Search for Product Files"](https://ladsweb.modaps.eosdis.nasa.gov/search/order/2/MYD08_D3--61) to select for a specific region, time period, and a set of data files. It is available from 4 July 2002 to near present time. 


#### Note regarding the GOES data:

Currently, analysis of GOES is not part of this trajectory package. It was done separately in the past and might be included in the future. The original codes to do so are provided. It needs GOES pixel data and trajectory text files. Download GOES pixel data from [here](https://satcorps.larc.nasa.gov/prod/goes-west/visst-pixel-netcdf). If you need GOES gridded data, you can download them from [visst website](https://www-angler.larc.nasa.gov/prod/goes-west/visst-grid-netcdf).
Other useful links to access GOES data: Customized GOES netcdf files can be ordered from: [NCEI website](https://www.ncei.noaa.gov/access/search/dataset-search?pageNum=3) or [NCDC website](https://www.ncdc.noaa.gov/airs-web/search)

#
### Generating Lagrangian trajectories: 

Two methods of creating trajectories have been tested:

-  Lagrangian forward trajectories are created isobarically from an initial level of 950 hPa for 72 hours, using ECMWF ERA5 winds and trajectory code developed at UW. (Reference for this code: [Eastman and Wood, 2016](https://doi.org/10.1175/JAS-D-15-0193.1))
- Both forward and backward trajectories are based on GFS analysis data at 500 m and for isobaric vertical motion (which is constant pressure surface). You can use the [HYSPLIT website](https://www.ready.noaa.gov/HYSPLIT_traj.php) and click on Compute archive trajectories.

#
### Possible issues:

- If you get this error:

TypeError: Using a DataArray object to construct a variable is ambiguous, please extract the data using the .data property.


it is due to the newest version of xarray library. Downgrade this library to the 0.18 version.


- If you get an error about proj_lib in the Basemap library, add these lines at the top of /MERRA2/add_to_trajectory.py before loading Basemap:


``` conda_file_dir = conda.__file__
conda_dir = conda_file_dir.split('lib')[0]
proj_lib = os.path.join(os.path.join(conda_dir, 'share'), 'proj')
os.environ["PROJ_LIB"] = proj_lib
from mpl_toolkits.basemap import Basemap 

```
#
### The analyses are featured in:
Erfani, E., Wood, R., Blossey, P., Doherty, S., & Eastman, R. (2022). Using a phase space of environmental variables to drive an ensemble of cloud-resolving simulations of low marine clouds, American Geophysical Union (AGU) Fall Meeting, Chicago, IL, 12-16 Dec. 2022, (poster). http://dx.doi.org/10.13140/RG.2.2.25888.66562/1

Erfani, E., Blossey, P., Wood, R., Mohrmann, J., Doherty, S. J., Wyant, M., & O, K. (2022). Simulating aerosol lifecycle impacts on the subtropical stratocumulus-to-cumulus transition using large-eddy simulations. Journal of Geophysical Research: Atmospheres, 127, e2022JD037258. https://doi.org/10.1029/2022JD037258

#
### Disclaimer:
The software is provided "as is", without warranty of any kind. In no event shall the authors be liable for any claim, damages, or other liability.
