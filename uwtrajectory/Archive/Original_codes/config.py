########################################
###  Created by Johannes Mohrmann    ###
### Modified by Ehsan Erfani         ###
########################################

from pathlib import Path

##### output path:
trajectory_netcdf_dir = r'/home/disk/eos3/erfani/Data/CSET/Lagrangian_project/trajectory_files/'
# Note: Within this directory, create a folder named "old". 
# The code will attempt to move the existing outputs to old folder to avoild overwriting. 

#### Region:
region='nep'

#### input data paths
data_dir = Path('/home/disk/eos4/jkcm/Data')

MODIS_pbl_dayfile   = data_dir / 'CSET' / 'Ryan/Daily_1x1_JHISTO_CTH_c6_day_v2_calboxes_top10_Interp_hif_zb_2011-2016_corrected.nc'
MODIS_pbl_nightfile = data_dir / 'CSET' / 'Ryan/Daily_1x1_JHISTO_CTH_c6_night_v2_calboxes_top10_Interp_hif_zb_2011-2016.nc'
AMSR_fmt            = '/home/disk/eos9/jkcm/Data/amsr/rss/all/amsr_unified_*-*.nc'
CERES_fmt           = '/home/disk/eos9/jkcm/Data/ceres/proc/split/CERES_SYN1deg-1H_Terra-Aqua-MODIS_Ed4.*.nc'  

GOES_file_fmt       = data_dir / 'CSET/GOES/VISST_pixel/G15V03.0.NH.{}{:03}.{}*.NC'

if region=='nep':
    MERRA_dir = data_dir / 'MERRA' / '3h' / 'more_vertical'
    MERRA_fmt = "MERRA2_400.inst3_3d_aer_Nv.{:%Y%m%d}.nc4.nc4"
    
    ERA_source     = data_dir / 'CSET' / 'ERA5'
    ERA_ens_source = ERA_source / 'ensemble'
    ERA_ens_fmt    = "ERA5.enda.pres.NEP.{:%Y-%m-%d}.nc"
    ERA_fmt        = "ERA5.pres.NEP.{:%Y-%m-%d}.nc"
    ERA_sfc_fmt    = "ERA5.sfc.NEP.{:%Y-%m-%d}.nc"
    
    SSMI_file_fmt  = '/home/disk/eos9/jkcm/Data/ssmi/all/ssmi_unified_{}*.nc'

    
elif region=='barbados':
    #todo fill in some variables here
    pass

elif region=='sea':
    MERRA_dir = data_dir / 'MERRA' / 'sea' / 'new'
    MERRA_fmt = "MERRA2_400.inst3_3d_aer_Nv.{:%Y%m%d}.SUB.nc"
    
    ERA_source    = data_dir / 'ORACLES/ERA5'
    ERA_fmt       = "ERA5.pres.SEA.{:%Y-%m-%d}.nc"
    ERA_sfc_fmt   = "ERA5.sfc.SEA.{:%Y-%m-%d}.nc"
    SSMI_file_fmt = '/home/disk/eos9/jkcm/Data/ssmi/new/all/ssmi_unified_{}*.nc'
    
#### HYSPLIT CONFIG
# working directory for HYSPLIT to write CONTROL file. need write access
HYSPLIT_working_dir = data_dir / 'HYSPLIT/working' 
# write directory for HYSPLIT output files. need write access (can be anywhere)
HYSPLIT_tdump_dir   = data_dir / 'HYSPLIT/tdump/cset_profiles'
# pathname for HYSPLIT executable. need execute access
HYSPLIT_call        = '/home/disk/p/jkcm/hysplit/trunk/exec/hyts_std'
# write directory for saving plots
plot_dir = '~'
# read directory for HYSPLIT data. This shouldn't need changing unless you're downloading analysis I don't have
HYSPLIT_source_dir  = data_dir / 'HYSPLIT/source'
# write directory for geostationary imagery
imagery_dir         = data_dir / 'Minnis'





# project_dir = r'/home/disk/eos4/jkcm/Data/CSET/Lagrangian_project'
# trajectory_dir = os.path.join(project_dir, 'Trajectories')
#trajectory_netcdf_dir = r'/home/disk/eos4/jkcm/Data/CSET/Lagrangian_project/trajectory_files/'
# GOES_source = '/home/disk/eos4/jkcm/Data/CSET/GOES/VISST_pixel'
# GOES_trajectories = '/home/disk/eos4/jkcm/Data/CSET/GOES/flight_trajectories/data'
# GOES_flights = '/home/disk/eos4/jkcm/Data/CSET/GOES/flightpath/GOES_netcdf'
# dropsonde_dir = '/home/disk/eos4/jkcm/Data/CSET/AVAPS/NETCDF'
# latlon_range = {'lat': (15, 50), 'lon': (-160, -110)}
# HYSPLIT_workdir = '/home/disk/eos4/jkcm/Data/HYSPLIT/working'  # storing CONTROL
# HYSPLIT_call = '/home/disk/p/jkcm/hysplit/trunk/exec/hyts_std'  # to run HYSPLIT
# HYSPLIT_source = '/home/disk/eos4/jkcm/Data/HYSPLIT/source'
# ERA_ens_source = r'/home/disk/eos4/jkcm/Data/CSET/ERA5/ensemble'
# ERA_ens_temp_source = r'/home/disk/eos4/jkcm/Data/CSET/ERA5/ens_temp'
# # MERRA_source = r'/home/disk/eos4/jkcm/Data/CSET/MERRA'
# MERRA_source = r'/home/disk/eos4/jkcm/Data/MERRA/3h'
# base_date = dt.datetime(2015, 7, 1, 0, 0, 0, tzinfo=pytz.UTC)
# CSET_flight_dir = r'/home/disk/eos4/jkcm/Data/CSET/flight_data'
# sausage_dir = '/home/disk/eos4/jkcm/Data/CSET/sausage'
# plot_dir = r'/home/disk/p/jkcm/plots/lagrangian_paper_figures'
# flight_trajs = '/home/disk/eos4/jkcm/Data/CSET/Trajectories'
