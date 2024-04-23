############################################################
### Created by Johannes Mohrmann                         ###
### Modified by Ehsan Erfani (2022):                     ###
### Changes in parameters and functions to facilitate    ###
### automation and to generalize the code for projects   ###
### other than CSET                                      ###
############################################################

from pathlib import Path

##### Parameters to change:
data_dir = Path('/home/disk/eos3/erfani/Data')
ryan_dir = '/home/disk/eos10/rmeast/efold/climatologies/'
camp     = 'NEP'  # Campaign name (CSET, Ship_tracks, NEP, ...)

traj_txt_dir = data_dir / camp /  'Trajectories/'  # Path to txt files containing trajectory information
traj_prefix  = 'analysis.ERA5'   # prefix of trajectory file name

trajectory_netcdf_dir = data_dir / camp / 'Lagrangian_project/'   ## Output path
# Note: Within this directory, create a folder named "old". 
# The code will attempt to move the existing outputs to old folder to avoild overwriting. 

region = 'nep'   ## Region. Currently: North East Pacific

#### input data paths
if region == 'nep':
    ERA_source     = data_dir / camp / 'ERA5'
    ERA_ens_source = ERA_source / 'ensemble'
    ERA_ens_fmt    = "ERA5.enda.pres.NEP.{:%Y-%m-%d}.nc"
    ERA_fmt        = "ERA5.pres.NEP.{:%Y-%m-%d}.nc"
    ERA_sfc_fmt    = "ERA5.sfc.NEP.{:%Y-%m-%d}.nc"

    MERRA_dir      = data_dir / camp / 'MERRA' / '3h' / 'more_vertical'
    MERRA_fmt      = "MERRA2_400.inst3_3d_aer_Nv.{:%Y%m%d}.nc4.nc4"
    ######################
    ## Modified by Ehsan Erfani to include file format for years after 2020:
    MERRA_fmt_20s  = "MERRA2_401.inst3_3d_aer_Nv.{:%Y%m%d}.nc4.nc4"
    ######################
    CERES_fmt          = str(data_dir / camp / 'ceres/proc/split/CERES_SYN1deg-1H_Terra-Aqua-MODIS_Ed4.*.nc')

    SSMI_file_fmt  = str(data_dir / camp / 'ssmi/all/ssmi_unified_{}*.nc')

    AMSR_fmt           = str(data_dir / camp / 'amsr/all/amsr_unified_*-*.nc')

    # AMSR_Tb_fmt        = [[This should be set in the AMSR_Tb code]]

    MODIS_pbl_dayfile  = ryan_dir + 'Daily_1x1_JHISTO_CTH_c6_day_v2_calboxes_top10_Interp_hif_zb_2012-2021.nc'
    MODIS_pbl_nightfile= ryan_dir + 'Daily_1x1_JHISTO_CTH_c6_night_v2_calboxes_top10_Interp_hif_zb_2012-2021.nc'
    
    GOES_file_fmt      = data_dir / camp / 'GOES/VISST_pixel/G15V03.0.NH.{}{:03}.{}*.NC'  ## This is not used for now
    
elif region == 'barbados':
    #todo: fill in some variables here
    pass

elif region == 'sea':
    #todo: fill in some variables here
    pass

#######################################
## Available satellites for SSMI:
sats = ['f16' ,'f17', 'f18']  # Shiptrack 
# sats = ['f15', 'f16', 'f17', 'f18']  # CSET, which was in 2015 and included f15

#######################################
## 2019 Ship Track case names and dates

# "ST01" is reference number  and is prefix of the output name.
# "1.0"  is trajectory number and is suffix of the output name.
all_cases = {
    1: {'TLC_name': 'ST01_1.0',       
        'trajectories': [0, 1]}   # Don't change this.
        }
# start_dates = prescribe_date  ## deprecated and moved to the main code

#############################################################################
#### No need to change below this. These are parts of the original code. ####
#############################################################################

# #cset flight case names
# all_cases = {
#     1: {'ALC_name': 'ALC_RF02B-RF03CD',
#         'TLC_name': 'TLC_RF02-RF03_1.0-1.5-2.0',            #opt 1.0, fine
#         'trajectories': [0, 1]},
#     2: {'ALC_name': 'ALC_RF02C-RF03AB',
#         'TLC_name': 'TLC_RF02-RF03_0.5-1.0',                #opt 1.0, fine
#         'trajectories': [0, 1]},
#     3: {'ALC_name': 'ALC_RF04A-RF05CDE',
#         'TLC_name': 'TLC_RF04-RF05_2.0-2.3-2.5-3.0',            #opt 2.0. check
#         'trajectories': [0, 1]},
#     4: {'ALC_name': 'ALC_RF04BC-RF05AB',
#         'TLC_name': 'TLC_RF04-RF05_1.0-2.0',                #opt 2.0, ok
#         'trajectories': [0, 1]},
#     5: {'ALC_name': 'ALC_RF06A-RF07BCDE',
#         'TLC_name': 'TLC_RF06-RF07_3.5-4.0-4.3-4.6-5.0',        #opt 3.0, check 3.5
#         'trajectories': [0, 1]},
#     6: {'ALC_name': 'ALC_RF06BC-RF07A',
#         'TLC_name': 'TLC_RF06-RF07_1.6-2.0-2.3-2.6-3.0',    #opt 1.6, check
#         'trajectories': [0, 1]},
#     7: {'ALC_name': 'ALC_RF08A-RF09DEF',
#         'TLC_name': 'TLC_RF08-RF09_4.0-4.5-5.0',
#         'trajectories': [0, 1]},
#     8: {'ALC_name': 'ALC_RF08B-RF09BC',
#         'TLC_name': 'TLC_RF08-RF09_3.0-3.5', 
#         'trajectories': [0, 1]},
#     9: {'ALC_name': 'ALC_RF08CD-RF09A',
#         'TLC_name': 'TLC_RF08-RF09_1.5-2.0', 
#         'trajectories': [0, 1]},
#     10: {'ALC_name': 'ALC_RF10A-RF11DE',
#         'TLC_name': 'TLC_RF10-RF11_5.5-6.0',                #opt 5.0, removed 
#         'trajectories': [0, 1]},
#     11: {'ALC_name': 'ALC_RF10BC-RF11BC',
#         'TLC_name': 'TLC_RF10-RF11_3.0-3.5-4.0-5.0',        #opt 5.0, fine
#         'trajectories': [0, 1]},
#     12: {'ALC_name': 'ALC_RF10D-RF11A',
#         'TLC_name': 'TLC_RF10-RF11_1.0-1.5',                #opt 1.0, ok
#         'trajectories': [0, 1]},
#     13: {'ALC_name': 'ALC_RF12A-RF13E',
#         'TLC_name': 'TLC_RF12-RF13_4.5',                    #opt 5.0, removed
#         'trajectories': [0, 1]},
#     14: {'ALC_name': 'ALC_RF12B-RF13CD',
#         'TLC_name': 'TLC_RF12-RF13_3.0-3.5',                #added 3.0, ok
#         'trajectories': [0, 1]},
#     15: {'ALC_name': 'ALC_RF12C-RF13B',
#         'TLC_name': 'TLC_RF12-RF13_2.5-3.0',                
#         'trajectories': [0, 1]},
#     16: {'ALC_name': 'ALC_RF14A-RF15CDE',
#         'TLC_name': 'TLC_RF14-RF15_3.5-4.0',            
#         'trajectories': [0, 1]},
#     17: {'ALC_name': 'ALC_RF14B-RF15B',
#         'TLC_name': 'TLC_RF14-RF15_3.0',
#         'trajectories': [0, 1]},    
#     18: {'ALC_name': 'ALC_RF14CD-RF15A',
#         'TLC_name': 'TLC_RF14-RF15_1.0-2.0', 
#         'trajectories': [0, 1]}
# }


#### HYSPLIT CONFIG
# working directory for HYSPLIT to write CONTROL file. need write access
HYSPLIT_working_dir = data_dir / camp / 'HYSPLIT/working' 
# write directory for HYSPLIT output files. need write access (can be anywhere)
HYSPLIT_tdump_dir   = data_dir / camp / 'HYSPLIT/tdump/profiles'
# pathname for HYSPLIT executable. need execute access
HYSPLIT_call        = data_dir / camp / 'HYSPLIT/trunk/exec/hyts_std'
# write directory for saving plots
plot_dir = '~'
# read directory for HYSPLIT data. This shouldn't need changing unless you're downloading analysis I don't have
HYSPLIT_source_dir  = data_dir / camp / 'HYSPLIT/source'
# write directory for geostationary imagery
imagery_dir         = data_dir / camp / 'Minnis'





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
