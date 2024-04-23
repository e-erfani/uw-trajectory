from ftplib import FTP
import sys
import os
import glob
import numpy as np
import re
import xarray as xr
import datetime as dt

def get_ftp_files():
    server_url = 'ftp.remss.com'
    ftp = FTP(server_url)
    ftp.login(user=r'jkcm@uw.edu')#, passwd=r'jkcm@uw.edu')
    for year in ['2014', '2016']:
        ftp.cwd(f'/ascat/metopa/bmaps_v02.1/y{year}/')
        local_basedir = f'/home/disk/eos9/jkcm/Data/ascat/rss/{year}'
        for i in ftp.nlst():
            ftp.cwd(f'/ascat/metopa/bmaps_v02.1/y{year}/')
            ftp.cwd(i)
            if not os.path.exists(os.path.join(local_basedir, i)):
                os.makedirs(os.path.join(local_basedir, i))
            filenames = ftp.nlst()
            for filename in filenames:
                local_filename = os.path.join(local_basedir, i, filename)
                file = open(local_filename, 'wb')
                ftp.retrbinary('RETR '+ filename, file.write)
    ftp.close()
    
    
def convert_files_to_netcdf(files, dates, data_reader):
        datasets = [data_reader(i) for i in files]
        try:
            assert [np.allclose(i.variables['latitude'], datasets[0].variables['latitude']) for i in datasets]
            assert [np.allclose(i.variables['longitude'], datasets[0].variables['longitude']) for i in datasets]
        except AssertionError as e:
            print('caught')
            print(len(files))
            for i in datasets:
                print(i.variables['latitude'])
            raise e
        ds = xr.Dataset()
        ds['time'] = dates
        ds['latitude'] = datasets[0].variables['latitude']
        ds['longitude'] = datasets[0].variables['longitude']
        ds['orbit_segment'] = np.arange(2)
        for k, v in datasets[0].variables.items():
            if len(v.shape) == 3:
                concat_var = np.array([i.variables[k] for i in datasets])
                concat_var[concat_var<v.valid_min] = np.nan
                concat_var[concat_var>v.valid_max] = np.nan
                name = k if not k=='time' else 'UTCtime'
                ds[name] = (('time', 'orbit_segment', 'latitude', 'longitude'), concat_var)
            attrs = dict(long_name = v.long_name, units=v.units, valid_min=v.valid_min, valid_max=v.valid_max)
            if type(v.valid_min) == bool:
                attrs['valid_min'] = int(attrs['valid_min'])
                attrs['valid_max'] = int(attrs['valid_max'])
            ds[name].attrs = attrs
        ds.attrs['creation date'] = str(dt.datetime.utcnow())
        return ds    
    
    
def convert_ascat_to_netcdf():
    data_reader = make_ascat_data_reader()
    basefolder = '/home/disk/eos9/jkcm/Data/ascat/rss/'
    for year in ['2014', '2015', '2016']:
        months = [f'm{i+1:02}' for i in range(12)]
        for month in months:
            save_name = f'/home/disk/eos9/jkcm/Data/ascat/rss/all/ascat_unified_{year}-{month[1:]}.nc'
            mdir = os.path.join(basefolder, year, month)
            files = [os.path.join(mdir, f) for f in os.listdir(mdir) if re.match(r'ascat_[0-9]{8}_v02.1.gz', f)]
            dates = [dt.datetime.strptime(os.path.basename(i)[6:14], '%Y%m%d') for i in files]
            ds = convert_files_to_netcdf(files, dates, data+reader)
            ds.attrs['comments'] = "netcdf created by jkcm@uw.edu, adapted from bytemaps from Remote Sensing Systems. " +\
                                "http://remss.com/missions/ascat/"
            comp = dict(zlib=True, complevel=2)
            ds.to_netcdf(save_name, engine='h5netcdf', encoding={var: comp for var in ds.data_vars})

            
def make_ascat_data_reader():
    from ascat_daily import ASCATDaily
    def read_data(filename):
        dataset = ASCATDaily(filename, missing=np.nan)
        if not dataset.variables: sys.exit('file not found')
        return dataset
    return read_data
    
    
def make_amsr_data_reader():
    from amsr_proc.amsr2_daily import AMSR2daily
    def read_data(filename):
        dataset = AMSR2daily(filename, missing=np.nan)
        if not dataset.variables: 
            print(filename)
            sys.exit('file not found')
        return dataset
    return read_data

def make_ssmi_data_reader():
    from ssmi_proc.ssmi_daily_v7 import SSMIdaily
    def read_data(filename):
        dataset = SSMIdaily(filename, missing=np.nan)
        if not dataset.variables: 
            print(filename)
            sys.exit('file not found')
        return dataset
    return read_data
    
def convert_amsr_to_netcdf():
    data_reader = make_amsr_data_reader()
    basefolder = '/home/disk/eos9/jkcm/Data/amsr/rss/'
    for year in [2014, 2015, 2016]:
        for month in [f'{i+1:02}' for i in range(12)]:
            save_name = os.path.join(basefolder, 'all', f'amsr_unified_{year}-{month}.nc')
            files = sorted(glob.glob(os.path.join(basefolder, f'f34_{year}{month}[0-9][0-9]v8.gz')))
            dates = [dt.datetime.strptime(os.path.basename(i)[4:12], '%Y%m%d') for i in files]
            ds = convert_files_to_netcdf(files, dates, data_reader)
            ds.attrs['comments'] = "netcdf created by jkcm@uw.edu, adapted from bytemaps from Remote Sensing Systems. " +\
                                "http://remss.com/missions/amsr/"
            comp = dict(zlib=True, complevel=2)
            ds.to_netcdf(save_name, engine='h5netcdf', encoding={var: comp for var in ds.data_vars})

def convert_ssmi_to_netcdf():
    data_reader = make_ssmi_data_reader()
    basefolder = '/home/disk/eos9/jkcm/Data/ssmi/new'
    for year in [2017]:
        for month in ['08']:
            for sat,ver in zip(['f15', 'f16', 'f17', 'f18'], ['v06', 'v07', 'v07', 'v08']):
                print(sat)
                save_name = os.path.join(basefolder, 'all', f'ssmi_unified_{sat}_{year}-{month}.nc')
                print(os.path.join(basefolder, sat))
                files = sorted(glob.glob(os.path.join(basefolder, sat, '*', f'y{year}', f'm{month}', f'[fF][0-9][0-9]_{year}{month}[0-9][0-9]v[0-9].gz')))
#                 print(files)
                dates = [dt.datetime.strptime(os.path.basename(i)[4:12], '%Y%m%d') for i in files]
#                 print(dates)
                ds = convert_files_to_netcdf(files, dates, data_reader)
                ds.attrs['comments'] = "netcdf created by jkcm@uw.edu, adapted from bytemaps from Remote Sensing Systems. " +\
                                    "http://remss.com/missions/ssmi/"
                comp = dict(zlib=True, complevel=2)
                ds.to_netcdf(save_name, engine='h5netcdf', encoding={var: comp for var in ds.data_vars})            
            
            
if __name__ == "__main__":
#     print("I exist!")
    convert_ssmi_to_netcdf()