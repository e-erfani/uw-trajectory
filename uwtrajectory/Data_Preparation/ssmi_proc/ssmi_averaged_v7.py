from bytemaps import sys
from bytemaps import Dataset
from bytemaps import Verify


class SSMIaveraged(Dataset):
    """ Read averaged SSMI bytemaps. """
    """
    Public data:
        filename = name of data file
        missing = fill value used for missing data;
                  if None, then fill with byte codes (251-255)
        dimensions = dictionary of dimensions for each coordinate
        variables = dictionary of data for each variable
    """

    def __init__(self, filename, missing=None):
        """
        Required arguments:
            filename = name of data file to be read (string)
                
        Optional arguments:
            missing = fill value for missing data,
                      default is the value used in verify file
        """       
        self.filename = filename
        self.missing = missing
        Dataset.__init__(self)

    # Dataset:

    def _attributes(self):
        return ['coordinates','long_name','units','valid_min','valid_max']

    def _coordinates(self):
        return ('variable','latitude','longitude')

    def _shape(self):
        return (4,720,1440)

    def _variables(self):
        return ['wspd_mf','vapor','cloud','rain',
                'longitude','latitude','land','ice','nodata']                

    # _default_get():
    
    def _get_index(self,var):
        return {'wspd_mf' : 0,
                'vapor' : 1,
                'cloud' : 2,
                'rain' : 3,
                }[var]

    def _get_scale(self,var):
        return {'wspd_mf' : 0.2,
                'vapor' : 0.3,
                'cloud' : 0.01,
                'rain' : 0.1,
                }[var]

    def _get_offset(self,var):
        return {'cloud' : -0.05,
                }[var]

    # _get_ attributes:
    
    def _get_long_name(self,var):
        return {'wspd_mf' : '10m Surface Wind Speed',
                'vapor' : 'Columnar Water Vapor',
                'cloud' : 'Cloud Liquid Water',
                'rain' : 'Surface Rain Rate',
                'longitude' : 'Grid Cell Center Longitude',
                'latitude' : 'Grid Cell Center Latitude',
                'land' : 'Is this land?',
                'ice' : 'Is this ice?',
                'nodata' : 'Is there no data?',
                }[var]

    def _get_units(self,var):
        return {'wspd_mf' : 'm/s',
                'vapor' : 'mm',
                'cloud' : 'mm',
                'rain' : 'mm/hr',
                'longitude' : 'degrees east',
                'latitude' : 'degrees north',
                'land' : 'True or False',
                'ice' : 'True or False',
                'nodata' : 'True or False',
                'land' : 'True or False',
                'ice' : 'True or False',
                'nodata' : 'True or False',
                }[var]

    def _get_valid_min(self,var):
        return {'wspd_mf' : 0.0,
                'vapor' : 0.0,
                'cloud' : -0.05,
                'rain' : 0.0,
                'longitude' : 0.0,
                'latitude' : -90.0,
                'land' : False,
                'ice' : False,
                'nodata' : False,
                }[var]

    def _get_valid_max(self,var):
        return {'wspd_mf' : 50.0,
                'vapor' : 75.0,
                'cloud' : 2.45,
                'rain' : 25.0,
                'longitude' : 360.0,
                'latitude' : 90.0,
                'land' : True,
                'ice' : True,
                'nodata' : True,
                }[var]


class ThreedayVerify(Verify):
    """ Contains info for verification. """
    
    def __init__(self,dataset):
        self.filename = 'Verify_SSMI_v7.txt'
        self.ilon1 = 170
        self.ilon2 = 175
        self.ilat1 = 274
        self.ilat2 = 278        
        self.variables = ['wspd_mf','vapor','cloud','rain']
        self.startline = {'wspd_mf' : 102,
                          'vapor' : 109,
                          'cloud' : 116,
                          'rain' : 123 }
        Verify.__init__(self,dataset)
        

class WeeklyVerify(Verify):
    """ Contains info for verification. """
    
    def __init__(self,dataset):
        self.filename = 'Verify_SSMI_v7.txt'
        self.ilon1 = 170
        self.ilon2 = 175
        self.ilat1 = 274
        self.ilat2 = 278        
        self.variables = ['wspd_mf','vapor','cloud','rain']
        self.startline = {'wspd_mf' : 133,
                          'vapor' : 140,
                          'cloud' : 147,
                          'rain' : 154 }
        Verify.__init__(self,dataset)
        

class MonthlyVerify(Verify):
    """ Contains info for verification. """
    
    def __init__(self,dataset):
        self.filename = 'Verify_SSMI_v7.txt'
        self.ilon1 = 170
        self.ilon2 = 175
        self.ilat1 = 274
        self.ilat2 = 278        
        self.variables = ['wspd_mf','vapor','cloud','rain']
        self.startline = {'wspd_mf' : 165,
                          'vapor' : 172,
                          'cloud' : 179,
                          'rain' : 186 }
        Verify.__init__(self,dataset)
        

if __name__ == '__main__':
    """ Automated testing. """    

    # read 3-day averaged:
    ssmi = SSMIaveraged('f10_19950120v7_d3d.gz')
    if not ssmi.variables: sys.exit('file not found')

    # verify 3-day:
    verify = ThreedayVerify(ssmi)
    if verify.success: print('successful verification for 3-day')
    else: sys.exit('verification failed for 3-day')
    print('')

    # read weekly averaged:
    ssmi = SSMIaveraged('f10_19950121v7.gz')
    if not ssmi.variables: sys.exit('file not found')

    # verify weekly:
    verify = WeeklyVerify(ssmi)
    if verify.success: print('successful verification for weekly')
    else: sys.exit('verification failed for weekly')     
    print('')
    
    # read monthly averaged:
    ssmi = SSMIaveraged('f10_199501v7.gz')
    if not ssmi.variables: sys.exit('file not found')
    
    # verify:
    verify = MonthlyVerify(ssmi)
    if verify.success: print('successful verification for monthly')
    else: sys.exit('verification failed for monthly')      
    print('')
    
    print('all tests completed successfully')
    print('')
