#	This routine reads version-7 RSS SSMIS uncertainty files 
#
# 	filename  name of data file with path in form satname_yyyymmddv7.unc.gz
#	where satname      = name of satellite (ssmis)
#	         yyyy      = year
#		  	   mm      = month
#		       dd      = day of month
#           missing = fill value used for missing data;
#                     if None, then fill with byte codes (251-255)

#	the output values correspond to:
#	time_unc	time of measurement in fractional hours GMT within the file
#	wspd_mf_unc	input-induced uncertainty for the 10m surface wind, in meters/second
#	vapor_unc	input-induced uncertainty for columnar or integrated water vapor in millimeters
#	cloud_unc	input-induced uncertainty for cloud liquid water in millimeters
#	rain_unc	input-induced uncertainty for rain rate in millimeters/hour
#       longitude	Grid Cell Center Longitude', LON = 0.25*x_grid_location - 0.125 degrees east
#       latitude	Grid Cell Center Latitude',  LAT = 0.25*y_grid_location - 90.125
#       land		Is this land?
#       ice			Is this ice?
#       nodata		Is there no data
#
#   if you need help with this python read code, contact http:\\www.remss.com\support


from bytemaps import sys
from bytemaps import Dataset
from bytemaps import Verify


class SSMISuncertainty(Dataset):
    """ Read daily SSMIS uncertainty bytemaps. """
    """
    Public data:
        filename = name of uncertainty data file
        missing  = fill value used for missing data;
                  if None, then fill with byte codes (251-255)
        dimensions = dictionary of dimensions for each coordinate
        variables  = dictionary of data for each variable
    """

    def __init__(self, filename, missing=None):
        """
        Required arguments:
            filename = name of uncertainty data file to be read (string)
                
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
        return ('orbit_segment','variable','latitude','longitude')

    def _shape(self):
        return (2,5,720,1440)

    def _variables(self):
        return ['time_unc','wspd_mf_unc','vapor_unc','cloud_unc','rain_unc',
                'longitude','latitude','land','ice','nodata']                

    # _default_get():
    
    def _get_index(self,var):
        return {'time_unc' : 0,
                'wspd_mf_unc' : 1,
                'vapor_unc' : 2,
                'cloud_unc' : 3,
                'rain_unc' : 4,
                }[var]

    def _get_scale(self,var):
        return {'time_unc' : 0.1,
                'wspd_mf_unc' : 0.01,
                'vapor_unc' : 0.01,
                'cloud_unc' : 0.001,
                'rain_unc' : 0.002,
                }[var]


    # _get_ attributes:
    
    def _get_long_name(self,var):
        return {'time_unc' : 'Fractional Hour GMT for uncertainty file',
                'wspd_mf_unc' : 'Input-induced Uncertainty for 10m Surface Wind Speed',
                'vapor_unc' : 'Input-induced Uncertainty for Columnar Water Vapor',
                'cloud_unc' : 'Input-induced Uncertainty for Cloud Liquid Water',
                'rain_unc' : 'Input-induced Uncertainty for Surface Rain Rate',
                'longitude' : 'Grid Cell Center Longitude',
                'latitude' : 'Grid Cell Center Latitude',                
                'land' : 'Is this land?',
                'ice' : 'Is this ice?',
                'nodata' : 'Is there no data?',
                }[var]

    def _get_units(self,var):
        return {'time_unc' : 'Fractional Hour GMT',
                'wspd_mf_unc' : 'm/s',
                'vapor_unc' : 'mm',
                'cloud_unc' : 'mm',
                'rain_unc' : 'mm/hr',
                'longitude' : 'degrees east',
                'latitude' : 'degrees north',
                'land' : 'True or False',
                'ice' : 'True or False',
                'nodata' : 'True or False',
                }[var]

    def _get_valid_min(self,var):
        return {'time_unc' : 0.0,
                'wspd_mf_unc' : 0.0,
                'vapor_unc' : 0.0,
                'cloud_unc' : 0.0,
                'rain_unc' : 0.0,
                'longitude' : 0.0,
                'latitude' : -90.0,
                'land' : False,
                'ice' : False,
                'nodata' : False,
                }[var]

    def _get_valid_max(self,var):
        return {'time_unc' : 24.0,
                'wspd_mf_unc' : 2.5,
                'vapor_unc' : 2.5,
                'cloud_unc' : 0.25,
                'rain_unc' : 0.5,
                'longitude' : 360.0,
                'latitude' : 90.0,
                'land' : True,
                'ice' : True,
                'nodata' : True,
                }[var]
    

class UncVerify(Verify):
    """ Contains info for verification. """
    
    def __init__(self,dataset):
        self.filename = 'Verify_F17_SSMIS_V7_uncertainty.txt'
        self.ilon1 = 170
        self.ilon2 = 175
        self.ilat1 = 274
        self.ilat2 = 278        
        self.iasc = 1
        self.variables = ['time_unc','wspd_mf_unc','vapor_unc','cloud_unc','rain_unc']
        self.startline = {'time_unc' : 61,
                          'wspd_mf_unc' : 76,
                          'vapor_unc' : 92,
                          'cloud_unc' : 107,
                          'rain_unc' : 124 }
        Verify.__init__(self,dataset)
        

if __name__ == '__main__':
    """ Automated testing. """



    # read uncertainties:    
    ssmi = SSMISuncertainty('f17_20140221v7.unc.gz')
    if not ssmi.variables: sys.exit('problem reading file')
    
    # verify daily:
    verify = UncVerify(ssmi)
    if verify.success: print 'successful verification for uncertainty file'
    else: sys.exit('verification failed for daily')
    print

    print 'all tests completed successfully'
    
