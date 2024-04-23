from amsr2_daily import AMSR2daily

def read_data(filename='f34_20121206v8.gz'):
    dataset = AMSR2daily(filename, missing=missing)
    if not dataset.variables: sys.exit('problem reading file')
    return dataset

ilon = (169,174)
ilat = (273,277)
iasc = 0
avar = 'vapor'
missing = -999.

#----------------------------------------------------------------------------

def show_dimensions(ds):
    print
    print 'Dimensions'
    for dim in ds.dimensions:
        print ' '*4, dim, ':', ds.dimensions[dim]

def show_variables(ds):
    print
    print 'Variables:'
    for var in ds.variables:
        print ' '*4, var, ':', ds.variables[var].long_name

def show_validrange(ds):
    print
    print 'Valid min and max and units:'
    for var in ds.variables:
        print ' '*4, var, ':', \
              ds.variables[var].valid_min, 'to', \
              ds.variables[var].valid_max,\
              '(',ds.variables[var].units,')'

def show_somedata(ds):
    print
    print 'Show some data for:',avar
    print 'index range: (' + str(iasc) + ', ' + \
          str(ilat[0]) + ':' + str(ilat[1]) + ' ,' + \
          str(ilon[0]) + ':' + str(ilon[1]) + ')'
    print ds.variables[avar][iasc, ilat[0]:ilat[1]+1, ilon[0]:ilon[1]+1]

#----------------------------------------------------------------------------

if __name__ == '__main__':
    import sys
    dataset = read_data()
    show_dimensions(dataset)
    show_variables(dataset)
    show_validrange(dataset)
    show_somedata(dataset)
    print
    print 'done'
