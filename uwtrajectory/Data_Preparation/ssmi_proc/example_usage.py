from ssmis_daily_v7 import SSMISdaily

def read_data(filename='f17_20090120v7.gz'):
    dataset = SSMISdaily(filename, missing=missing)
    if not dataset.variables: sys.exit('file not found')
    return dataset

ilon = (169,174)
ilat = (273,277)
iasc = 0
avar = 'vapor'
missing = -999.

#----------------------------------------------------------------------------

def show_dimensions(ds):
    print('')
    print('Dimensions')
    for dim in ds.dimensions:
        aline = ' '.join([' '*3, dim, ':', str(ds.dimensions[dim])])
        print(aline)

def show_variables(ds):
    print('')
    print('Variables:')
    for var in ds.variables:
        aline = ' '.join([' '*3, var, ':', ds.variables[var].long_name])
        print(aline)

def show_validrange(ds):
    print('')
    print('Valid min and max and units:')
    for var in ds.variables:
        aline = ' '.join([' '*3, var, ':',
                str(ds.variables[var].valid_min), 'to',
                str(ds.variables[var].valid_max),
                '(',ds.variables[var].units,')'])
        print(aline)

def show_somedata(ds):
    print('')
    print(' '.join(['Show some data for:',avar]))
    aline = ''.join(['index range: (', str(iasc), ', ',
            str(ilat[0]), ':', str(ilat[1]), ', ',
            str(ilon[0]), ':', str(ilon[1]), ')'])
    print(aline)
    for klat in [ilat[0]+i for i in range(ilat[1]-ilat[0]+1)]:
        aline = ' '.join([str(itm) for itm in ds.variables[avar][ \
                 iasc, klat, ilon[0]:ilon[1]+1 ] ])
        print(aline)

#----------------------------------------------------------------------------

if __name__ == '__main__':    
    import sys
    dataset = read_data()
    show_dimensions(dataset)
    show_variables(dataset)
    show_validrange(dataset)
    show_somedata(dataset)
    print('')
    print('done')
    print('')
