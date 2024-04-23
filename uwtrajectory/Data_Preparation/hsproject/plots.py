# -*- coding: utf-8 -*-
"""
plots.py: plotting routines for CSET GOES data
"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import dates
import netCDF4 as nc4
from read_tdump import read_tdump
import os
import pandas as pd
from util import get_trajectory_data, datetime64_to_datetime,\
                 create_cloud_droplet_concentration_array, closest_val,\
                 TrajectoryError
print "pandas version: {}".format(pd.__version__)


def pcolormesh_cloud_particle_size():
    """
    Plots the cloud_particle_size variable from G15V03.0.NH.2015182.1600.PX.08K
    using pcolormesh and an appropriate colormap. Includes a title, colorbar,
    and axis labels.
    """
    dataset = nc4.Dataset('data/G15V03.0.NH.2015182.1600.PX.08K.NC')
    cps = dataset.variables['cloud_particle_size']
    lat = dataset.variables['latitude']
    lon = dataset.variables['longitude']
    # allows computer to access the proper dataset and variables for colormap
    plt.figure()
    # tells computer to make figure
    plt.pcolor(lon[:], lat[:], cps[:], cmap='RdBu')
    # tells the computer to create the map in pcolor and notes which
    # variables are displayed on which axes
    plt.title('Cloud Particle Size: 7/1/15 16:00Z')
    # makes title
    plt.axis([-180, -90, 0, 60])
    # constrains the area shown in the colormap
    plt.colorbar()
    # makes a scale which equates numbers to colors. now it's ready to show.
    plt.show()
    # shows colormap


def pcolormesh_g15_image(
        filename, varname, cmap=None, latlim=None, lonlim=None):
    """
    Assumes filename is a file containing goes 15 pixel image data, formatted
    like G15V03.0.NH.2015182.1600.PX.08K.NC. Takes in this filename and a
    variable name from the file that is assumed to have (image_y, image_x)
    coordinates. Optionally takes in a colormap string and latitude and
    longitude limits.
    Plots the variable name given from the filename given using the indicated
    colormap, if any (default otherwise), and for the given region, if
    indicated (all data otherwise).
    """
    from mpl_toolkits.basemap import Basemap
    # imports Basemap, letting the computer to overlay a colormap on a basemap
    dataset = nc4.Dataset(filename)
    var = dataset.variables[varname]
    lat = dataset.variables['latitude']
    lon = dataset.variables['longitude']
    # allows computer to access the proper dataset and variables for colormap
    plt.figure()
    if lonlim is None:
        lonmin = lon[:].min()
    else:
        lonmin = lonlim[0]
    if latlim is None:
        latmin = lat[:].min()
    else:
        latmin = latlim[0]
    if lonlim is None:
        lonmax = lon[:].max()
    else:
        lonmax = lonlim[1]
    if latlim is None:
        latmax = lat[:].max()
    else:
        latmax = latlim[1]
    m = Basemap(llcrnrlon=lonmin, llcrnrlat=latmin,
                urcrnrlon=lonmax, urcrnrlat=latmax, projection='merc')
    # this section of code determines whether the user of the code has changed
    # the longitude and latitude limits of the colormap. if not, the code
    # automatically sets these limits at the extreme values of the data. sets
    # projection to mercator
    m.drawmapboundary(fill_color='white')
    m.drawcoastlines(color='aqua')
    m.drawparallels(np.arange(-90, 91, 10), color='aqua')
    m.drawmeridians(np.arange(-180, 180, 15), color='aqua')
    # draws coastlines, and latitude/longitude lines on the map
    if cmap is None:
        cmap = 'RdBu'
    else:
        cmap = cmap
    # determines whether the user has changed the colormap type and if not,
    # sets it as red-blue.
    m.pcolor(lon[:], lat[:], var[:], cmap=cmap, latlon=True)
    # tells the computer to create the map in pcolor and notes which
    # variables are displayed on which axes
    plt.colorbar()
    # makes a scale which equates numbers to colors. now it's ready to show.
    plt.show()
    # shows colormap


def basemap():
    from mpl_toolkits.basemap import Basemap
    # imports Basemap, letting the computer to overlay a colormap on a basemap
    plt.figure()
    # tells computer to make figure
    lonmin = -180
    latmin = 0
    lonmax = -90
    latmax = 60
    # sets the longitude and latitude limits of the colormap
    m = Basemap(llcrnrlon=lonmin, llcrnrlat=latmin,
                urcrnrlon=lonmax, urcrnrlat=latmax,
                projection='merc')
    # sets the edges of the colormap at the long and lat limits set above, sets
    # projection to mercator
    m.drawmapboundary(fill_color='white')
    m.drawcoastlines(color='aqua', linewidth=1)
    m.drawparallels(np.arange(-90, 91, 10), color='aqua')
    m.drawmeridians(np.arange(-180, 180, 15), color='aqua')
    # draws coastlines and latitude/longitude lines onto the colormap in the
    # specified area
    dataset = nc4.Dataset('data/G15V03.0.NH.2015182.1600.PX.08K.NC')
    cps = dataset.variables['cloud_particle_size']
    lat = dataset.variables['latitude']
    lon = dataset.variables['longitude']
    # allows computer to access the proper dataset and variables for colormap
    m.pcolor(lon[:], lat[:], cps[:], cmap='hot', latlon=True)
    # tells the computer to create the map in pcolor and notes which
    # variables are displayed on which axes
    plt.title('Cloud Particle Size: 7/1/15 16:00Z')
    # makes title
    plt.colorbar()
    # makes a scale which equates numbers to colors. now it's ready to show.
    plt.show()
    # shows colormap


def is_night(lat, lon, dataset):
    i_lat = closest_val(lat, dataset.variables['grid_lat'][:])
    i_lon = closest_val(lon, dataset.variables['grid_lon'][:])
    return dataset.variables['solar_zenith_angle'][i_lat, i_lon] > 70.


def get_goes_data(lat_traj, lon_traj, date_traj, goes_foldername, vis_varnames,
                  ir_varnames, delta_index=20, beginning_time=0, end_time=23):

    # a list that stores numpy arrays of float data that comes from each netcdf
    # file
    out_of_bounds_flag = 0
    out_data = {varname: [] for varname in (vis_varnames + ir_varnames)}
    long_names = {}
    units = {}

    # cycles through all of the times since the times is the same length as the
    # lat_traj and lon_traj, to take data from each of the trajectory points
    # an add that to a data list and add all of the tiems to a list of times
    for x, time_item in enumerate(date_traj):
        # finds the day in the year which is necesary for the name of the
        # netcdf file that must be inputted
        day_in_year = time_item.timetuple().tm_yday
        # the file name for the netcdf file to get data from
        out_filename = 'G15V03.0.NH.{}{:03d}.{:02d}00.PX.08K.NC'.format(
                       time_item.year, day_in_year, time_item.hour)
        # try except statement incase the file does not exist
        try:

            # retreives data from the netcdf file
            dataset = nc4.Dataset(goes_foldername + out_filename)

            # adds the time for the file to the time_list
            night = is_night(lat_traj[x], lon_traj[x], dataset)
            # retreives the lat and lon coordinates that correspond to
            # points in the variable array
            lat = dataset.variables['latitude'][:]
            lon = dataset.variables['longitude'][:]

            # adds the data from the array in the specified lat_traj and
            # lon_traj into the data_list
            if night:
                varnames = ir_varnames
            else:
                varnames = vis_varnames + ir_varnames
            for varname in varnames:
                try:
                    append_data = get_trajectory_data(
                            dataset.variables[varname], lat, lon,
                            lat_traj[x:x+1], lon_traj[x:x+1],
                            delta_index).flatten()
                except TrajectoryError:
                        out_of_bounds_flag = 1
                        pass
                append_data[append_data ==
                            dataset.variables[varname]._FillValue] = np.NaN
                out_data[varname].append(append_data)
                if varname not in long_names:
                    long_names[varname] = dataset.variables[varname].long_name
                if varname not in units:
                    try:
                        units[varname] = dataset.variables[varname].units
                    except AttributeError:
                        units[varname] = 'unitless'
            if night:
                for varname in vis_varnames:
                    out_data[varname].append(
                        np.empty(((delta_index*2+1)**2,))*np.NaN)
            dataset.close()
            # in case the file does not exist
        except RuntimeError:
            print 'Caught RuntimeError'
            # adds the time so that a time is still represented in the
            # graph even though it should display as empty
            # adds an array of non-numbers to the data_list to show that
            # there was no data inputted at this time
            for varname in (vis_varnames + ir_varnames):
                out_data[varname].append(
                    np.empty(((delta_index*2+1)**2,))*np.NaN)

    if len(units.keys()) == 0:
        # no GOES files were found
        raise TrajectoryError('no GOES files found')
    # takes all of the data from the data_list and puts it into a numpy array
    # for ease of using numpy methods
    for varname in (vis_varnames + ir_varnames):
        out_data[varname] = np.asarray(out_data[varname])
    out_data['time'] = np.asarray(date_traj)
    print out_data['time']
    out_data['Nd'] = (
        1.4067 * 10**(4) * (out_data['cloud_visible_optical_depth'] ** 0.5) /
        (out_data['cloud_particle_size'] ** 2.5))
    long_names['Nd'] = 'cloud droplet number concentration'
    units['Nd'] = '1/cm^3'
    out_data['lat_traj'] = lat_traj
    out_data['lon_traj'] = lon_traj
    for k in units.keys():
        units[k] = str(units[k])
        long_names[k] = str(long_names[k])

    return out_data, long_names, units, out_of_bounds_flag


def get_goes_from_multitrajectory(trajectory_filename, goes_foldername,
                                  vis_varnames, ir_varnames, delta_index=20,
                                  beginning_time=0, end_time=23, lt=None,
                                  keylist=None):
    # retreives the information from the trajectory_filename to get all of the
    # trajectories necessary for inputting data from the netcdf files
    T = read_tdump(os.path.join(goes_foldername, trajectory_filename))
    t = T.groupby('tnum')
    out_data_list = []
    long_names_list = []
    units_list = []
    for tnum in t.groups.keys():
        if keylist is not None:
            if tnum not in keylist:
                print('skipping trajectory {}'.format(tnum))
                continue
        print('getting trajectory {}'.format(tnum))
        if lt is not None:
            lt.update(overwrite=False)
        traj = t.get_group(tnum)

        # inputs the lat_traj and lon traj from the trajectory_file for the lat and
        # lon coordinates that data would be taken from in each netcdf file which
        # each represent a time
        lat_traj = np.array(traj['lat'])
        lon_traj = np.array(traj['lon'])
        # retrieves all of the times to sample data
        date_traj = np.array([datetime64_to_datetime(item) for item in
                             np.array(traj['lat'].index)])
        try:
            out_data, long_names, units, oob_flag = get_goes_data(
                lat_traj, lon_traj, date_traj, goes_foldername,
                vis_varnames, ir_varnames, delta_index, beginning_time, end_time)
            if oob_flag:
                print 'Trajectory out of bounds: {}, tnum {}'.format(
                    trajectory_filename, tnum)
        except TrajectoryError as e:
            print '{}: tnum {}: {}'.format(trajectory_filename, tnum, e)
            continue
        out_data_list.append(out_data)
        long_names_list.append(long_names)
        units_list.append(units)

    if len(out_data_list) == 0:
        raise RuntimeError('No GOES data found for trajectory file')
#    if len(out_data_list) == 1:
#        return out_data_list[0], long_names_list[0], units_list[0]
    else:
        return out_data_list, long_names_list, units_list


def get_goes_from_flightfile(flight_filename, goes_foldername,
                             vis_varnames, ir_varnames, delta_index=20,
                             beginning_time=0, end_time=23, lt=None):
    with nc4.Dataset(flight_filename) as dataset:
        lats = dataset.variables['GGLAT'][:]
        lons = dataset.variables['GGLON'][:]
        date_var = dataset.variables['Time']
        dates = nc4.num2date(date_var[:], units=date_var.units)

        ds = pd.DataFrame(index=dates, data={'lat': lats, 'lon': lons})
        ds_lowres = ds.resample('15Min').mean()
        lat_traj = ds_lowres.lat.values
        lon_traj = ds_lowres.lon.values
        date_traj = np.array([datetime64_to_datetime(item) for item in
                             np.array(ds_lowres.index)])

        out_data, long_names, units, oob_flag = get_goes_data(
            lat_traj, lon_traj, date_traj, goes_foldername,
            vis_varnames, ir_varnames, delta_index, beginning_time, end_time)
        if oob_flag:
            print 'Trajectory out of bounds: {}, tnum {}'.format(
                    flight_filename)
            raise IndexError('trajectory out of bounds')
    return out_data, long_names, units


def goes_histogram(var_name, trajectory_filename, goes_foldername,
                   bin_spacing=None, max_height=None, vmin=None,
                   vmax=None, title=None, reduction_function=np.nanmedian,
                   delta_index=20, beginning_time=0, end_time=23):
    """
    Plot a 2D histogram of cloud top height in a small region around a
    trajectory specified in trajectory_filename (a HYSPLIT tdump file), using
    data from GOES pixel netCDF files contained in the folder specified by
    goes_foldername.

    Parameters
    ----------
    var_name: string
        The name in the netCDF file of the variable being plotted
    bin_spacing: int, optional
        The width of each bin on the y (data) axis, in the same units as the
        data. If no bin spacing is given, the computer automatically creates 20
        equally spaced bins on the histogram
    max_height: float, optional
        The height of the top of the histogram, in the same units as the data.
        If no height is given, the max height is set at the value of the
        highest data point
    vmin: float, optional
        The lowest data value used in calculating the mean or median of the
        data, in the same units as the data. If no vmin is given, there is
        assumed to be no lower limit to valid data values
    vmax: float, optional
        The highest data value used in calculating the mean or median of the
        data, in the same units as the data. If no vmax is given, there is
        assumed to be no upper limit to valid data values
    title: string, optional
        The title of the plot, defaults to var_name
    reduction_function: func, optional
        The function applied to the data to reduce it to 1-D time series,
        defaults to median.
    delta_index: int, optional
        The number of pixels out in each direction from the trajectory point
        that are used in the histogram, defaults to 20

    Returns
    -------
    time: ndarray
        The 1-D array of times (datetime objects) covered by the data and
        histogram
    data: ndarray
        The 2-D array of data used in the histogram

    """

    # for traj_num in range(0,23):
    tnum = 1
    # retreives the information from the trajectory_filename to get all of the
    # trajectories necessary for inputting data from the netcdf files
    tdump = read_tdump(os.path.join(goes_foldername, trajectory_filename))
    # inputs the lat_traj and lon traj from the trajectory_file for the lat and
    # lon coordinates that data would be taken from in each netcdf file which
    # each represent a time
    lat_traj = np.array(tdump[tdump['tnum'] == tnum]['lat'])
    lon_traj = np.array(tdump[tdump['tnum'] == tnum]['lon'])
    # retrieves all of the times to sample data
    time = np.array([datetime64_to_datetime(item)
                    for item in
                    np.array(tdump[tdump['tnum'] == tnum]['lat'].index)])
    # a list that stores all of the times taken from the trajectory_file
    time_list = []
    # a list that stores numpy arrays of float data that comes from each netcdf
    # file
    data_list = []
    # cycles through all of the times since the times is the same length as the
    # lat_traj and lon_traj, to take data from each of the trajectory points
    # an add that to a data list and add all of the tiems to a list of times
    for x in range(0, len(lat_traj)):
        # gets a datetime object that is retreived earlier from the trajectory
        # data
        time_item = time[x]
        print(time_item.hour)
        if((beginning_time <= end_time and time_item.hour <= end_time and
            time_item.hour >= beginning_time) or (beginning_time > end_time and
            (time_item.hour <= end_time or time_item.hour >= beginning_time))):
            # finds the day in the year which is necesary for the name of the
            # netcdf file that must be inputted
            day_in_year = time_item.timetuple().tm_yday
            # the file name for the netcdf file to get data from
            out_filename = 'G15V03.0.NH.{}{:03d}.{:02d}{:02d}.PX.08K.NC'.format(
                           time_item.year, day_in_year, time_item.hour,
                           time_item.minute)
            print(out_filename)
            # try except statement incase the file does not exist
            try:
                # retreives data from the netcdf file
                dataset = nc4.Dataset(goes_foldername + out_filename)
                # adds the time for the file to the time_list
                time_list.append(time_item)
                # retreives the lat and lon coordinates that correspond to
                # points in the variable array
                lat = dataset.variables['latitude']
                lon = dataset.variables['longitude']
                # adds the data from the array in the specified lat_traj and
                # lon_traj into the data_list
                if var_name == 'cloud_droplet_concentration':
                    array = (create_cloud_droplet_concentration_array
                                     ((dataset.variables['cloud_particle_size']),
                                      (dataset.variables
                                      ['cloud_visible_optical_depth']), lat[:],
                                      lon[:], lat_traj[x:x+1], lon_traj[x:x+1],
                                      delta_index).flatten())
                    data_list.append(array)
                else:
                    # gets the necesary information specified by the var_name
                    # for the type of data that the user wants to access
                    array = dataset.variables[var_name][:].filled(fill_value=np.nan)
                    print('printing array')
                    print(array)
                    print('ending print')
                    data_list.append(get_trajectory_data(array[:], lat[:],
                                                         lon[:],
                                                         lat_traj[x:x+1],
                                                         lon_traj[x:x+1],
                                                         delta_index).flatten()
                                     )
                dataset.close()
                # in case the file does not exist
            except RuntimeError:
                # adds the time so that a time is still represented in the
                # graph even though it should display as empty
                time_list.append(time_item)
                # adds an array of non-numbers to the data_list to show that
                # there was no data inputted at this time
                data_list.append(np.empty([(delta_index*2+1)**2])*np.NaN)
        else:
            # adds the time so that a time is still represented in the
            # graph even though it should display as empty
            time_list.append(time_item)
            # adds an array of non-numbers to the data_list to show that
            # there was no data inputted at this time
            data_list.append(np.empty([(delta_index*2+1)**2])*np.NaN)
    # takes all of the data from the data_list and puts it into a numpy array
    # for ease of using numpy methods
    data = np.asarray(data_list)
    # takes all of the dates in the time_list and makes them numbers
    times = dates.date2num(time_list)
    # takes the times and makes it into a numpy array
    times = np.asarray(times)
    # gets max and min of the times for setting graph width and for managing
    # the spacing of the ticks
    max_time = times.max()
    min_time = times.min()
    # takes all of the times and addes another time to the end of the array
    # which is 1/24 of a day(1 hour) after the last time so that the bin edges
    # can be inbetween each time
    x_bin_edges = np.append(times[:], (times[len(times) - 1] + 1. / 24))
    # subtracts 1/48 of a day so that the bin edges are in between every time
    # value which takes slot every hour
    x_bin_edges = x_bin_edges - 1 / 48.
    # the format for the time when displayed by the ticks on the graph
    hfmt = dates.DateFormatter('%m/%d %H:%M')
    # creates a figure with the pyplot to help with displaying ticks
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    # the number of time tick labels - 1 so this currently displays 11 ticks
    # and a corresponding amount of time labels.
    number_of_time_labels = 10
    # sets the tick marks by making the lowest time be the first tick and makes
    # the last tick time equal to max_time by making the point where tick marks
    # are no longer added at max_time + the distance between each tick which is
    # set at (max_time - min_time) / number_of_time_labels
    ax.xaxis.set_ticks(np.arange(min_time,
                       (max_time + (max_time - min_time) /
                        number_of_time_labels),
                       (max_time - min_time) / number_of_time_labels))
    # formats the time in each of the ticks with the desired format set above
    ax.xaxis.set_major_formatter(hfmt)
    # rotates the ticks to be verticle so that the tick labels are easier to
    # read
    plt.xticks(rotation='vertical')
    # creates an empty data array(temp_data_array) which will eventually be
    # filled by the subset of the data within the limits set by vmin and vmax
    temp_data_array = np.zeros_like(data) * np.NaN
    # finds which data points are between vmin and vmax and puts them in the
    # array valid_indices. If vmin is none, there is no lower bound on the data
    # and likewise if vmax is none. If both vmin and vmax are none, all the
    # data goes in valid_indices
    valid_indices = np.zeros(temp_data_array.shape, dtype=np.bool)
    if vmin is None and vmax is None:
        valid_indices[:] = True
    elif vmin is None:
        valid_indices[data < vmax] = True
    elif vmax is None:
        valid_indices[data > vmin] = True
    else:
        valid_indices[(data > vmin) & (data < vmax)] = True
        # fills temp_data_array with the data in valid_indices
    temp_data_array[valid_indices] = data[valid_indices]
    # runs a function chosen in argument reduction_function on each time
    # and then stores these values in stats_array
    stat_array = reduction_function(temp_data_array, axis=1)
    # tests what fraction of the data points within a certain delta index of
    # the trajectory points are within valid_indices for each time and if that
    # fraction is under 3/8, the output of the reduction function is nan
    num_hits = np.sum(valid_indices, axis=1)
    stat_array[num_hits < (2*delta_index+1)**2 * 3./8] = np.NaN
    # creates a copy of times array
    times1d = times
    # puts together the times and data arrays so that each data point
    # corresponds with the right time
    times_out, data_out = np.broadcast_arrays(times[:, None], data)
    # takes out every value in times and data that is nan
    times = times_out[~np.isnan(data)]
    data = data_out[~np.isnan(data)]
    # automatically sets the y value of the top of the graph equal to the
    # maximum value in the data if no max_height argument has been entered
    if max_height is None:
        max_height = data.max()
    # automatically creates 20 equally spaced bins in the y direction if no
    # bin_spacing argument has been entered
    if bin_spacing is None:
        bin_spacing = max_height/20
    # creates the histogram with the time on the x-axis, the data value on the
    # y-axis, and the number of data points in each bin using colors (colormap
    # type is set to Blues so the colors are different shades of blue). The bin
    # edges and maximum height are set using the user-entered or automatically
    # generated x_bin_edges, bin_spacing, and max_height values
    plt.hist2d(times, data, bins=[x_bin_edges, (max_height / bin_spacing)],
               range=[[min_time, max_time], [0, max_height]], cmap='Blues')
    # puts a red x on the graph at the value of the reduction function for each
    # time (no x is drawn if the value in stats_array is nan)
    plt.plot(times1d, stat_array, 'rx')
    # creates a colorbar for the histogram
    plt.colorbar()
    # puts a title on the histogram (it automatically sets the title to the
    # variable name [minus underscores] if the user has not entered a title and
    # otherwise uses the title the user has entered)
    if title is None:
        plt.title(var_name.replace('_', ' '))
    else:
        plt.title(title)
    # prints out the data in stats_array and the corresponding times in the
    # console. If a value in stats_array is nan, the computer prints
    # "insufficient data"
    for i in range(0, len(stat_array)):
        if np.isnan(stat_array[i]):
            print(dates.num2date(times1d[i]), 'insufficient data')
        else:
            print(dates.num2date(times1d[i]), stat_array[i])
    # shows the histogram
    plt.show()
    return times_out, data_out

if __name__ == '__main__':
    # pcolormesh_cloud_particle_size()
    trajectory_filename = ('analysis.UW_HYSPLIT_GFS.201507071830.airmass_tra'
                           'jectories.txt')
    # goes_foldername = '/Users/Dylan/Desktop/data/'
    goes_foldername = '/Users/kylebretherton/cset/data/'
    arg_dict = {
        'liquid+ice water path': {
            'args': ['cloud_lwp_iwp', trajectory_filename, goes_foldername],
            'kwargs': {
                'vmax': 500,
                'max_height': 200,
                'bin_spacing': 5,
                'reduction_function': np.nanmean,
                'delta_index': 20
            }
        },
        'cloud top height': {
            'args': ['cloud_top_height', trajectory_filename, goes_foldername],
            'kwargs': {
                'vmax': 3,
                'max_height': 3,
                'bin_spacing': .1,
                'reduction_function': np.nanmedian,
                'delta_index': 20
            }
        },
        'cloud droplet concentration': {
            'args': ['cloud_droplet_concentration', trajectory_filename,
                     goes_foldername],
            'kwargs': {
                'vmax': 500,
                'max_height': 500,
                'reduction_function': np.nanmedian,
                'delta_index': 12,
                'beginning_time': 16,
                'end_time': 2
            }
        },
        'cloud visible optical depth': {
            'args': ['cloud_visible_optical_depth', trajectory_filename,
                     goes_foldername],
            'kwargs': {
                'vmax': 150,
                'vmin': 0.1,
                'max_height': 50,
                'reduction_function': np.nanmedian,
                'delta_index': 20,
                'beginning_time': 16,
                'end_time': 2
            }
        }
    }
    # plot_choice = 'cloud droplet concentration'
    goes_histogram('cloud_visible_optical_depth', trajectory_filename,
                  goes_foldername, beginning_time=16, end_time=2)
    # goes_histogram(*arg_dict[plot_choice]['args'],
    #                **arg_dict[plot_choice]['kwargs'])
