# -*- coding: utf-8 -*-
"""
util.py: utility functions to use in scripts.
"""
import numpy as np
from scipy import interpolate
from datetime import datetime, timedelta
import pytz
import math
from numba import autojit


class TrajectoryError(Exception):
    pass


def find_nearest_index(value, iterable):
    """
    Takes in a value and an iterable (like a list or array).

    Returns the index of the iterable whose corresponding value is closest to
    the given value.

    Raises TypeError if iterable is not iterable, or if comparison (> and <) is
    not defined for value's type.
    """
    closest_value_index = 0
    closest_value = np.abs(value - iterable[0])
    for x in range(1, len(iterable)):
        if(np.abs(iterable[x] - value) < closest_value):
            closest_value_index = x
            closest_value = np.abs(value - iterable[x])
    return closest_value_index


@autojit
def find_index_limits(axis, limits):
    """
    Takes in a 1-D monotonically increasing array and a tuple
    (min_value, max_value).

    Returns (min_index, max_index) which are the lowest
    and highest indicies in the 1-D array that are within the (inclusive)
    limits provided.

    Raises ValueError if there are no values on the axis within the requested
    limits.
    """
    high = -1
    low = 0
    # iterates through the numbers in the array to determine
    # which indices are closest to the limits while still within them.
    # the variables "high" and "low" are used to keep track of these values.
    for number in axis:
        if number <= limits[1]:
            high = high+1
            # for every number below the upper limit, we add one to the
            # variable "high" so that once we've iterated through all the
            # values, "high" will equal the highest index below the upper
            # limit
        else:
            high = high
            break
        # for computational efficiency, this stops the computer from iterating
        # through values once it has gone past the upper limit.
        if number < limits[0]:
            low = low+1
            # this works the same way,adding one to "low" for each value below
            # the lower limit so that once we've iterated through all the
            # values, "low" will equal the loweset index number above the lower
            # limit
        else:
            low = low
    return (low, high)


def closest_index(lat_traj, lon_traj, lat, lon):
    dist = ((lat - lat_traj)**2 + (lon - lon_traj)**2)**(0.5)
    return np.unravel_index(np.argmin(dist), dist.shape)


def closest_val(x, L):
    '''
    FROM JCMGIBBON'S ATMOS PACKAGE
    Finds the index value in an iterable closest to a desired value.

    Parameters
    ----------
    x : object
        The desired value.
    L : iterable
        The iterable in which to search for the desired value.

    Returns
    -------
    index : int
        The index of the closest value to x in L.

    Notes
    -----
    Assumes x and the entries of L are of comparable types.

    Raises
    ------
    ValueError:
        if L is empty
    '''
    # Make sure the iterable is nonempty
    if len(L) == 0:
        raise ValueError('L must not be empty')
    if isinstance(L, np.ndarray):
        # use faster numpy methods if using a numpy array
        return (np.abs(L-x)).argmin()
    # for a general iterable (like a list) we need general Python
    # start by assuming the first item is closest
    min_index = 0
    min_diff = abs(L[0] - x)
    i = 1
    while i < len(L):
        # check if each other item is closer than our current closest
        diff = abs(L[i] - x)
        if diff < min_diff:
            # if it is, set it as the new closest
            min_index = i
            min_diff = diff
        i += 1
    return min_index


def get_region(array, lat, lon, latlim, lonlim):
    """
    Takes in a 2D [lat, lon] regularly gridded data array, and corresponding
    1-D latitude and longitude arrays. Optionally takes in a tuple
    (lat_min, lat_max) and one with (lon_min, lon_max). If no limits are given,
    they are assumed to be infinite.

    Returns subarray, sublat, sublon where each array is the subset of the
    given arrays within the indicated latitude and longitude limits.

    Raises TypeError if any of array, lat, or lon are not ndarrays.
    Raises ValueError if the shape of lat and lon are not appropriate for the
    shape of the data array.
    """
    if not isinstance(array, np.ndarray):
        raise ValueError('array must be an ndarray')
    if latlim is None:
        min_lat_index = 1
        max_lat_index = len(lat) - 1
    else:
        # find the correct min and max lat indices
        min_lat_index, max_lat_index = find_index_limits(lat, latlim)
    if lonlim is None:
        min_lon_index = 1
        max_lon_index = len(lon) - 1
    else:
        # find the correc min and max lon indices
        min_lon_index, max_lon_index = find_index_limits(lon, lonlim)
    # now that we have the min and max indices, return the subsets of our
    # arrays that are within those index limits
        subarray = array[min_lat_index:max_lat_index+1,
                         min_lon_index:max_lon_index+1]
        sublat = lat[min_lat_index:max_lat_index+1]
        sublon = lon[min_lon_index:max_lon_index+1]
    return subarray, sublat, sublon


def get_trajectory_data(array, lat, lon, lat_traj, lon_traj,
                        delta_index):
    """
    Takes in a 2D data array, corresponding 2-D latitude and longitude arrays,
    1-D trajectory latitude and longitude arrays, and an optional number of
    indices away from the trajectory to sample data (default 0).

    Returns a 2D array whose first axis corresponds to the indices in lat_traj
    and lon_traj, and second axis is an arbitrary data point index.
    return_array[i,:] should contain all the data from the input array in a
    latitude-longitude box around the nearest point to
    (lat_traj[i], lon_traj[i]). This box should include that point and all
    points +/- delta_index from that point (so the side length of the box is
    2*delta_index + 1 indices).

    Raises IndexError if delta_index is negative.
    Raises TypeError if array, lat, or lon are not ndarrays, or if delta_index
    is not an integer.
    Raises ValueError if data is requested that is outside the range of the
    data array.
    """
    # the width of the window of values obtained at each time
    box_width = delta_index * 2 + 1
    # since there should be an index for every value
    # of lat_traj/lon_traj, length should be length of lat_traj
    data_array = np.empty([len(lat_traj), box_width * box_width])
    for time in range(0, len(lat_traj)):
        # position in the middle of a box in the lat and lon in lowest index
        # for the index of the current center point of the square where data is
        # obtained from in the array. crnt = current
        crnt_lat_index, crnt_lon_index = closest_index(
            lat_traj[time], lon_traj[time], lat, lon)
        if min(crnt_lat_index, crnt_lon_index) < delta_index or \
                crnt_lat_index > array.shape[0]-delta_index or \
                crnt_lon_index > array.shape[1]-delta_index:
            # we are too close to the edge of the GOES data to retrieve our box
            raise TrajectoryError('trajectory out of GOES bounds')
        data_array[time][:] = (array[(crnt_lat_index - delta_index):
                               (crnt_lat_index + delta_index+1),
                               (crnt_lon_index - delta_index):
                               (crnt_lon_index + delta_index+1)]).flatten()
    return data_array


def trajectory_stats(array, lat, lon, lat_traj, lon_traj, delta_index=0):
    """
    Takes in a 2D data array, corresponding 1-D latitude and longitude arrays,
    1-D trajectory latitude and longitude arrays, and an optional number of
    indices away from the trajectory to sample data (default 0).

    Returns mean, std which are both 1-D arrays whose axis corresponds to the
    indices in lat_traj and lon_traj, and contain the mean and standard
    deviation of the data in a lat-lon box of side length (2*delta_index + 1)
    indices centered on the point nearest to (lat_traj[i], lon_traj[i]).

    Raises IndexError if delta_index is negative.
    Raises TypeError if array, lat, or lon are not ndarrays, or if delta_index
    is not an integer.
    Raises ValueError if data is requested that is outside the range of the
    data array.
    """
    mean = np.empty([len(lat_traj)])
    std = np.empty([len(lat_traj)])
    # creates variables which will eventually track the mean and std of the
    # dataset being used
    side_length = delta_index * 2 + 1
    # variable that states the side length of the box which contains the data
    for traj_index in range(0, len(lat_traj)):
        current_mean = 0
        # creates a variable which is used to solve for the mean of the data
        crnt_lat_index = find_nearest_index(lat_traj[traj_index], lat)
        crnt_lon_index = find_nearest_index(lon_traj[traj_index], lon)
        # finds the lat/lon coords of the point in the center of the box
        for lat_index in range(-delta_index, delta_index + 1):
            for lon_index in range(-delta_index, delta_index + 1):
                # basically dictates the points that the computer will look at
                # and tells the computer to cycle through these values
                current_mean += array[crnt_lat_index +
                                      lat_index][crnt_lon_index + lon_index]
        current_mean = current_mean / (side_length**2)
        # while cycling through the values, the computer finds the total of the
        # array values and divides them by the number of points in the array to
        # get the mean
        mean[traj_index] = current_mean
        # sets mean equal to the final value of current_mean (this should
        # equal the the mean of the data)
        current_deviation = 0
        # creates a variable which is used to solve for the std of the data

        for lat_index in range(-delta_index, delta_index + 1):
            for lon_index in range(-delta_index, delta_index + 1):
                array_value = array[crnt_lat_index +
                                    lat_index][crnt_lon_index + lon_index]
        # makes computer call all of the array values and cycle through them
                current_deviation += (array_value - current_mean)**2
        # while cycling through the values, the computer adds up all the
        # squared deviations of the array points from the mean (the first
        # part of the method to calculate std)
        current_deviation = (current_deviation / (side_length**2))**(1/2)
        # the computer then divides the value from the previous step by the
        # number of points in the array and then square-roots the result
        # (the remaining steps to find std)
        std[traj_index] = current_deviation
        # sets std equal to the final value of current_deviation (this should
        # equal the std of the data)
    return mean, std


def datetime64_to_datetime(dt64):
    if isinstance(dt64, np.datetime64):
        return datetime.utcfromtimestamp(
            (dt64 - np.datetime64('1970-01-01T00:00:00Z')) /
            np.timedelta64(1, 's')).replace(tzinfo=pytz.UTC)
    elif isinstance(dt64, datetime):
        return dt64
    else:
        raise ValueError('dt64 must be np.datetime64 object')


def interpolate_to_regular_grid(array, lat_in, lon_in, lat_out, lon_out):
    """
    Takes in a 2D data array, corresponding 2-D latitude and longitude arrays,
    and 1-D output latitude and longitude arrays. Assumes the regular grid
    defined by lat_out, lon_out is contained within the input data.

    Returns a 2D array containing data from the input array linearly
    interpolated to a regular lat, lon grid based on lat_out and lon_out.

    Raises TypeError if any of the arguments are not ndnarrays.
    Raises ValueError if any array shapes do not conform to this docstring.
    """
    lon_in.flatten()  # turned into 1-D array for interpolation
    lat_in.flatten()
    array.flatten()
    interp_func = interpolate.interp2d(lat_in, lon_in,
                                       array, kind='linear')
    return interp_func(lat_out, lon_out)


def create_cloud_droplet_concentration_array(particle_radius,
                                             cloud_optical_depth, lat, lon,
                                             lat_traj, lon_traj, delta_index):
    box_width = delta_index * 2 + 1
    data_array = np.empty([len(lat_traj), box_width * box_width])
    np.ma.set_fill_value(data_array, np.NaN)
    for time in range(0, len(lat_traj)):
        crnt_lat_index, crnt_lon_index = closest_index(
            lat_traj[time], lon_traj[time], lat, lon)
        data_array[time][:] = ((1.4067 * 10**(4) * (cloud_optical_depth[
            crnt_lat_index-delta_index:crnt_lat_index+delta_index+1,
            crnt_lon_index-delta_index:crnt_lon_index+delta_index+1]) ** 0.5) /
            (particle_radius[
                crnt_lat_index-delta_index:crnt_lat_index+delta_index+1,
                crnt_lon_index-delta_index:crnt_lon_index+delta_index+1]
             )**2.5).flatten()
    return data_array
