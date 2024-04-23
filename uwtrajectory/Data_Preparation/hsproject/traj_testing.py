# -*- coding: utf-8 -*-
"""
Created on Mon Jun 22 15:06:12 2015

@author: mcgibbon
"""
from goes import plot_goes as goesplot
from read_tdump import read_tdump
import numpy as np
from datetime import datetime
import pytz
import os


def datetime64_to_datetime(dt64):
    if isinstance(dt64, np.datetime64):
        return datetime.utcfromtimestamp(
            (dt64 - np.datetime64('1970-01-01T00:00:00Z')) /
            np.timedelta64(1, 's')).replace(tzinfo=pytz.UTC)
    elif isinstance(dt64, datetime):
        return dt64
    else:
        raise ValueError('dt64 must be np.datetime64 object')


if __name__ == '__main__':
    goes_dir = '/home/disk/eos4/jkcm/Data/GOES/2015/'
    output_dir_base = 'traj_20150617H18_8degx12deg_'
    tdump_filename = 'tdump20150617H18'
    for tnum in range(1, 6):
        tdump = read_tdump(tdump_filename)
        lat = np.array(tdump[tdump['tnum'] == tnum]['lat'])
        lon = np.array(tdump[tdump['tnum'] == tnum]['lon'])
        time = np.array([datetime64_to_datetime(item)
                         for item in
                         np.array(tdump[tdump['tnum'] == tnum]['lat'].index)])
        output_dir = output_dir_base + str(tnum)
        try:
            os.mkdir(output_dir)
        except OSError:
            pass
        goesplot.animate_goes(
            goes_dir, output_dir, datetime_start=time[0],
            datetime_end=time[-1], follow_lon=lon,
            follow_lat=lat, follow_datetimes=time, lon_length=12.,
            lat_length=8., parallelsep=2.5, meridiansep=2.5, follow=True,
            step_vis=1, step_IR=1, figsize=(8., 6.),  mode='cset')
        goesplot.smooth_folder(output_dir, goesplot.cset_regex_str + '.png')
        os.system(
            "ffmpeg -framerate {framerate:d} -y -pattern_type glob "
            "-i '{output_dir}/*.png' -pix_fmt yuv420p -c:v libx264 "
            "{output_dir}.mp4".format(framerate=12, output_dir=output_dir))
