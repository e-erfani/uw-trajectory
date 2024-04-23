# -*- coding: utf-8 -*-
"""
Created on Thu Aug 18 11:46:43 2016
Contains all code for extracting GOES along a trajectory
@author: jkcm

Code mostly written by jmcgibbon, modified by jkcm
"""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import os
from hsproject.plots import get_goes_from_multitrajectory
from hsproject.util import TrajectoryError
from colormaps import viridis
import xray
import datetime as dt
import pytz
import utils
import sys
from LoopTimer import LoopTimer

pairs = (
    (6, '1deg'),
   (13, '2deg'),
   (27, '4deg'),
   (65, '10deg')
)

max_high_cloud_pct = 0.3
 # 7 21 all 1deg

#date_list = [dt.datetime(2015, 7, 01) + dt.timedelta(x) for x in range(0, 60)]
#date_list = [dt.datetime(2015, 07, 01) + dt.timedelta(x) for x in range(0, 1)]
#date_list = [dt.datetime(2015, 7, 01), dt.datetime(2015, 8, 03)]
date_list = [dt.datetime(2015, 7, 17)]
keylist = [1]
keylist = None

for i, d in enumerate(date_list):
    date_list[i] = d.replace(tzinfo=pytz.UTC)

# trajectory_folder = '/home/disk/eos4/jkcm/Data/CSET/Lagrangian_project/Trajectories/'
goes_folder = '/home/disk/eos4/jkcm/Data/CSET/GOES/VISST_pixel'
plot_folder = '/home/disk/eos4/jkcm/Data/CSET/Lagrangian_project/GOES_plots/new'
netcdf_folder = '/home/disk/eos4/jkcm/Data/CSET/Lagrangian_project/GOES_netcdf/new'



vis_varnames = (
    'cloud_visible_optical_depth',
    'cloud_particle_size', 'cloud_lwp_iwp',
    'broadband_shortwave_albedo'
)
ir_varnames = (
    #        'temperature_sw',
    'latitude', 'longitude',
    #        'cloud_effective_temperature',
           'broadband_longwave_flux',
    #        'temperature_sir', 'temperature_ir', 'cloud_ir_emittance',
    'pixel_skin_temperature',
    'cloud_top_pressure', 'cloud_top_height',
    'cloud_top_temperature',
)
defined_ylims = {
    'Nd': (0, 500),
    'cloud_lwp_iwp': (0, 400),
    'cloud_visible_optical_depth': (0, 15),
    'pixel_skin_temperature': (288, 298),
}
pct_varnames = (
    'cloud_visible_optical_depth',
    'cloud_particle_size',
    'cloud_lwp_iwp',
    'cloud_top_pressure',
    'cloud_top_height',
    'cloud_top_temperature',
    'Nd',
)
mean_only_varnames = (
    'pixel_skin_temperature',
)
image_only_varnames = (
    'visible_count',
    'latitude',
    'longitude',
    'cloud_phase',
)
cloud_only_varnames = (
    'cloud_top_pressure',
    'cloud_top_height',
    'cloud_top_temperature',
    'cloud_visible_optical_depth',
    'cloud_particle_size',
    'cloud_lwp_iwp',
    'cloud_top_pressure',
    'cloud_top_height',
    'cloud_top_temperature',
    'Nd',
)

# lt = LoopTimer(len(date_list)*90)
# print('Restarting from {:%Y-%m-%d}'.format(date_list[0]))
# for di, date in enumerate(date_list):
#     for delta_index, label in pairs:
#         print('Working on date {}, {}'.format(date, label))
#         filename = os.path.join(utils.trajectory_dir,
#                                 'tdump'+date.strftime('%Y%m%dH%H%M'))
# for all trajectories...

files = glob.glob('/home/disk/eos4/jkcm/Data/CSET/Trajectories/*/analysis.UW_HYSPLIT_GFS.*.airmass_trajectories_500m_-48.txt')

#above was the old loop. new loop:
for filename in files:

#        if all((os.path.isfile(
#                os.path.join(
#                    plot_folder, flightname + '-' + flight_label + '-' +
#                    varname + '-' + label + '.png')))
#                for varname in
#                (vis_varnames + ir_varnames + ('cloud_fraction',))):
#            continue
        try:
            data_list, long_names_list, units_list = get_goes_from_multitrajectory(
                filename,
                goes_folder, vis_varnames + ('visible_count',),
                ir_varnames + ('cloud_phase',),
                delta_index=delta_index, lt=lt,
                keylist=keylist)
        except RuntimeError as e:
            print('caught in datelist loop')
            continue

        print('Finished GOES read')
        sys.stdout.flush()

        for i, (data, long_names, units) in enumerate(zip(
                data_list, long_names_list, units_list)):
            print('Working on trajectory {}'.format(i))
            if i != 0:
                print ('skipping {}'.foramt(i))
                continue
            sys.stdout.flush()

            date_name = "{:%Y%m%d}".format(date)
            traj_label = "{:02d}".format(i)

            long_names['cloud_particle_size'] = 'effective particle radius'

            out_data = xray.Dataset()
            out_data.attrs['Title'] = 'Pixel box cloud products'
            out_data.attrs['institution'] = (
                'Department of Atmospheric Sciences, University of Washington')
            out_data.attrs['contact'] = 'jkcm@uw.edu'
            out_data.attrs['VISST'] = (
                "NASA-Langley cloud and radiation products are produced using "
                "the VISST (Visible Infrared Solar-infrared Split-Window "
                "Technique), SIST (Solar-infrared Infrared Split-Window "
                "Technique) and SINT (Solar-infrared Infrared Near-Infrared "
                "Technique). REFERENCE:Minnis, P., et al. , 2008: Near-real "
                "time cloud retrievals from operational and research "
                "meteorological satellites. Proc. SPIE Europe Remote Sens. "
                "2008, Cardiff, Wales, ID, 15-18 September, 7107, No. 2, "
                "8pp.http://www-angler.larc.nasa.gov/site/doc-library/99-Minni"
                "s.etal.SPIE.abs.08.pdf. Additional references on calibration,"
                " validation can be found at http://www-pm.larc.nasa.gov -> "
                "publication")
            out_data.attrs['comments'] = (
                "This dataset contains statistics of VISST cloud and radiation"
                " products for CSET airmass trajectories. At each hour, the "
                "VISST pixel closest to the airmass is identified, and a box "
                "of a predefined pixel size is taken around this pixel. All "
                "statistics are performed inside that pixel box. Cloud "
                "statistics are performed using liquid cloud pixels only.")
            out_data.attrs['references'] = (
                "Minnis, P., et al. , 2008: Near-real "
                "time cloud retrievals from operational and research "
                "meteorological satellites. Proc. SPIE Europe Remote Sens. "
                "2008, Cardiff, Wales, ID, 15-18 September, 7107, No. 2, "
                "8pp.http://www-angler.larc.nasa.gov/site/doc-library/99-Minni"
                "s.etal.SPIE.abs.08.pdf. Additional references on calibration,"
                " validation can be found at http://www-pm.larc.nasa.gov -> "
                "publication")

            out_data.coords['time'] = (('time',), data['time'])
            out_data.coords['reference_time'] = data['time'][0]

            time_axis = matplotlib.dates.date2num(data.pop('time'))
            x_bin_edges = np.append(
                time_axis, time_axis[-1:] + 1. / 24) - 1. / 48
            for varname in (
                    vis_varnames + ir_varnames +
                    ('Nd', 'cloud_phase', 'visible_count')):
                if varname == 'cloud_lwp_iwp':
                    data['cloud_lwp'] = data['cloud_lwp_iwp']
                    varname = 'cloud_lwp'
                    units['cloud_lwp'] = units['cloud_lwp_iwp']
                    long_names['cloud_lwp'] = 'Liquid Water Path'
                fig, ax = plt.subplots(3)
                try:
                    ax[0].set_title(' '.join((
                        date_name.upper(), traj_label, long_names[varname],
                        '($' + units[varname].replace('_', ' ') +
                        '$)\nin warm low cloud')))
                except KeyError:
                    # no long_name or units for this var means no data
                    continue
                image = np.zeros((len(time_axis), delta_index * 2 + 1,
                                  delta_index * 2 + 1)) * np.nan
                for i in range(image.shape[0]):
                    image[i, :, :] = data[varname][i, :].reshape(
                        (delta_index * 2 + 1, delta_index * 2 + 1))
                out_data[varname] = (
                    ('time', 'y', 'x'), image,
                    {'units': units[varname],
                     'long_name': long_names[varname]})
                mean = np.zeros(time_axis.shape) * np.nan
                median = np.zeros(time_axis.shape) * np.nan
                for i in range(len(time_axis)):
                    too_much_ice_cloud = np.sum(
                        (data['cloud_phase'][i, :] == 2) |
                        (data['cloud_top_height'][i, :] > 4.))/float(
                        data['cloud_phase'][i, :].size) > max_high_cloud_pct
                    if varname in cloud_only_varnames:
                        use_me = (
                            (data['cloud_phase'][i, :] == 1) &
                            (data['cloud_top_height'][i, :] < 4.) &
                            (data['cloud_top_temperature'][i, :] > 273.15))
                    else:
                        use_me = np.ones(
                            data['cloud_phase'][i, :].shape, dtype=np.bool)
                    if use_me.any() and (not too_much_ice_cloud):
                        mean[i] = np.nanmean(data[varname][i, :][use_me])
                        median[i] = np.nanmedian(data[varname][i, :][use_me])
                    else:
                        mean[i], median[i] = np.nan, np.nan
                ax[0].plot(time_axis, mean,
                           'c-', label='mean')
                if ((varname in mean_only_varnames) or
                        (varname in pct_varnames)):
                    out_data[varname + '_mean'] = (
                        ('time',), mean, {
                            'units': units[varname],
                            'long_name': 'mean ' + long_names[varname]})
                if varname in pct_varnames:
                    out_data[varname + '_median'] = (
                        ('time',), median, {
                            'units': units[varname],
                            'long_name': 'median ' + long_names[varname]})
                ax[0].plot(time_axis, median,
                           'b-', label='median')
                percentile = {'25': [], '75': []}
                for p in percentile.keys():
                    for j in range(len(time_axis)):
                        too_much_ice_cloud = np.sum(
                            (data['cloud_phase'][i, :] == 2) |
                            (data['cloud_top_height'][i, :] > 4.))/float(
                            data['cloud_phase'][j, :].size) > \
                                max_high_cloud_pct
                        if varname in cloud_only_varnames:
                            use_me = (
                                (data['cloud_phase'][j, :] == 1) &
                                (data['cloud_top_height'][j, :] < 4.) &
                                (data['cloud_top_temperature'][j, :] > 273.15))
                        else:
                            use_me = np.ones(
                                data['cloud_phase'][i, :].shape, dtype=np.bool)
                        use_me = use_me & ~np.isnan(data[varname][j, :])
                        if use_me.any() and not too_much_ice_cloud:
                            percentile[p].append(
                                np.percentile(data[varname][j, :][
                                    use_me], int(p)))
                        else:
                            percentile[p].append(np.nan)
                    percentile[p] = np.asarray(percentile[p])
                    ax[0].plot(time_axis, percentile[p], 'b--',
                               label=p + 'th %ile',)
                    if varname in pct_varnames:
                        out_data[varname + '_{}th_percentile'.format(p)] = (
                            ('time',), percentile[p], {
                                'units': units[varname],
                                'long_name': ('{}th percentile '.format(p) +
                                              long_names[varname])})
                ax[0].legend(loc='center left', bbox_to_anchor=(1., 0.5))
                hist_time, hist_data = np.broadcast_arrays(time_axis[:, None],
                                                           data[varname])
                # don't need to check if in cloud_only_varnames for this one
                # because we plot all-pixel and cloud-only histograms
                use_me = ((data['cloud_phase'][:, :] == 1) &
                          (data['cloud_top_height'][:, :] < 4.) &
                          (data['cloud_top_temperature'][:, :] > 273.15))
                use_me = use_me & (~np.isnan(hist_data))
                hist_time_masked = hist_time[use_me]
                hist_data_masked = hist_data[use_me]
                if varname in defined_ylims and len(hist_time) > 0:
                    val_min, val_max = defined_ylims[varname]
                    val_min_masked, val_max_masked = defined_ylims[varname]
                else:
                    try:
                        val_min = np.nanmin(hist_data)
                        val_max = np.nanmax(hist_data)
                        val_min_masked = np.nanmin(hist_data_masked)
                        val_max_masked = np.nanmax(hist_data_masked)
                    except ValueError:
                        val_min = np.nan
                        val_max = np.nan
                        val_min_masked = np.nan
                        val_max_masked = np.nan
                if (~np.isnan(val_min) and ~np.isnan(val_max) and
                        val_max != val_min):
                    _, _, _, im = ax[1].hist2d(
                        hist_time[~np.isnan(hist_data)],
                        hist_data[~np.isnan(hist_data)],
                        bins=(x_bin_edges, 15),
                        cmap=viridis, range=(
                            (time_axis[0], time_axis[-1]),
                            (val_min, val_max)))
                    cax = fig.add_axes([0.8, 0.14, 0.03, 0.45])
                    cbar = plt.colorbar(im, cax=cax)
                    cbar.set_ticklabels([])
                    _, _, _, im = ax[2].hist2d(
                        hist_time_masked, hist_data_masked,
                        bins=(x_bin_edges, 15),
                        cmap=viridis, range=(
                            (time_axis[0], time_axis[-1]),
                            (val_min_masked, val_max_masked)))
                for j in range(3):
                    ax[j].xaxis_date()
                    ax[j].set_xlim(time_axis[0], time_axis[-1])
                ax[0].set_xticklabels([])
                ax[1].set_xticklabels([])
                ax[1].set_title('Unmasked')
                ax[2].set_title('Warm low cloud pixels only')
                ax[2].xaxis.set_major_formatter(
                    matplotlib.dates.DateFormatter('%b %d\n%H:%M'))
                ax[2].set_xlabel('Time (UTC)')
                plt.tight_layout(rect=(0, 0, 0.8, 1))
                fig.savefig(
                    os.path.join(
                        plot_folder, date_name + '-' + traj_label + '-' +
                        varname + '-' + label + '.png'))
                plt.close(fig)
            out_data['trajectory_latitude'] = (
                ('time',), data['lat_traj'], {
                    'units': 'degrees_north',
                    'long_name': 'trajectory latitude',
                    'valid_min': -90.,
                    'valid_max': 90.})
            out_data['trajectory_longitude'] = (
                ('time',), data['lon_traj'], {
                    'units': 'degrees_east',
                    'long_name': 'trajectory longitude',
                    'valid_min': -180.,
                    'valid_max': 180.})
            # plot cloud fraction, treated separately
            filename = os.path.join(
                plot_folder, date_name + '-' + traj_label +
                '-cloud_fraction-' + label + '.png')
#            if not os.path.isfile(filename):  # this was stupid.
            fig, ax = plt.subplots(1)
            ax.set_title(' '.join((
                date_name.upper(), traj_label,
                'cloud fraction ($unitless$)')))
            cloud_fraction = np.zeros(
                data['cloud_phase'].shape[:1]) * np.nan
            for i in range(len(cloud_fraction)):
                cloud_fraction[i] = 1 - np.sum(
                    data['cloud_phase'][i, :] == 4) / float(np.sum(
                        (data['cloud_phase'][i, :] != 3) &  # no retrieval
                        (data['cloud_phase'][i, :] != 5)))  # bad data
            out_data['cloud_fraction'] = (
                ('time',), cloud_fraction, {
                    'long_name': 'cloud fraction'})
            ax.plot(time_axis, cloud_fraction, 'b-')
            ax.xaxis_date()
            ax.set_xlim(time_axis[0], time_axis[-1])
            ax.xaxis.set_major_formatter(
                matplotlib.dates.DateFormatter('%b %d\n%H:%M'))
            ax.set_xlabel('Time (UTC)')
            ax.set_ylim(0, 1)

            plt.tight_layout()
            fig.savefig(
                filename)
            plt.close(fig)
            cold_cloud_fraction = np.zeros(data['cloud_phase'].shape[:1])
            warm_cloud_fraction = np.zeros(data['cloud_phase'].shape[:1])
            warm_low_cloud_fraction = np.zeros(data['cloud_phase'].shape[:1])
            unmasked_pixel_fraction = np.zeros(data['cloud_phase'].shape[:1])

            for i in range(len(cold_cloud_fraction)):
                valid_pixel_count = float(np.sum(
                        (data['cloud_phase'][i, :] != 3) &  # no retrieval
                        (data['cloud_phase'][i, :] != 5) &  # bad data
                        (data['cloud_phase'][i, :] != 6) &  # suspected water
                        (data['cloud_phase'][i, :] != 7) &  # suspected ice
                        (data['cloud_phase'][i, :] != 13)))  # cleaned_data
                cold_cloud_fraction[i] = np.sum(
                    data['cloud_phase'][i, :] == 2) / valid_pixel_count
                warm_cloud_fraction[i] = np.sum(
                    data['cloud_phase'][i, :] == 1) / float(np.sum(  # liquid
                        (data['cloud_phase'][i, :] == 4) |  # clear
                        (data['cloud_phase'][i, :] == 1)))  # water
                warm_low_cloud = (
                    (data['cloud_phase'][i, :] == 1) &
                    (data['cloud_top_height'][i, :] < 4.) &
                    (data['cloud_top_temperature'][i, :] > 273.15))
                warm_low_cloud_fraction[i] = np.sum(
                    warm_low_cloud) / float(np.sum(
                        (data['cloud_phase'][i, :] == 4) |  # clear
                        warm_low_cloud))
                unmasked_pixel_fraction[i] = np.sum(warm_low_cloud)/float(
                    data['cloud_phase'].shape[1])
            for fraction_data, f_label in (
                    (warm_cloud_fraction, 'warm cloud fraction\n'
                     'with random overlap assumption'),
                    (cold_cloud_fraction, 'cold cloud fraction'),
                    (warm_low_cloud_fraction, 'warm low cloud fraction\n'
                     'with random overlap assumption'),
                    (unmasked_pixel_fraction,
                     'warm low cloud pixel fraction')):
                varname = f_label.replace(
                    ' ', '_').replace(
                        '\nwith_random_overlap_assumption', '')
                filename = os.path.join(
                    plot_folder, date_name + '-' + traj_label +
                    '-' + varname +
                    '-' + label + '.png')
                fig, ax = plt.subplots(1)
                ax.set_title(' '.join((
                    date_name.upper(), traj_label,
                    '{} ($unitless$)'.format(varname.replace('_', ' ')))))
                ax.plot(time_axis, fraction_data, 'b-')
                ax.xaxis_date()
                ax.set_xlim(time_axis[0], time_axis[-1])
                ax.xaxis.set_major_formatter(
                    matplotlib.dates.DateFormatter('%b %d\n%H:%M'))
                ax.set_xlabel('Time (UTC)')
                ax.set_ylim(0, 1)
                plt.tight_layout()
                fig.savefig(filename)
                plt.close(fig)
                out_data[varname] = (
                    ('time',), fraction_data, {'long_name': f_label})
            out_data['cloud_phase'].attrs['value_0'] = 'snow'
            out_data['cloud_phase'].attrs['value_1'] = 'water'
            out_data['cloud_phase'].attrs['value_2'] = 'ice'
            out_data['cloud_phase'].attrs['value_3'] = 'no retrieval'
            out_data['cloud_phase'].attrs['value_4'] = 'clear'
            out_data['cloud_phase'].attrs['value_5'] = 'bad data'
            out_data['cloud_phase'].attrs['value_6'] = 'suspected water'
            out_data['cloud_phase'].attrs['value_7'] = 'suspected ice'
            out_data['cloud_phase'].attrs['value_13'] = 'cleaned data'
            use_me = ((data['cloud_phase'][:, :] == 1) &
                      (data['cloud_top_height'][:, :] < 4.) &
                      (data['cloud_top_temperature'][:, :] > 273.15))
            image = np.zeros((len(time_axis), delta_index * 2 + 1,
                              delta_index * 2 + 1)) * np.nan
            for i in range(image.shape[0]):
                image[i, :, :] = use_me[i, :].reshape(
                    (delta_index * 2 + 1, delta_index * 2 + 1))
            out_data['warm_low_cloud_flag'] = (
                ('time', 'y', 'x'), image,
                {'description': 'True corresponds to warm low cloud'})
            out_data.to_netcdf(
                os.path.join(
                    netcdf_folder, date_name + '-' + traj_label + '-' +
                    label + '.nc'))
