import glob
import os
import argparse
import pandas as pd

#from flare_mission_sim.plot_mc_outputs_new_heatmap import *
from plot_mc_outputs_new_heatmap import *
#import flare_mission_sim as fm
import __init__ as fm
#from flare_mission_sim import flare_mission_sim as ms
import flare_mission_sim as ms

# Get the timestamp from the command line
parser = argparse.ArgumentParser(description='Make plots of a simulation run.')
parser.add_argument('--timestamp', '-t', default=None, type=str,
                    help='time stamp denoting a simulation run in the format YYYYMMDD_hhmm')
parser.add_argument('--create_heatmaps', default=False, action='store_true', help='create heatmap plots')
args = parser.parse_args()
timestamp = args.timestamp
if timestamp is None:
    raise ValueError('A timestamp with format "YYYYMMDD_hhmm" must be set.')

# find all of the files
files = glob.glob(os.path.join(fm.out_dir, 'xflare_observing_fractions_*.txt'))


def output_filepath(dataname, extension):
    filename = f'{dataname}_{fm.MISSION_NAME}_{timestamp}.{extension}'
    return os.path.join(fm.out_dir, filename)


#xflare_observing_fractions = np.loadtxt(output_filepath('xflare_observing_fractions', 'txt'))
#mflare_observing_fractions = np.loadtxt(output_filepath('mflare_observing_fractions', 'txt'))
#xflare_observing_totals = np.loadtxt(output_filepath('xflare_observing_totals', 'txt'))
#mflare_observing_totals = np.loadtxt(output_filepath('mflare_observing_totals', 'txt'))


#plot_flare_observing_total_histograms(xflare_observing_totals, mflare_observing_totals)
#plot_observing_efficiency_histogram(xflare_observing_fractions, mflare_observing_fractions)

# load flare and active region data
#ar_data = ms.load_ar_data()
#flare_data = ms.load_flare_data_from_hdf5()

# remove all C class flares
#flare_data = flare_data[flare_data['goes_class'].str.contains('X|M')]

#flare_data = ms.fix_flare_data(flare_data, ar_data)
#flare_data = ms.shift_by_solar_cycle(flare_data)
#ar_data = ms.shift_by_solar_cycle(ar_data)

#flare_data = flare_data[fm.PHASE_E[0]: fm.PHASE_E[1]]

#xtotal_flare_count = np.sum(flare_data['goes_class'].str.contains('X'))
#mtotal_flare_count = np.sum(flare_data['goes_class'].str.contains('M'))


#plot_pointing_efficiency_histogram(xflare_observing_totals / xtotal_flare_count * 100.0, mflare_observing_totals / mtotal_flare_count * 100.0)

# Plot where the mission observes a flare
heatmap_shape = (2200, 2200)

filepath = output_filepath('phase_e_flares', 'csv')
print(f'Loading {filepath}')
phase_e_flares = pd.read_csv(filepath)

filepath = output_filepath('flares_in_view', 'csv')
print(f'Loading {filepath}')
flares_in_view = pd.read_csv(filepath)

plot_map_of_observed_flares(phase_e_flares, flares_in_view, heatmap_shape, 'test')


if args.create_heatmaps:

    # Mission pointing heatmaps
    filepath = output_filepath('mission_pointings', 'csv')
    print(f'Loading {filepath}')
    mission_pointings = pd.read_csv(filepath)
    heatmap_mission_pointings_central = ms.make_heatmap(mission_pointings, fm.central_fov, heatmap_shape)
    plot_heatmaps_and_summaries(heatmap_mission_pointings_central,
                                'heat_mission_pointing_heat_using_central_FOV',
                                'mission pointing heat using central FOV')
    heatmap_mission_pointings_full = ms.make_heatmap(mission_pointings, fm.full_fov, heatmap_shape)
    plot_heatmaps_and_summaries(heatmap_mission_pointings_full,
                                'heat_mission_pointing_heat_using_full_FOV',
                                'mission pointing heat using full FOV')

    # Plot a FOV for all the flares in Phase E whether we saw them or not
    filepath = output_filepath('phase_e_flares', 'csv')
    print(f'Loading {filepath}')
    phase_e_flares = pd.read_csv(filepath)
    heatmap_phase_e_flares_central = ms.make_heatmap(phase_e_flares, fm.central_fov, heatmap_shape)
    plot_heatmaps_and_summaries(heatmap_phase_e_flares_central,
                                'heat_all_phase_E_flares_represented_by_a_central_FOV',
                                'all phase E flares (each represented by a central FOV')

    # Plot a FOV for all the flares observed by the mission
    filepath = output_filepath('flares_in_view', 'csv')
    print(f'Loading {filepath}')
    flares_in_view = pd.read_csv(filepath)
    heatmap_flares_in_view_full = ms.make_heatmap(flares_in_view, fm.full_fov, heatmap_shape)
    plot_heatmaps_and_summaries(heatmap_flares_in_view_full,
                                'heat_all_flares_observed_by_mission_each_represented_by_a_full_FOV',
                                'all flares observed by mission (each represented by a full FOV)')
