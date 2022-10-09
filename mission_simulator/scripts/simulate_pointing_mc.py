import datetime
import argparse
import numpy as np
import pandas as pd

#import flare_mission_sim as fm
import __init__ as fm
#from flare_mission_sim import flare_mission_sim as ms
import flare_mission_sim as ms

# Get the number of simulations to run from the command line
parser = argparse.ArgumentParser(description='Simulate the observatory pointing.')
parser.add_argument('-n', default=1, type=int,
                    help='number of simulations to run')
parser.add_argument('--save_heatmap_info', default=False, action='store_true',
                    help='save output to create heatmaps')
parser.add_argument('-use_random_delay', default=False, action='store_true',
                    help='''If set to True, use a randomized repoint delay. Otherwise, specific
                    ground station pass times from an orbital file are used to determine repoint
                    times. Default is False.''')
parser.add_argument('-target_method', default='mcintosh', type=str,
                    choices=['mcintosh', 'mcintosh_bloomfield', 'mcintosh_bornmann_shaw', 'old_mcintosh', 'hale', 'random', 'flare_index'],
                    help='select methodology for choosing AR targets')
parser.add_argument('--delay_distribution', default='12hr', type=str,
                    choices=['12hr', '24hr', '1hr', '36hr'],
                    help='Select from defined delay distributions')
parser.add_argument('--ignore_eclipses', default=False, action='store_true',
                    help='Ignore eclipses/SAA when determining if a flare has been observed')
parser.add_argument('--rank_threshold', default=0, type=float,
                    help='Minimum threshold for the numerical ranking, otherwise the target fallback is used')
parser.add_argument('--fallback_target', default=0, type=int,
                    help='The fallback target to use, which should typically be <= 0 (i.e., not a NOAA number)')

args = parser.parse_args()

# Set up the logging
log = fm.setup_logging()

# Load the delay parameters
delay_params_dict = {'12hr': [10., 2., 2.],
                     '24hr': [20., 4., 2.],
                     '1hr' : [0.0, 1.0, 0.5],
                     '36hr': [34., 2., 2.]}
delay_params = delay_params_dict[args.delay_distribution]

# Read in orbit events
eclipse, saa, contacts, polar = ms.read_orbit_events()

# Shift orbit data back to solar cycle 24, where we have real flare and AR data
eclipse = ms.shift_by_solar_cycle(eclipse, n=-1)
saa = ms.shift_by_solar_cycle(saa, n=-1)
contacts = ms.shift_by_solar_cycle(contacts, n=-1)
eclipse_and_saa = ms.concatenate_range_series(eclipse, saa)

# read in major flare watches
mfw = ms.load_major_flare_watches(tshift=False)

# load flare and active region data
ar_data = ms.load_ar_data()

# Returns a Pandas dataframe.  Note that the index is the start time of the flare.
flare_data = ms.load_flare_data_from_hdf5()

# remove all C class flares
flare_data = flare_data[flare_data['goes_class'].str.contains('X|M')]
# fix bad flare positions by cross-referencing with SFF and manual databases
flare_data = ms.fix_flare_positions_sff(flare_data)
flare_data = ms.fix_flare_positions_manual(flare_data)
#flare_data = ms.shift_by_solar_cycle(flare_data)
#ar_data = ms.shift_by_solar_cycle(ar_data)

#flare_data = flare_data[fm.PHASE_E[0]: fm.PHASE_E[1]]

# generate AR targets as defined by the SOC
if not args.target_method == 'random':
    ar_targets = ms.generate_ar_targets(ar_data[fm.PHASE_E[0]: fm.PHASE_E[1]], mfw=mfw,
                                        target_method=args.target_method, rank_threshold=args.rank_threshold,
                                        fallback_target=args.fallback_target,
                                        allowed_days_of_week = [1,2,3,4,5,6,7])
    

locations = []
xflare_observing_fractions = []
mflare_observing_fractions = []

xflare_observing_totals = []
mflare_observing_totals = []

xflare_ar_matches_target_ar_totals = []
mflare_ar_matches_target_ar_totals = []

total_xflare_count = np.sum(
    flare_data['goes_class'][fm.PHASE_E[0]: fm.PHASE_E[1]].str.contains('X'))

total_mflare_count = np.sum(
    flare_data['goes_class'][fm.PHASE_E[0]: fm.PHASE_E[1]].str.contains('M'))

# Define the size of the pointing heatmap.  Dimensions are understood to be in arcseconds
pointing_heatmap_full = np.zeros((2200, 2200))
pointing_heatmap_central = np.zeros_like(pointing_heatmap_full)
pointing_heatmap_time_step = pd.Timedelta(24, 'hour')
nsamples = 0

for i in range(0, args.n):
    log.info('Performing sample {:n} out of {:n}.'.format(i+1, args.n))

    # for the random method only, recalculate the ar targets each iteration
    if args.target_method == 'random':
        ar_targets = ms.generate_ar_targets(ar_data[fm.PHASE_E[0]: fm.PHASE_E[1]], mfw=mfw,
                                        target_method=args.target_method, rank_threshold=args.rank_threshold,
                                        fallback_target=args.fallback_target,
                                        allowed_days_of_week = [1,2,3,4,5,6,7])
        
    # Shuffle and concatenate eclipse and saa data
    eclipse_and_saa2, contacts2 = ms.shuffle_concatenated_series(eclipse_and_saa, contacts)

    # Generate pointing change commands
    pointing_commands = ms.generate_pointing_commands(ar_targets, contacts2, monte_carlo=True,
                                                      use_random_delay=args.use_random_delay,
                                                      delay_params=delay_params)

    # Compare pointing to flare positions and observable times
    flare_to_pointing = ms.get_flare_observations(flare_data, pointing_commands,
                                                  ar_data,
                                                  eclipse_and_saa2 if not args.ignore_eclipses else None)

    # Get the pointing locations
    # Once you have the pointing commands, use get_pointing to interpolate between the start and
    # end time of the pointing, tracking the motion of the target AR.
    # Use "pointing commands" to get how pointing is selected
    # Pass that to "get_pointing" along with the already known "ar_series" data
    # There is a function "plot_pointing" that might be useful - use it to create masks that
    # can be added up to create a "heat map" of pointing.
    if args.save_heatmap_info:
        log.info('Generating mission pointings')
        mission_pointings = ms.get_mission_pointings(pointing_commands, ar_data,
                                                     cadence=pointing_heatmap_time_step)

    # Count the flares in view, and the pointing for the flares in view
    flares_in_view_count, flares_in_view_pointing, flares_in_clean_view_pointing = ms.count_flares_in_view(flare_to_pointing)

    # How many times did the flare AR match the target AR?
    flare_ar_matches_target_ar = ms.count_flare_ar_matches_target_ar(flare_to_pointing)

    total_xflares_observed = flares_in_view_count['on detector and observed']['X']
    fraction_of_xflares_observed = total_xflares_observed / total_xflare_count

    total_mflares_observed = flares_in_view_count['on detector and observed']['M']
    fraction_of_mflares_observed = total_mflares_observed / total_mflare_count

    xflare_observing_fractions.append(fraction_of_xflares_observed)
    mflare_observing_fractions.append(fraction_of_mflares_observed)
    xflare_observing_totals.append(total_xflares_observed)
    mflare_observing_totals.append(total_mflares_observed)

timestamp = datetime.datetime.now()
timestamp_s = timestamp.strftime(fm.now_str_format)
log.info('Use the timestamp "{:s}" to load the output simulation data.'.format(timestamp_s))

np.savetxt(fm.in_output_dir('xflare_observing_fractions', 'txt', timestamp=timestamp), np.array(xflare_observing_fractions))
np.savetxt(fm.in_output_dir('mflare_observing_fractions', 'txt', timestamp=timestamp), np.array(mflare_observing_fractions))
np.savetxt(fm.in_output_dir('xflare_observing_totals', 'txt', timestamp=timestamp), np.array(xflare_observing_totals), fmt='%i')
np.savetxt(fm.in_output_dir('mflare_observing_totals', 'txt', timestamp=timestamp), np.array(mflare_observing_totals), fmt='%i')

# The heatmaps found for all the simulation runs
if args.save_heatmap_info:
    mission_pointings.to_csv(fm.in_output_dir('mission_pointings', 'csv', timestamp=timestamp))

# All the flare pointings found over all the simulation runs
# Rename the position columns so that the resulting CSV file works with the pointing heatmap code
# TODO: check that this renaming does not break anything else!
flares_in_view_pointing.rename(columns={"pointing_hpc_x": "hpc_x", "pointing_hpc_y": "hpc_y"}, inplace=True)
flares_in_view_pointing.to_csv(fm.in_output_dir('flares_in_view', 'csv', timestamp=timestamp))

# All the flares in the central (aka clean) FOV found over all the simulation runs
# Rename the position columns so that the resulting CSV file works with the pointing heatmap code
# TODO: check that this renaming does not break anything else!
flares_in_clean_view_pointing.rename(columns={"pointing_hpc_x": "hpc_x", "pointing_hpc_y": "hpc_y"}, inplace=True)
flares_in_clean_view_pointing.to_csv(fm.in_output_dir('flares_in_clean_view', 'csv', timestamp=timestamp))

# All the flares in the Phase E timerange
phase_e_flares = flare_data[fm.PHASE_E[0]: fm.PHASE_E[1]]
phase_e_flares.to_csv(fm.in_output_dir('phase_e_flares', 'csv', timestamp=timestamp))

