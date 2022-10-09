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
parser.add_argument('-use_random_delay', default=False, action='store_true',
                    help='''If set to True, use a randomized repoint delay. Otherwise, specific
                    ground station pass times from an orbital file are used to determine repoint
                    times. Default is False.''')
parser.add_argument('-target_method', default='mcintosh', type=str,
                    choices=['mcintosh', 'mcintosh_bloomfield', 'mcintosh_bornmann_shaw', 'old_mcintosh', 'hale', 'random'],
                    help='select methodology for choosing AR targets')
parser.add_argument('--delay_distribution', default='12hr', type=str,
                    choices=['12hr', '24hr', '1hr', '36hr'],
                    help='Select from defined delay distributions')
parser.add_argument('--ignore_eclipses', default=False, action='store_true',
                    help='Ignore eclipses/SAA when determining if a flare has been observed')
parser.add_argument('--fallback_target', default=0, type=int,
                    help='The fallback target to use, which should typically be <= 0 (i.e., not a NOAA number)')
parser.add_argument('--rank_thresholds', default='np.arange(0,1,0.1)', type=str,
                    help='A string of Python code without whitespace that specifies the list of thresholds to use')

args = parser.parse_args()

# Set up the logging
log = fm.setup_logging()

# Interpret the rank_thresholds argument
exec(f'rank_thresholds = {args.rank_thresholds}')
if isinstance(rank_thresholds, str):
    exec(f'rank_thresholds = {rank_thresholds}')

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

total_xflare_count = np.sum(
    flare_data['goes_class'][fm.PHASE_E[0]: fm.PHASE_E[1]].str.contains('X'))

total_mflare_count = np.sum(
    flare_data['goes_class'][fm.PHASE_E[0]: fm.PHASE_E[1]].str.contains('M'))

M_mean = []
M_stdev = []
M_median = []
X_mean = []
X_stdev = []
X_median = []
MX_mean = []
MX_stdev = []
MX_median = []
EW_mean = []
EW_stdev = []
EW_median = []


for rank_threshold in rank_thresholds:

    log.info(f"Strategy = {args.target_method}, "
             f"Rank threshold = {rank_threshold}, "
             f"Fallback target = {args.fallback_target}")

    # generate AR targets as defined by the SOC
    ar_targets = ms.generate_ar_targets(ar_data[fm.PHASE_E[0]: fm.PHASE_E[1]], mfw=mfw,
                                        target_method=args.target_method,
                                        rank_threshold=rank_threshold,
                                        fallback_target=args.fallback_target)

    xflare_observing_fractions = []
    mflare_observing_fractions = []

    xflare_observing_totals = []
    mflare_observing_totals = []
    mx_observing_totals = []
    ew_ratios = []

    for i in range(0, args.n):
        log.info('Performing sample {:n} out of {:n}.'.format(i+1, args.n))

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


        # Count the flares in view, and the pointing for the flares in view
        flares_in_view_count, flares_in_view_pointing, flares_in_clean_view_pointing = ms.count_flares_in_view(flare_to_pointing)

        total_xflares_observed = flares_in_view_count['on detector and observed']['X']
        total_mflares_observed = flares_in_view_count['on detector and observed']['M']

        xflare_observing_totals.append(total_xflares_observed)
        mflare_observing_totals.append(total_mflares_observed)
        mx_observing_totals.append(total_xflares_observed + total_mflares_observed)
        east_mx = np.sum(flares_in_view_pointing['flare_hpc_x'] < 0)
        west_mx = np.sum(flares_in_view_pointing['flare_hpc_x'] > 0)
        ew_ratios.append(east_mx / west_mx)


    M_mean.append(np.mean(mflare_observing_totals) / total_mflare_count)
    M_stdev.append(np.std(mflare_observing_totals, ddof=1) / total_mflare_count)
    M_median.append(np.median(mflare_observing_totals) / total_mflare_count)
    X_mean.append(np.mean(xflare_observing_totals) / total_xflare_count)
    X_stdev.append(np.std(xflare_observing_totals, ddof=1) / total_xflare_count)
    X_median.append(np.median(xflare_observing_totals) / total_xflare_count)
    MX_mean.append(np.mean(mx_observing_totals) / (total_xflare_count + total_mflare_count))
    MX_stdev.append(np.std(mx_observing_totals, ddof=1) / (total_xflare_count + total_mflare_count))
    MX_median.append(np.median(mx_observing_totals) / (total_xflare_count + total_mflare_count))
    EW_mean.append(np.mean(ew_ratios))
    EW_stdev.append(np.std(ew_ratios, ddof=1))
    EW_median.append(np.median(ew_ratios))
    
timestamp = datetime.datetime.now()
timestamp_s = timestamp.strftime(fm.now_str_format)
log.info('Use the timestamp "{:s}" to load the output scan data.'.format(timestamp_s))

scan = pd.DataFrame({'rank_threshold': rank_thresholds,
                     'M_mean': M_mean,
                     'M_stdev': M_stdev,
                     'M_median': M_median,
                     'X_mean': X_mean,
                     'X_stdev': X_stdev,
                     'X_median': X_median,
                     'MX_mean': MX_mean,
                     'MX_stdev': MX_stdev,
                     'MX_median': MX_median,
                     'EW_mean': EW_mean,
                     'EW_stdev': EW_stdev,
                     'EW_median': EW_median,
                     })
scan.to_csv(fm.in_output_dir('rank_threshold_scan', 'csv', timestamp=timestamp))
