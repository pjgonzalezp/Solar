import numpy as np

from flare_mission_sim import flare_mission_sim as ms
import flare_mission_sim as fm

# read in orbit events
eclipse, saa, contacts, polar = ms.read_orbit_events()

not_observing = ms.concatenate_range_series(eclipse, saa)
if polar is not None:
    not_observing = ms.concatenate_range_series(not_observing, polar)

# read in major flare watches
mfw = ms.load_major_flare_watches(tshift=False)

# load flare and active region data
ar_data = ms.load_ar_data()
flare_data = ms.load_flare_data_from_hdf5()

# remove all C class flares
flare_data = flare_data[flare_data['goes_class'].str.contains('X|M')]

flare_data = ms.fix_flare_data(flare_data, ar_data)
flare_data = ms.shift_by_solar_cycle(flare_data)
ar_data = ms.shift_by_solar_cycle(ar_data)

flare_data = flare_data[fm.PHASE_E[0]: fm.PHASE_E[1]]
#ar_data = ar_data[fm.PHASE_E[0]: fm.PHASE_E[1]]

# generate SOC decision only for the ARs during phase E
soc_decisions = ms.generate_pointing_decisions(ar_data[fm.PHASE_E[0]: fm.PHASE_E[1]], mfw=mfw)
soc_decisions.to_csv(fm.in_output_dir('soc_decisions', 'csv'))

# generate pointing change commands
pointing_commands = ms.generate_pointing_commands(soc_decisions, contacts)
pointing_commands.to_csv(fm.in_output_dir('pointing_changes', 'csv'))

# compare pointing to flare positions and observable times
flare_to_pointing = ms.get_flare_observations(flare_data, pointing_commands,
                                              ar_data, not_observing)
flare_to_pointing.to_csv(fm.in_output_dir('flare_to_pointing', 'csv'))

flares_in_view = ms.count_flares_in_view(flare_to_pointing)
flares_in_view.to_csv(fm.in_output_dir('flares_in_view_summary', 'csv'))

for this_class in ['X', 'M']:
    total_flare_count = np.sum(flare_data[fm.PHASE_E[0]:fm.PHASE_E[1]]['goes_class'].str.contains(this_class))
    print("The total number of {0} flares during this time period is {1}".format(this_class, total_flare_count))
    print("The pointing efficiency is {0:0.2f}%".format(flares_in_view['in clean fov'][this_class] / total_flare_count * 100))
