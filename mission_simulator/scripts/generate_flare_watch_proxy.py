import numpy as np
from flare_mission_sim import major_flare_watch_proxy

#from flare_mission_sim import flare_mission_sim as ms
import flare_mission_sim as ms
#import flare_mission_sim as fm
import __init__ as ms

ar_data = ms.load_ar_data()
flare_data = ms.load_flare_data_from_hdf5()

# remove all C class flares
flare_data = flare_data[flare_data['goes_class'].str.contains('X|M')]

flare_data = ms.fix_flare_data(flare_data, ar_data)
flare_data = ms.shift_by_solar_cycle(flare_data)
ar_data = ms.shift_by_solar_cycle(ar_data)

flare_data = flare_data[fm.PHASE_E[0]: fm.PHASE_E[1]]

mfw_series = major_flare_watch_proxy.generate_mfw_proxy(flare_data, ar_data)

