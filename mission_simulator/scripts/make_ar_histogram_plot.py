import numpy as np
import datetime

import matplotlib.pyplot as plt

#from flare_mission_sim import flare_mission_sim as ms
import flare_mission_sim as ms
#import flare_mission_sim as fm
import __init__ as fm

# load flare and active region data
ar_data = ms.load_ar_data()
ar_data = ms.shift_by_solar_cycle(ar_data)

# smartex phase E
phase_e = [datetime.datetime(2025, 8, 1), datetime.datetime(2027, 8, 1)]

this_ar_data = ar_data[phase_e[0]:phase_e[1]]

bins = np.arange(0, np.max(this_ar_data['number']), 1)
hist, bins = np.histogram(this_ar_data['number'].values, bins=bins)
center = (bins[:-1] + bins[1:]) / 2 - 0.5

plt.figure(1, figsize=(12, 4))

plt.subplot(1, 2, 1)

plt.bar(center, hist, align='center', width=1)
plt.title(f'{phase_e[0]:%Y/%m/%d} to {phase_e[1]:%Y/%m/%d}')
plt.ylabel('Number of Days')
plt.xlabel('Number of Active regions')
plt.xticks(bins)
plt.grid(False)

plt.subplot(1, 2, 2)
plt.pie(hist, autopct='%.1f', labels=bins[:-1])

plt.savefig(fm.in_output_dir('active_region_histogram', 'pdf'))
