import datetime
import os.path

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

#from flare_mission_sim import flare_mission_sim as ms
import flare_mission_sim as ms
#import flare_mission_sim as fm
import __init__ as fm
# smartex phase E
#phase_e = [datetime.datetime(2025, 8, 1), datetime.datetime(2027, 8, 1)]

# foxsi phase E
#phase_e = [datetime.datetime(2022, 8, 1), datetime.datetime(2024, 8, 1)]

# lunanet launch
#phase_e = [datetime.datetime(2021, 12, 1), datetime.datetime(2022, 12, 1)]

# fierce phase e extension

phase_e = [datetime.datetime(2025, 8, 1), datetime.datetime(2027, 8, 1)]
threshold_phase_e = [datetime.datetime(2025, 8, 1), datetime.datetime(2027, 2, 1)]



datafile_url = os.path.join(fm.data_dir, 'soda_predict.csv')
soda_predict = pd.read_csv(datafile_url, index_col=0, parse_dates=True,
                           names=('times', 'index'))

flares = ms.load_flare_data_from_hdf5()
flares = ms.shift_by_solar_cycle(flares)

xlim = [datetime.datetime(2020, 1, 1),
        datetime.datetime(2030, 1, 1)]

class_list = ['X', 'M', 'C']

class_index = {this_class: flares['goes_class'].str.contains(this_class)
               for this_class in class_list}

flare_count = pd.DataFrame({'X': flares['goes_class'][class_index['X']].resample('M').count(),
                            'M': flares['goes_class'][class_index['M']].resample('M').count(),
                            'C': flares['goes_class'][class_index['C']].resample('M').count()})
flare_count.fillna(0, inplace=True)

times = flare_count['C'].index

total_flares_during_phase_e = {'X': np.sum(flare_count['X'][phase_e[0]:phase_e[1]]),
                               'M': np.sum(flare_count['M'][phase_e[0]:phase_e[1]]),
                               'C': np.sum(flare_count['C'][phase_e[0]:phase_e[1]])}

f = plt.figure()
ax = f.add_subplot(111)
ax.bar(times, flare_count['C'].values, width=1, label='C flares')
ax.bar(times, flare_count['M'].values, width=1, label='M flares', bottom=flare_count['C'].values)
ax.bar(times, flare_count['X'].values, width=1, label='X flares', bottom=flare_count['C'].values + flare_count['M'].values)
ax.set_ylabel('Number of flares (per month)')
plt.legend(loc=2)
plt.axvline(phase_e[0], color='red', label='launch')
plt.axvspan(phase_e[0], phase_e[1], alpha=0.2, color='grey', label='phase e')
plt.axvline(threshold_phase_e[1], color='grey', label='threshold mission end')
plt.legend(loc=1)
plt.title(f'Total flares in Phase E: {int(total_flares_during_phase_e["X"])}X, '
          f'{int(total_flares_during_phase_e["M"])}M, '
          f'{int(total_flares_during_phase_e["C"])}C')

ax2 = ax.twinx()
# adding a weird shift to make it line up with the prediction of the peak, need to
# get correct soda prediction index
ax2.plot(soda_predict.index + datetime.timedelta(days=220), soda_predict.values, label='SODA index', linestyle='--', color='grey')
ax2.set_ylim(0, 150)
ax2.set_ylabel('SODA index')

# now plot the prediction of the peak with errors
soda_x = datetime.datetime(2025, 1, 1) + datetime.timedelta(days=0.2*365.5)
soda_y = np.array([140])
ax2.plot([[soda_x-datetime.timedelta(days=1.5*365.5)], [soda_x+datetime.timedelta(days=1.5*365.5)]], [140, 140], color='grey')
ax2.plot([soda_x, soda_x], [110, 170], color='grey')

plt.legend(loc=2)
plt.xlim(xlim[0], xlim[1])
plt.savefig(fm.in_output_dir('mission_solar_cycle', 'pdf', timestamp=False))
