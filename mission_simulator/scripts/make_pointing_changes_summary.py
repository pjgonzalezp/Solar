import matplotlib.pyplot as plt
import pandas as pd
import os
#import flare_mission_sim as fm
import __init__ as fm
import astropy.units as u
import numpy as np

f = os.path.join(fm.out_dir, 'pointing_changes_fierce_20190909_1411.csv')
data = pd.read_csv(f, parse_dates=True, index_col=0)

dt_pointing = u.Quantity((data.index[1:] - data.index[0:-1]).total_seconds(), 's').to('day')

binwidth = 0.2
bins = np.arange(min(dt_pointing.value), max(dt_pointing.value) + binwidth, binwidth)

plt.figure(1, figsize=(12, 4))

plt.subplot(1, 2, 1)
plt.hist(dt_pointing.value, bins=bins, density=True)
mn = np.mean(dt_pointing.value)
plt.axvline(mn, color='black', linestyle='dashed', label=f'mean={mn:.1f}')
plt.xlabel('Dt pointing changes [days]')
plt.ylabel('Number')
plt.legend()

plt.subplot(1, 2, 2)
plt.hist(dt_pointing.value, bins=bins, density=True, cumulative=True)
mn = np.mean(dt_pointing.value)
plt.axvline(mn, color='black', linestyle='dashed', label=f'mean={mn:.1f}')
plt.xlabel('Dt pointing changes [days]')
plt.ylabel('Number')
plt.legend()

plt.savefig(fm.in_output_dir('pointing_changes', 'pdf'))
plt.close()
