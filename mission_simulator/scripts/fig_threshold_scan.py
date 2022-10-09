# Run the following script.  Here N=100, which takes ~350 minutes on my computer.
"""
%run scripts/scan_rank_thresholds.py\
    -target_method=mcintosh_bloomfield\
    -n=100\
    -use_random_delay\
    --delay_distribution=12hr\
    --ignore_eclipses\
    --fallback_target=-110\
    --rank_thresholds='np.concatenate([np.arange(0,0.4,0.02),np.arange(0.4,1.4,0.05),np.arange(1.5,1.8,0.1)])'
"""

import matplotlib.pyplot as plt
import pandas as pd

# Use the timestamp that is in the last line of printed output from the earlier script
data = pd.read_csv(f'out/rank_threshold_scan_foxsi_20220906_1640.csv')

# The script produces a E/W ratio, so this converts it back to a East fraction and West fraction
ratio = data['EW_mean']
east = ratio / (1 + ratio)
east = east.fillna(1)
west = 1 / (1 + ratio)

plt.figure(figsize=(6, 6))

plt.plot(data['rank_threshold'], data['MX_mean'], label='Total M+X observed')
plt.plot(data['rank_threshold'], east * data['MX_mean'], label='East M+X observed')
plt.plot(data['rank_threshold'], west * data['MX_mean'], label='West M+X observed')

plt.grid()

plt.title('McIntosh Bloomfield, 12-hr delay, no eclipses, fallback is NE')
plt.xlabel('Threshold for flare productivity (flares/day)')
plt.ylabel('Fraction of all M+X flares')
plt.legend()

plt.savefig('threshold_scan.pdf')
