import argparse

import pandas as pd
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser(description='Make plots of a rank-threshold scan')
parser.add_argument('--timestamp', default=None, type=str,
                    help='time stamp denoting a scan run in the format YYYYMMDD_hhmm')
parser.add_argument('--title', default=None, type=str,
                    help='title')
args = parser.parse_args()

data = pd.read_csv(f'out/rank_threshold_scan_foxsi_{args.timestamp}.csv')

fig, ax = plt.subplots(1, 2, figsize=(8, 4), constrained_layout=True)

ax[0].errorbar(data['rank_threshold'], data['MX_mean'], yerr=data['MX_stdev'], fmt='o')
ax[0].set_xlabel('Threshold for ranking')
ax[0].set_ylabel('M+X capture fraction')

ax[1].errorbar(data['rank_threshold'], data['EW_mean'], yerr=data['EW_stdev'], fmt='o')
ax[1].set_xlabel('Threshold for ranking')
ax[1].set_ylabel('E/W ratio')

fig.suptitle(args.title)

plt.show()
