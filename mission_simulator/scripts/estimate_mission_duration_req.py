import datetime
import numpy as np
from scipy.stats import binom

#from flare_mission_sim import flare_mission_sim as ms
import flare_mission_sim as ms
#import flare_mission_sim as fm
import __init__ as fm
import matplotlib.pyplot as plt

flares = ms.load_flare_data_from_hdf5()
flares = ms.shift_by_solar_cycle(flares)

averaging_interval = [datetime.datetime(2021, 12, 1), datetime.datetime(2022, 12, 1)]
dt = averaging_interval[1] - averaging_interval[0]
num_months = int(np.ceil(dt.days / 30.5))

my_flares = flares[averaging_interval[0]:averaging_interval[1]]

class_list = ['X', 'M', 'C']

class_index = {this_class: my_flares['goes_class'].str.contains(this_class)
               for this_class in class_list}

num_per_month = {this_class: my_flares['goes_class'][class_index[this_class]].resample('M').count().values
                 for this_class in class_list}

chance_of_month_with_none = {}
cum_prob_of_no_flares = np.zeros((len(class_list), 12))
for i, this_class in enumerate(class_list):
    x = num_per_month[this_class]
    x_pad = np.zeros(num_months)
    x_pad[0:len(x)] = x
    num_months_with_zero = np.sum(x_pad == 0)
    chance_of_month_with_none.update({this_class: num_months_with_zero / num_months})

    for j in np.arange(1, 13):
        cum_prob_of_no_flares[i, j-1] = binom(j, num_months_with_zero / num_months).pmf(j)

for i, this_class in enumerate(class_list):
    plt.plot(range(1, 13), 1 - cum_prob_of_no_flares[i, :], label=this_class)

plt.legend()
plt.xlabel("Duration of phase e in months")
plt.ylabel('Confidence of getting at least one event')
plt.savefig(fm.in_output_dir('estimate_mission_duration', 'pdf'))
