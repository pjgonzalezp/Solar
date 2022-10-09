import numpy as np
import matplotlib.pyplot as plt
import pandas
from scipy import io

import flare_mission_sim as fm
from flare_mission_sim import flare_mission_sim as ms
from flare_mission_sim import util

import seaborn as sns
sns.set_context('paper')
sns.set_style('ticks')

def plot_solar_cycle():

    f107_data = io.readsav('f107_AA.sav')
    f107_time = f107_data['f107_aa_t_years']
    f107_flux = f107_data['f107_aa_val_sfu']
    
    flares = ms.load_flare_data_from_hdf5()
    flares2 = flares[flares['goes_class'].str.contains('M|X')]
    flares_x = flares[flares['goes_class'].str.contains('X')]
    flares_m = flares[flares['goes_class'].str.contains('M')]

    years = np.linspace(1996,2017,22, dtype='int')

    num_xflares = []
    num_mflares = []
    
    for y in years:
        y2 = str(y)
        num_xflares.append(len(flares_x[y2]))
        num_mflares.append(len(flares_m[y2]))
        

    plt.figure(1,figsize=(12,5))
    plt.tick_params(direction='in', labelsize=14)

    plt.axvspan(2011.08, 2015.0, alpha=0.5, color = 'powderblue')
    plt.text(2011.2, 300,'Mission simulation',fontsize=14)

    plt.text(1997.4, 170, 'Cycle 23', fontsize=12)
    plt.text(2009.4, 100, 'Cycle 24', fontsize=12)
    
    plt.bar(years + 0.5,num_mflares, 0.98, label = 'M flares')
    plt.bar(years + 0.5,num_xflares, 0.98, bottom = num_mflares, label = 'X flares')
    plt.ylabel('Number of flares', fontsize=14)
    plt.xlabel('Year', fontsize=14)
    plt.xlim([1997,2018])

    plt.locator_params(axis='x', nbins=11)

    
    plt.legend(fontsize=14)
    
    ax = plt.gca()
    ax2 = ax.twinx()

    sunspots = pandas.read_csv('/Users/ainglis/physics/mission_simulator/SN_y_tot_V2.0.csv', sep = ';', index_col=0, header=None)
    ax2.plot(sunspots.index, sunspots[1],linewidth=2, color='red')
    ax2.yaxis.label.set_color('red')
    plt.ylabel('Average monthly Sunspot number', color ='red', fontsize=14)
    ax2.yaxis.label.set_color('red')
    ax2.tick_params(axis='y', colors='red', labelsize=14)

    plt.savefig('mission_simulator_solar_cycle.pdf')
    
    plt.show()

    

