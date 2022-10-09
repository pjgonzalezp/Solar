import pandas as pd
import os
import numpy as np
import matplotlib.pyplot as plt
import flare_mission_sim as fm
from flare_mission_sim import util
import copy

def mcintosh_rank_order():

    bornmann_shaw_productivities =  pd.read_csv(os.path.join(fm.data_dir, 'mcintosh_flare_rates_bornmann_shaw_1994.csv'), index_col=0)
    bornmann_shaw_productivities.columns = bornmann_shaw_productivities.columns.str.strip()
    #rename column to avoid confusion
    bornmann_shaw_productivities.rename(columns = {'flare_productivity':'bornmann_flare_productivity'}, inplace=True)

    bloomfield_productivities = pd.read_csv(os.path.join(fm.data_dir, 'mcintosh_flare_rates_bloomfield_2012.csv'), index_col=0)
    bloomfield_productivities.columns = bloomfield_productivities.columns.str.strip()
    

    #merge the two dataframes into one

    combined_productivities = pd.concat([bloomfield_productivities, bornmann_shaw_productivities], axis = 1)
    combined_productivities.sort_values(by = 'flare_productivity', ascending = False, inplace=True)

    prod = []
    for ind, c in combined_productivities.iterrows():
        prod.append(util.mcintosh_class_productivity(ind))
    combined_productivities['average_productivities'] = prod
    
    
    plt.figure(1, figsize = (10,8))
    plt.subplots_adjust(left = 0.05, right = 0.95, bottom = 0.06, top = 0.95)
    plt.plot(combined_productivities['flare_productivity'], combined_productivities.index,'o-', label='Bloomfield et al. 2012')
    plt.plot(combined_productivities['bornmann_flare_productivity'], combined_productivities.index,'o-',label = 'Bornmann and Shaw 1994')
   # plt.plot(combined_productivities['average_productivities'], combined_productivities.index, 'o-',label = 'Averaging method')
    plt.xlabel('Flare productivity (flares /day)', fontsize=14)
    plt.legend()
    ax = plt.gca()
    ax.set_xticks(np.arange(0,2.75,0.25))
    plt.grid()
    plt.savefig('mcintosh_productivity_comparison_values.pdf')
    
   # plt.show()

    
    bloomfield_rank_order = bloomfield_productivities.rank(ascending = False)
    bloomfield_rank_order.sort_values(by='flare_productivity', inplace=True)

    bornmann_rank_order = bornmann_shaw_productivities.rank(ascending = False)
    bornmann_rank_order.sort_values(by='bornmann_flare_productivity', inplace=True)

    averaging_rank_order = combined_productivities['average_productivities'].rank(ascending = False)
    averaging_rank_order.sort_values(inplace = True)

    

    plt.figure(2,figsize=(10,8))
    plt.subplots_adjust(left = 0.05, right = 0.95, bottom = 0.06, top = 0.95)
    plt.plot(bloomfield_rank_order['flare_productivity'],bloomfield_rank_order.index, 'o-', label='Bloomfield et al. 2012')
    plt.plot(bornmann_rank_order['bornmann_flare_productivity'],bornmann_rank_order.index, 'o-',label = 'Bornmann and Shaw 1994')
  #  plt.plot(averaging_rank_order.values, averaging_rank_order.index, 'o-', label = 'Averaging method')
    #plt.plot(['EKC','AXX'],[1,53],'--')
    plt.xlabel('Flare productivity rank (1 = most productive)', fontsize=14)
    plt.xlim([55,0])
    plt.legend()
    ax = plt.gca()
    ax.set_xticks(np.arange(0,55,5))
    plt.grid()

    plt.savefig('mcintosh_productivity_comparison_rank_order.pdf')
   # plt.show()

