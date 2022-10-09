import pandas
import numpy as np
import flare_mission_sim as fm
from flare_mission_sim import flare_mission_sim as ms
import matplotlib.pyplot as plt

def trace_flare_stats():
    flare_data = ms.load_flare_data_from_hdf5()

    # remove all C class flares
    flare_data = flare_data[flare_data['goes_class'].str.contains('X|M')]
    # fix bad flare positions by cross-referencing with SFF and manual databases
    flare_data = ms.fix_flare_positions_sff(flare_data)
    flare_data = ms.fix_flare_positions_manual(flare_data)


    trace_flares = pandas.read_csv('/Users/ainglis/physics/mission_simulator/TRACE_M_and_X_flares.csv',index_col=0)

    num1999 = len(trace_flares['1999':'2000'])
    num2000 = len(trace_flares['2000':'2001'])
    num2001 = len(trace_flares['2001':'2002'])
    num2002 = len(trace_flares['2002':'2003'])
    num2003 = len(trace_flares['2003':'2004'])
    num2004 = len(trace_flares['2004':'2005'])
    num2005 = len(trace_flares['2005':'2006'])
    num2006 = len(trace_flares['2006':'2007'])

    numgoes1999 = len(flare_data['1999'])
    numgoes2000 = len(flare_data['2000'])
    numgoes2001 = len(flare_data['2001'])
    numgoes2002 = len(flare_data['2002'])
    numgoes2003 = len(flare_data['2003'])
    numgoes2004 = len(flare_data['2004'])
    numgoes2005 = len(flare_data['2005'])
    numgoes2006 = len(flare_data['2006'])

    frac_observed = [num1999 / numgoes1999, num2000 / numgoes2000, num2001 / numgoes2001, num2002 / numgoes2002,
                         num2003 / numgoes2003, num2004 / numgoes2004, num2005 / numgoes2005]#, num2006 / numgoes2006]


    print(frac_observed)

    plt.figure(1, figsize=(8,3.5))
    plt.subplots_adjust(bottom=0.13)
  #  plt.subplot(1,1,1)
    x = [1999, 2000, 2001, 2002, 2003, 2004, 2005]
    plt.plot(x, frac_observed, 'bo-')
    plt.xlabel('Year')
    plt.ylabel('Fraction of >M1 flares observed')
    plt.title('TRACE observations of flares')
    plt.grid()

 #   plt.subplot(2,1,2)
  #  plt.plot(x, [num1999, num2000, num2001, num2002, num2003, num2004, num2005, num2006],'bo-')
  #  plt.xlabel('Year')
  #  plt.ylabel('Number of flares')
  #  plt.grid()

    plt.savefig('trace_flare_performance.pdf')
    plt.show()
    

    

