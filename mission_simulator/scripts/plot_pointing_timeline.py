import pandas
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import seaborn
import datetime
from flare_mission_sim import flare_mission_sim as ms
import numpy as np

def plot_pointing_timeline(pointing_command_file):

    pointings = pandas.read_csv(pointing_command_file,index_col = 0)

    seg1 = ['2011-02-01','2011-05-01']

    segments = seg1

    xformat = mdates.DateFormatter('%Y%M%D')
    fig = plt.figure(1,figsize=(16,6))

    pointings_to_plot = pointings[segments[0]:segments[1]]

    times = []
    for t in pointings_to_plot.index.values:
        times.append(datetime.datetime.strptime(t,'%Y-%m-%d %H:%M:%S.%f'))
        
    
    plt.plot(times,pointings_to_plot['noaa'],drawstyle='steps-pre')
    plt.ylabel('NOAA AR')
    plt.grid()
    ax = plt.gca()
     #   ax.xaxis_date()
      #  ax.xaxis.set_major_formatter(xformat)
    fig.autofmt_xdate()
 #       ax.xaxis.set_major_locator(plt.MaxNLocator(6))

    plt.show()

    

def compare_pointing_timelines():

  #  mcintosh_file = '/Users/ainglis/physics/mission_simulator/mcintosh_pointing_commands.csv'
  #  mcintosh_file = '/Users/ainglis/physics/mission_simulator/out/number_of_pointings/list_of_pointing_commands_mcintosh.csv'
    mcintosh_file = '/Users/ainglis/physics/mission_simulator/out/number_of_pointings/list_of_pointing_commands_mcintosh_bloomfield.csv'
    mcintosh_pointings = pandas.read_csv(mcintosh_file,index_col = 0)

  #  hale_file = '/Users/ainglis/physics/mission_simulator/hale_pointing_commands.csv'
    hale_file = '/Users/ainglis/physics/mission_simulator/out/number_of_pointings/list_of_pointing_commands_hale.csv'
    hale_pointings = pandas.read_csv(hale_file,index_col = 0)

    #random_file = '/Users/ainglis/physics/mission_simulator/random_pointing_commands.csv'
    random_file = '/Users/ainglis/physics/mission_simulator/out/number_of_pointings/list_of_pointing_commands_random.csv'
    random_pointings = pandas.read_csv(random_file,index_col = 0)

    flare_index_file = '/Users/ainglis/physics/mission_simulator/out/number_of_pointings/list_of_pointing_commands_flare_index.csv'
    flare_index_pointings = pandas.read_csv(flare_index_file, index_col = 0)

    seg1 = ['2011-12-24','2012-03-05']
    seg2 = ['2012-01-01','2012-03-01']

    segments = seg1

    xformat = mdates.DateFormatter('%Y%M%D')
    fig = plt.figure(1,figsize=(16,12))
    plt.subplots_adjust(bottom=0.12, top = 0.95, left = 0.1, right = 0.95)

    mcintosh_pointings_to_plot = mcintosh_pointings[segments[0]:segments[1]]
    hale_pointings_to_plot = hale_pointings[segments[0]:segments[1]]
    random_pointings_to_plot = random_pointings[segments[0]:segments[1]]
    flare_index_pointings_to_plot = flare_index_pointings[segments[0]:segments[1]]

    num_mcintosh = len(mcintosh_pointings_to_plot[seg2[0]:seg2[1]])
    num_hale = len(hale_pointings_to_plot[seg2[0]:seg2[1]])
    num_random = len(random_pointings_to_plot[seg2[0]:seg2[1]])
    num_flare_index = len(flare_index_pointings_to_plot[seg2[0]:seg2[1]])
    
    mcintosh_times = []
    for t in mcintosh_pointings_to_plot.index.values:
        mcintosh_times.append(datetime.datetime.strptime(t,'%Y-%m-%d %H:%M:%S.%f'))

    hale_times = []
    for t in hale_pointings_to_plot.index.values:
        hale_times.append(datetime.datetime.strptime(t,'%Y-%m-%d %H:%M:%S.%f'))

    random_times = []
    for t in random_pointings_to_plot.index.values:
        random_times.append(datetime.datetime.strptime(t,'%Y-%m-%d %H:%M:%S.%f'))

    flare_index_times = []
    for t in flare_index_pointings_to_plot.index.values:
        flare_index_times.append(datetime.datetime.strptime(t,'%Y-%m-%d %H:%M:%S.%f'))



    # Returns a Pandas dataframe.  Note that the index is the start time of the flare.
    flare_data = ms.load_flare_data_from_hdf5()

    # remove all C class flares
    flare_data = flare_data[flare_data['goes_class'].str.contains('X|M')]
    # fix bad flare positions by cross-referencing with SFF and manual databases
    flare_data = ms.fix_flare_positions_sff(flare_data)
    flare_data = ms.fix_flare_positions_manual(flare_data)

    relevant_flares = flare_data['2012']
    relevant_flares = relevant_flares[relevant_flares['noaa'] != 0]
    relevant_flares = relevant_flares[seg2[0]:seg2[1]]
    

    plt.subplot(4,1,1)
    plt.plot(mcintosh_times,mcintosh_pointings_to_plot['ar_target'],drawstyle='steps-post',label='mcintosh',
                 color='blue')
    plt.plot(relevant_flares['peak_time'].values, relevant_flares['noaa'].values,'.',color='black', markersize=10)
    plt.gca().yaxis.set_major_locator(plt.MultipleLocator(10))
    plt.grid(which='both', linestyle='dashed')
    plt.ylabel('NOAA AR', fontsize=12)
    
    plt.legend(fontsize=12)
    plt.tick_params(labelsize=12)
    plt.xlim([datetime.datetime(2012,1,1),datetime.datetime(2012,3,1)])
    plt.text(datetime.datetime(2012,2,25),11410,'N = ' + str(num_mcintosh))

    plt.subplot(4,1,2)
    plt.plot(hale_times,hale_pointings_to_plot['ar_target'],drawstyle='steps-post',label='hale',
                 color='green')
    plt.plot(relevant_flares['peak_time'].values, relevant_flares['noaa'].values,'.',color='black', markersize=10)
    plt.gca().yaxis.set_major_locator(plt.MultipleLocator(10))
    plt.grid(linestyle='dashed')
    plt.ylabel('NOAA AR', fontsize=12)
    plt.legend(loc=2, fontsize=12)
    plt.tick_params(labelsize=12)
    plt.xlim([datetime.datetime(2012,1,1),datetime.datetime(2012,3,1)])
    plt.text(datetime.datetime(2012,2,25),11410,'N = ' + str(num_hale))

    plt.subplot(4,1,3)
    plt.plot(flare_index_times,flare_index_pointings_to_plot['ar_target'],drawstyle='steps-post',label='flare index',
                 color='purple')
    plt.gca().yaxis.set_major_locator(plt.MultipleLocator(10))
    plt.ylabel('NOAA AR', fontsize=12)
    plt.plot(relevant_flares['peak_time'].values, relevant_flares['noaa'].values,'.',color='black', markersize=10)
    plt.grid(linestyle='dashed')
    plt.xlim([datetime.datetime(2012,1,1),datetime.datetime(2012,3,1)])
    plt.text(datetime.datetime(2012,2,25),11410,'N = ' + str(num_flare_index))
    plt.legend(fontsize=12)
    plt.tick_params(labelsize=12)

    plt.subplot(4,1,4)
    plt.plot(random_times,random_pointings_to_plot['ar_target'],drawstyle='steps-post',label='random',
                 color='red')
    plt.gca().yaxis.set_major_locator(plt.MultipleLocator(10))
    plt.ylabel('NOAA AR', fontsize=12)
    plt.plot(relevant_flares['peak_time'].values, relevant_flares['noaa'].values,'.',color='black', markersize=10)
    plt.grid(linestyle='dashed')
    plt.xlim([datetime.datetime(2012,1,1),datetime.datetime(2012,3,1)])
    plt.text(datetime.datetime(2012,2,25),11410,'N = ' + str(num_random))
    plt.legend(fontsize=12)
    plt.tick_params(labelsize=12)
    
    ax = plt.gca()
     #   ax.xaxis_date()
      #  ax.xaxis.set_major_formatter(xformat)
    fig.autofmt_xdate()
 #       ax.xaxis.set_major_locator(plt.MaxNLocator(6))


    plt.savefig('pointing_timeline_4panel.pdf')
    plt.show()

    
    xformat = mdates.DateFormatter('%Y%M%D')
    fig2 = plt.figure(2,figsize=(16,6))
 
    plt.plot(mcintosh_times,mcintosh_pointings_to_plot['ar_target'],drawstyle='steps-pre',label='mcintosh',
                 color='blue')
    plt.plot(hale_times,hale_pointings_to_plot['ar_target'],drawstyle='steps-pre',label='hale',
                 color='green')
    plt.plot(flare_index_times,flare_index_pointings_to_plot['ar_target'],drawstyle='steps-pre',label='flare index',
                 color='purple')
    plt.plot(random_times,random_pointings_to_plot['ar_target'],drawstyle='steps-pre',label='random',
                 color='red')
    plt.xlim([datetime.datetime(2012,1,1),datetime.datetime(2012,3,1)])
    plt.gca().yaxis.set_major_locator(plt.MultipleLocator(10))
    plt.grid(which='both')
    plt.ylabel('NOAA AR')

    plt.legend(loc=2)

    plt.text(datetime.datetime(2012,2,25),11410,'N = ' + str(num_mcintosh),color='blue')
    plt.text(datetime.datetime(2012,2,25),11405,'N = ' + str(num_hale),color='green')
    plt.text(datetime.datetime(2012,2,25),11400,'N = ' + str(num_flare_index),color='purple')
    plt.text(datetime.datetime(2012,2,25),11395,'N = ' + str(num_random),color='red')
    
    plt.show()

    

def differences_in_target():
    #mcintosh_targets = pandas.read_csv('/Users/ainglis/physics/mission_simulator/out/number_of_pointings/list_of_ar_targets_mcintosh.csv', index_col=0)
    mcintosh_targets = pandas.read_csv('/Users/ainglis/physics/mission_simulator/out/number_of_pointings/list_of_ar_targets_mcintosh_bloomfield.csv', index_col=0)
    hale_targets = pandas.read_csv('/Users/ainglis/physics/mission_simulator/out/number_of_pointings/list_of_ar_targets_hale.csv',index_col=0)
    flare_index_targets = pandas.read_csv('/Users/ainglis/physics/mission_simulator/out/number_of_pointings/list_of_ar_targets_flare_index.csv',index_col=0)
    random_targets = pandas.read_csv('/Users/ainglis/physics/mission_simulator/out/number_of_pointings/list_of_ar_targets_random.csv',index_col=0)

    mcintosh_vs_hale = mcintosh_targets['ar_target'] == hale_targets['ar_target']
    mcintosh_vs_hale_same = np.sum(mcintosh_vs_hale)
    mcintosh_vs_hale_different = len(mcintosh_vs_hale) - mcintosh_vs_hale_same

    mcintosh_vs_flare_index = mcintosh_targets['ar_target'] == flare_index_targets['ar_target']
    mcintosh_vs_flare_index_same = np.sum(mcintosh_vs_flare_index)
    mcintosh_vs_flare_index_different = len(mcintosh_vs_flare_index) - mcintosh_vs_flare_index_same

    mcintosh_vs_random = mcintosh_targets['ar_target'] == random_targets['ar_target']
    mcintosh_vs_random_same = np.sum(mcintosh_vs_random)
    mcintosh_vs_random_different = len(mcintosh_vs_random) - mcintosh_vs_random_same

    hale_vs_flare_index = hale_targets['ar_target'] == flare_index_targets['ar_target']
    hale_vs_flare_index_same = np.sum(hale_vs_flare_index)
    hale_vs_flare_index_different = len(hale_vs_flare_index) - hale_vs_flare_index_same

    hale_vs_random = hale_targets['ar_target'] == random_targets['ar_target']
    hale_vs_random_same = np.sum(hale_vs_random)
    hale_vs_random_different = len(hale_vs_random) - hale_vs_random_same

    flare_index_vs_random = flare_index_targets['ar_target'] == random_targets['ar_target']
    flare_index_vs_random_same = np.sum(flare_index_vs_random)
    flare_index_vs_random_different = len(flare_index_vs_random) - flare_index_vs_random_same

    num_targets_each_day = []
    for i, t in enumerate(mcintosh_targets.index):
        num = np.unique((mcintosh_targets.iloc[i]['ar_target'],hale_targets.iloc[i]['ar_target'],
                           flare_index_targets.iloc[i]['ar_target'], random_targets.iloc[i]['ar_target']))
        num_targets_each_day.append(len(num))

    print(num_targets_each_day)

    plt.figure(1, figsize=(12,6))

    x = np.linspace(1,6,6)
    width = 0.35
    tick_labels = ['Mcintosh\n vs Hale', 'Mcintosh\n vs Flare Index',
                       'Hale vs\n Flare Index','Mcintosh\n vs Random', 'Hale vs\n Random', 'Flare index\n vs Random']
    ysame = [mcintosh_vs_hale_same, mcintosh_vs_flare_index_same,
                 hale_vs_flare_index_same, mcintosh_vs_random_same, hale_vs_random_same, flare_index_vs_random_same]
    ydiff = [mcintosh_vs_hale_different, mcintosh_vs_flare_index_different,
                 hale_vs_flare_index_different, mcintosh_vs_random_different, hale_vs_random_different, flare_index_vs_random_different]

    percent_same =  (np.array(ysame) / (np.array(ysame) + np.array(ydiff))) * 100.
    percent_diff = 100.0 - percent_same
    

    plt.bar(x - width/2, ysame, width, tick_label = tick_labels, label = 'Same target')
    plt.bar(x + width/2, ydiff, width, label = 'Different target')
    plt.tick_params(labelsize=12)
    for i, p in enumerate(percent_same):
        plt.text(x[i] - width + 0.05, ysame[i] - 50, str(np.round(p,1)) + '%', color='white')

    for i, p in enumerate(percent_diff):
        plt.text(x[i] + 0.05, ydiff[i] - 50, str(np.round(p,1)) + '%', color='black')
        
    plt.ylabel('Number of days')
    plt.legend()

    plt.savefig('pointing_differences.pdf')
    plt.show()

    
    plt.figure(2, figsize=(6,6))
    h = np.bincount(num_targets_each_day)

  
    plt.bar(np.linspace(1,4,4), h[1:])
    plt.ylabel('Number of days')
    plt.xlabel('Number of AR targets')
    plt.tick_params(labelsize=12)

    plt.show()
    
    
    
    

    
    
    
    


    

    
