
import pandas
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import astropy.units as u
from astropy.coordinates import SkyCoord
from sunpy.coordinates import frames
import datetime

def plot_flare_longitudes():

    plt.figure(1, figsize=(12,4))
    plt.subplots_adjust(bottom=0.1,top=0.9,left=0.05,right=0.95)

    plt.subplot(1,3,1)
    
   # all_flares = pandas.read_csv('/Users/ainglis/physics/mission_simulator/out/new_run_n1000_mcintosh_with_srs_corrections_12hr_noshift/phase_e_flares_foxsi_20200915_1649.csv',index_col=0)
   # flares = pandas.read_csv('/Users/ainglis/physics/mission_simulator/out/new_run_n1000_mcintosh_with_srs_corrections_12hr_noshift/flares_in_view_foxsi_20200915_1649.csv',index_col=0)

    all_flares = pandas.read_csv('/Users/ainglis/physics/mission_simulator/out/n1000_mcintosh_bloomfield_with_srs_corrections_12hr_noshift/phase_e_flares_foxsi_20201210_1624.csv',index_col=0)
    flares = pandas.read_csv('/Users/ainglis/physics/mission_simulator/out/n1000_mcintosh_bloomfield_with_srs_corrections_12hr_noshift/flares_in_view_foxsi_20201210_1624.csv',index_col=0)

    # convert to HG coordinates
    all_flares_lon_hg = []
    all_flares_lat_hg = []

    obs_flares_lon_hg = []
    obs_flares_lat_hg = []

    for i, f in all_flares.iterrows():
        c = SkyCoord(f['hpc_x'] * u.arcsec, f['hpc_y'] * u.arcsec, frame = 'helioprojective', obstime = f['peak_time'], observer='Earth')
        h = c.transform_to(frames.HeliographicStonyhurst)
        all_flares_lon_hg.append(h.lon.value)
        all_flares_lat_hg.append(h.lat.value)

        
    for i, f in flares.iterrows():
        c = SkyCoord(f['hpc_x'] * u.arcsec, f['hpc_y'] * u.arcsec, frame = 'helioprojective', obstime = i, observer='Earth')
        h = c.transform_to(frames.HeliographicStonyhurst)
        obs_flares_lon_hg.append(h.lon.value)
        obs_flares_lat_hg.append(h.lat.value)
        
        
    # now make the plots
    
    
  #  ind = np.where(flares['pointing_hpc_x'] > 0.0)
   # num_west = len(ind[0])
   # ind = np.where(flares['pointing_hpc_x'] < 0.0)
   # num_east = len(ind[0])
   # west_vs_east = num_west / num_east

    ind = np.where(np.array(obs_flares_lon_hg) > 0.0)
    num_west = len(ind[0])
    ind = np.where(np.array(obs_flares_lon_hg) < 0.0)
    num_east = len(ind[0])
    west_vs_east = num_west / num_east
    

    ind = np.where(np.array(all_flares_lon_hg) > 0.0)
    num_allflares_west = len(ind[0])
    ind = np.where(np.array(all_flares_lon_hg) < 0.0)
    num_allflares_east = len(ind[0])
    west_vs_east_allflares = num_allflares_west / num_allflares_east
    print(west_vs_east)
    print(west_vs_east_allflares)
    
    plt.hist(all_flares_lon_hg,bins=np.linspace(-90,90,19),label='All flares')
    plt.hist(obs_flares_lon_hg,bins=np.linspace(-90,90,19), label='Flares observed')
    plt.text(-38,45,'West / East (observed): ' + str(np.round(west_vs_east,3)))
    plt.text(-38,40,'West / East (all): ' + str(np.round(west_vs_east_allflares,3)))
    plt.xlabel('Longitude (degrees)')
    plt.ylabel('N')
    plt.title('McIntosh 12hr delay')
    plt.legend(fontsize=8)


    plt.subplot(1,3,2)
    
   # all_flares = pandas.read_csv('/Users/ainglis/physics/mission_simulator/out/new_run_n1000_hale_with_srs_corrections_12hr_noshift/phase_e_flares_foxsi_20200915_1733.csv',index_col=0)
    all_flares = pandas.read_csv('/Users/ainglis/physics/mission_simulator/out/n1000_mcintosh_bloomfield_with_srs_corrections_12hr_noshift/phase_e_flares_foxsi_20201210_1624.csv',index_col=0)
    flares = pandas.read_csv('/Users/ainglis/physics/mission_simulator/out/new_run_n1000_hale_with_srs_corrections_12hr_noshift/flares_in_view_foxsi_20200915_1733.csv',index_col=0)

     # convert to HG coordinates
    all_flares_lon_hg = []
    all_flares_lat_hg = []

    obs_flares_lon_hg = []
    obs_flares_lat_hg = []

    for i, f in all_flares.iterrows():
        c = SkyCoord(f['hpc_x'] * u.arcsec, f['hpc_y'] * u.arcsec, frame = 'helioprojective', obstime = f['peak_time'], observer='Earth')
        h = c.transform_to(frames.HeliographicStonyhurst)
        all_flares_lon_hg.append(h.lon.value)
        all_flares_lat_hg.append(h.lat.value)

        
    for i, f in flares.iterrows():
        c = SkyCoord(f['pointing_hpc_x'] * u.arcsec, f['pointing_hpc_y'] * u.arcsec, frame = 'helioprojective', obstime = i, observer='Earth')
        h = c.transform_to(frames.HeliographicStonyhurst)
        obs_flares_lon_hg.append(h.lon.value)
        obs_flares_lat_hg.append(h.lat.value)
        
    

    
    ind = np.where(np.array(obs_flares_lon_hg) > 0.0)
    num_west = len(ind[0])
    ind = np.where(np.array(obs_flares_lon_hg) < 0.0)
    num_east = len(ind[0])
    west_vs_east = num_west / num_east
    

    ind = np.where(np.array(all_flares_lon_hg) > 0.0)
    num_allflares_west = len(ind[0])
    ind = np.where(np.array(all_flares_lon_hg) < 0.0)
    num_allflares_east = len(ind[0])
    west_vs_east_allflares = num_allflares_west / num_allflares_east
    print(west_vs_east)
    print(west_vs_east_allflares)
    
    plt.hist(all_flares_lon_hg,bins=np.linspace(-90,90,19),label='All flares')
    plt.hist(obs_flares_lon_hg,bins=np.linspace(-90,90,19), label='Flares observed')
    plt.text(-38,45,'West / East (observed): ' + str(np.round(west_vs_east,3)))
    plt.text(-38,40,'West / East (all): ' + str(np.round(west_vs_east_allflares,3)))
    plt.xlabel('Longitude (degrees)')
    plt.ylabel('N')
    plt.title('Hale + spots 12hr delay')
    plt.legend(fontsize=8)


    plt.subplot(1,3,3)
    
   # all_flares = pandas.read_csv('/Users/ainglis/physics/mission_simulator/out/new_run_n1000_flare_index_with_srs_corrections_12hr_noshift/phase_e_flares_foxsi_20200915_1946.csv',index_col=0)
    all_flares = pandas.read_csv('/Users/ainglis/physics/mission_simulator/out/n1000_mcintosh_bloomfield_with_srs_corrections_12hr_noshift/phase_e_flares_foxsi_20201210_1624.csv',index_col=0)
    flares = pandas.read_csv('/Users/ainglis/physics/mission_simulator/out/new_run_n1000_flare_index_with_srs_corrections_12hr_noshift/flares_in_view_foxsi_20200915_1946.csv',index_col=0)

        # convert to HG coordinates
    all_flares_lon_hg = []
    all_flares_lat_hg = []

    obs_flares_lon_hg = []
    obs_flares_lat_hg = []

    for i, f in all_flares.iterrows():
        c = SkyCoord(f['hpc_x'] * u.arcsec, f['hpc_y'] * u.arcsec, frame = 'helioprojective', obstime = f['peak_time'], observer='Earth')
        h = c.transform_to(frames.HeliographicStonyhurst)
        all_flares_lon_hg.append(h.lon.value)
        all_flares_lat_hg.append(h.lat.value)

        
    for i, f in flares.iterrows():
        c = SkyCoord(f['pointing_hpc_x'] * u.arcsec, f['pointing_hpc_y'] * u.arcsec, frame = 'helioprojective', obstime = i, observer='Earth')
        h = c.transform_to(frames.HeliographicStonyhurst)
        obs_flares_lon_hg.append(h.lon.value)
        obs_flares_lat_hg.append(h.lat.value)
        
    
    ind = np.where(np.array(obs_flares_lon_hg) > 0.0)
    num_west = len(ind[0])
    ind = np.where(np.array(obs_flares_lon_hg) < 0.0)
    num_east = len(ind[0])
    west_vs_east = num_west / num_east
    

    ind = np.where(np.array(all_flares_lon_hg) > 0.0)
    num_allflares_west = len(ind[0])
    ind = np.where(np.array(all_flares_lon_hg) < 0.0)
    num_allflares_east = len(ind[0])
    west_vs_east_allflares = num_allflares_west / num_allflares_east
    print(west_vs_east)
    print(west_vs_east_allflares)
    
    plt.hist(all_flares_lon_hg,bins=np.linspace(-90,90,19),label='All flares')
    plt.hist(obs_flares_lon_hg,bins=np.linspace(-90,90,19), label='Flares observed')
    plt.text(-38,45,'West / East (observed): ' + str(np.round(west_vs_east,3)))
    plt.text(-38,40,'West / East (all): ' + str(np.round(west_vs_east_allflares,3)))
    plt.xlabel('Longitude (degrees)')
    plt.ylabel('N')
    plt.title('Flare index 12hr delay')
    plt.legend(fontsize=8)
    
    
    plt.savefig('flare_longitude_distributions.pdf')


def plot_flare_longitudes_new():
    
    plt.figure(1, figsize=(12,4))
    plt.subplots_adjust(bottom=0.1,top=0.9,left=0.05,right=0.95)

    all_flares = pandas.read_csv('/Users/ainglis/physics/mission_simulator/out/n1000_mcintosh_bloomfield_with_srs_corrections_12hr_noshift/phase_e_flares_foxsi_20201210_1624.csv',index_col=0)


    # convert to HG coordinates
    all_flares_lon_hg = []
    all_flares_lat_hg = []

    obs_flares_lon_hg = []
    obs_flares_lat_hg = []

    for i, f in all_flares.iterrows():
        c = SkyCoord(f['hpc_x'] * u.arcsec, f['hpc_y'] * u.arcsec, frame = 'helioprojective', obstime = f['peak_time'], observer='Earth')
        h = c.transform_to(frames.HeliographicStonyhurst)
        all_flares_lon_hg.append(h.lon.value)
        all_flares_lat_hg.append(h.lat.value)


    t1 = datetime.datetime.now()
        
    for n in range(0,100):
        print(n)
        
        flares = pandas.read_csv('/Users/ainglis/physics/mission_simulator/out/n100_mcintosh_bloomfield_12hr_for_fig8/flares_in_view' + str(n)+'_foxsi_20210224_1357.csv',index_col=0)
        
        for i, f in flares.iterrows():
            c = SkyCoord(f['hpc_x'] * u.arcsec, f['hpc_y'] * u.arcsec, frame = 'helioprojective', obstime = i, observer='Earth')
            h = c.transform_to(frames.HeliographicStonyhurst)
            obs_flares_lon_hg.append(h.lon.value)
            obs_flares_lat_hg.append(h.lat.value)

            
    ind = np.where(np.array(obs_flares_lon_hg) > 0.0)
    num_west = len(ind[0])
    ind = np.where(np.array(obs_flares_lon_hg) < 0.0)
    num_east = len(ind[0])
    west_vs_east = num_west / num_east
    

    ind = np.where(np.array(all_flares_lon_hg) > 0.0)
    num_allflares_west = len(ind[0])
    ind = np.where(np.array(all_flares_lon_hg) < 0.0)
    num_allflares_east = len(ind[0])
    west_vs_east_allflares = num_allflares_west / num_allflares_east
    print(west_vs_east)
    print(west_vs_east_allflares)


    plt.subplot(1,3,1)
    plt.hist(all_flares_lon_hg,bins=np.linspace(-90,90,19),label='All flares')
  #  y, x = np.histogram(obs_flares_lon_hg, bins=np.linspace(-90,90,19))
   # print(x)
   # print(y)
   # y2 = y/10.
   # print(y2)
    
   # plt.hist(y2,bins = x[:-1])
   # plt.hist(y,x, weights = np.zeros_like(y) * 10.0, label='Flares observed')
    
    plt.hist(obs_flares_lon_hg,bins=np.linspace(-90,90,19), weights = np.ones_like(obs_flares_lon_hg) * 0.01, label='Flares observed')
    plt.text(-38,45,'West / East (observed): ' + str(np.round(west_vs_east,3)))
    plt.text(-38,40,'West / East (all): ' + str(np.round(west_vs_east_allflares,3)))
    plt.xlabel('Longitude (degrees)')
    plt.ylabel('N')
    plt.title('McIntosh 12hr delay')
    plt.legend(fontsize=8)


    obs_flares_lon_hg = []
    obs_flares_lat_hg = []

    for n in range(0,100):
        print(n)
        
        flares = pandas.read_csv('/Users/ainglis/physics/mission_simulator/out/n100_hale_12hr_for_fig8/flares_in_view' + str(n)+'_foxsi_20210224_1413.csv',index_col=0)
        
        for i, f in flares.iterrows():
            c = SkyCoord(f['hpc_x'] * u.arcsec, f['hpc_y'] * u.arcsec, frame = 'helioprojective', obstime = i, observer='Earth')
            h = c.transform_to(frames.HeliographicStonyhurst)
            obs_flares_lon_hg.append(h.lon.value)
            obs_flares_lat_hg.append(h.lat.value)

            
    ind = np.where(np.array(obs_flares_lon_hg) > 0.0)
    num_west = len(ind[0])
    ind = np.where(np.array(obs_flares_lon_hg) < 0.0)
    num_east = len(ind[0])
    west_vs_east = num_west / num_east
    print(west_vs_east)
    print(west_vs_east_allflares)
    
    plt.subplot(1,3,2)
    plt.hist(all_flares_lon_hg,bins=np.linspace(-90,90,19),label='All flares')
  #  y, x = np.histogram(obs_flares_lon_hg, bins=np.linspace(-90,90,19))
   # print(x)
   # print(y)
   # y2 = y/10.
   # print(y2)
    
   # plt.hist(y2,bins = x[:-1])
   # plt.hist(y,x, weights = np.zeros_like(y) * 10.0, label='Flares observed')
    
    plt.hist(obs_flares_lon_hg,bins=np.linspace(-90,90,19), weights = np.ones_like(obs_flares_lon_hg) * 0.01, label='Flares observed')
    plt.text(-38,45,'West / East (observed): ' + str(np.round(west_vs_east,3)))
    plt.text(-38,40,'West / East (all): ' + str(np.round(west_vs_east_allflares,3)))
    plt.xlabel('Longitude (degrees)')
    plt.ylabel('N')
    plt.title('Hale 12hr delay')
    plt.legend(fontsize=8)


    obs_flares_lon_hg = []
    obs_flares_lat_hg = []


    for n in range(0,100):
        print(n)
        
        flares = pandas.read_csv('/Users/ainglis/physics/mission_simulator/out/n100_flare_index_12hr_for_fig8/flares_in_view' + str(n)+'_foxsi_20210224_1428.csv',index_col=0)
        
        for i, f in flares.iterrows():
            c = SkyCoord(f['hpc_x'] * u.arcsec, f['hpc_y'] * u.arcsec, frame = 'helioprojective', obstime = i, observer='Earth')
            h = c.transform_to(frames.HeliographicStonyhurst)
            obs_flares_lon_hg.append(h.lon.value)
            obs_flares_lat_hg.append(h.lat.value)

            
    ind = np.where(np.array(obs_flares_lon_hg) > 0.0)
    num_west = len(ind[0])
    ind = np.where(np.array(obs_flares_lon_hg) < 0.0)
    num_east = len(ind[0])
    west_vs_east = num_west / num_east
    print(west_vs_east)
    print(west_vs_east_allflares)
    
    plt.subplot(1,3,3)
    plt.hist(all_flares_lon_hg,bins=np.linspace(-90,90,19),label='All flares')
  #  y, x = np.histogram(obs_flares_lon_hg, bins=np.linspace(-90,90,19))
   # print(x)
   # print(y)
   # y2 = y/10.
   # print(y2)
    
   # plt.hist(y2,bins = x[:-1])
   # plt.hist(y,x, weights = np.zeros_like(y) * 10.0, label='Flares observed')
    
    plt.hist(obs_flares_lon_hg,bins=np.linspace(-90,90,19), weights = np.ones_like(obs_flares_lon_hg) * 0.01, label='Flares observed')
    plt.text(-38,45,'West / East (observed): ' + str(np.round(west_vs_east,3)))
    plt.text(-38,40,'West / East (all): ' + str(np.round(west_vs_east_allflares,3)))
    plt.xlabel('Longitude (degrees)')
    plt.ylabel('N')
    plt.title('Flare index 12hr delay')
    plt.legend(fontsize=8)


    
    t2 = datetime.datetime.now()

    total_time = t2 - t1
    print(total_time)

    plt.savefig('flare_longitude_distributions_new.pdf')

    
    plt.show()        
    
        

