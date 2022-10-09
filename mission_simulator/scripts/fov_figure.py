import seaborn as sns
import numpy as np
import os
import matplotlib.pyplot as plt
import pandas

def plot_fov_figure():

    base_path = '/Users/ainglis/physics/mission_simulator/out'
   # files = ['n100_mcintosh_with_srs_corrections_12hr_noshift_1arcmin_fov/mflare_observing_fractions_foxsi_20200928_1546.txt',
    #             'n100_mcintosh_with_srs_corrections_12hr_noshift_2arcmin_fov/mflare_observing_fractions_foxsi_20200928_1602.txt',
    #             'n100_mcintosh_with_srs_corrections_12hr_noshfit_3arcmin_fov/mflare_observing_fractions_foxsi_20200928_1459.txt',
    #             'n100_mcintosh_with_srs_corrections_12hr_noshift_5arcmin_fov/mflare_observing_fractions_foxsi_20200928_1515.txt',
    #             'n100_mcintosh_with_srs_corrections_12hr_noshift_7arcmin_fov/mflare_observing_fractions_foxsi_20200928_1530.txt',
    #             'n100_mcintosh_with_srs_corrections_12hr_noshift_9arcmin_fov/mflare_observing_fractions_foxsi_20200928_1626.txt',
    #             'n100_mcintosh_with_srs_corrections_12hr_noshift_11arcmin_fov/mflare_observing_fractions_foxsi_20200928_1640.txt',
    #             'n100_mcintosh_with_srs_corrections_12hr_noshift_13arcmin_fov/mflare_observing_fractions_foxsi_20200928_1654.txt',
    #             'n100_mcintosh_with_srs_corrections_12hr_noshift_15arcmin_fov/mflare_observing_fractions_foxsi_20200928_1707.txt',
    #             'n100_mcintosh_with_srs_corrections_12hr_noshfit_17.5arcmin_fov/mflare_observing_fractions_foxsi_20200929_1059.txt',
    #             'n100_mcintosh_with_srs_corrections_12hr_noshift_20arcmin_fov/mflare_observing_fractions_foxsi_20200928_1721.txt',
    #             'n100_mcintosh_with_srs_corrections_12hr_noshift_25arcmin_fov/mflare_observing_fractions_foxsi_20200929_1035.txt',
    #             'n100_mcintosh_with_srs_corrections_12hr_noshift_30arcmin_fov/mflare_observing_fractions_foxsi_20200929_1610.txt',
    #             'n100_mcintosh_with_srs_corrections_12hr_noshfit_35arcmin_fov/mflare_observing_fractions_foxsi_20200929_1625.txt']

    files = ['n100_mcintosh_bloomfield_with_srs_corrections_12hr_noshift_1arcmin_fov/mflare_observing_totals_foxsi_20210107_1231.txt',
                 'n100_mcintosh_bloomfield_with_srs_corrections_12hr_noshift_2arcmin_fov/mflare_observing_totals_foxsi_20210107_1249.txt',
                 'n100_mcintosh_bloomfield_with_srs_corrections_12hr_noshift_3arcmin_fov/mflare_observing_totals_foxsi_20210106_1437.txt',
                 'n100_mcintosh_bloomfield_with_srs_corrections_12hr_noshift_5arcmin_fov/mflare_observing_totals_foxsi_20210106_1454.txt',
                 'n100_mcintosh_bloomfield_with_srs_corrections_12hr_noshift_7arcmin_fov/mflare_observing_totals_foxsi_20210106_1516.txt',
                 'n100_mcintosh_bloomfield_with_srs_corrections_12hr_noshift_9arcmin_fov/mflare_observing_totals_foxsi_20210106_1531.txt',
                 'n100_mcintosh_bloomfield_with_srs_corrections_12hr_noshift_11arcmin_fov/mflare_observing_totals_foxsi_20210106_1548.txt',
                 'n100_mcintosh_bloomfield_with_srs_corrections_12hr_noshift_13arcmin_fov/mflare_observing_totals_foxsi_20210106_1604.txt',
                 'n100_mcintosh_bloomfield_with_srs_corrections_12hr_noshift_15arcmin_fov/mflare_observing_totals_foxsi_20210106_1620.txt',
                 'n100_mcintosh_bloomfield_with_srs_corrections_12hr_noshift_17.5arcmin_fov/mflare_observing_totals_foxsi_20210106_1638.txt',
                 'n100_mcintosh_bloomfield_with_srs_corrections_12hr_noshift_20arcmin_fov/mflare_observing_totals_foxsi_20210106_1652.txt',
                 'n100_mcintosh_bloomfield_with_srs_corrections_12hr_noshift_25arcmin_fov/mflare_observing_totals_foxsi_20210106_1705.txt',
                 'n100_mcintosh_bloomfield_with_srs_corrections_12hr_noshift_30arcmin_fov/mflare_observing_totals_foxsi_20210106_1746.txt',
                 'n100_mcintosh_bloomfield_with_srs_corrections_12hr_noshift_35arcmin_fov/mflare_observing_totals_foxsi_20210107_1047.txt']

    xfiles = ['n100_mcintosh_bloomfield_with_srs_corrections_12hr_noshift_1arcmin_fov/xflare_observing_totals_foxsi_20210107_1231.txt',
                 'n100_mcintosh_bloomfield_with_srs_corrections_12hr_noshift_2arcmin_fov/xflare_observing_totals_foxsi_20210107_1249.txt',
                 'n100_mcintosh_bloomfield_with_srs_corrections_12hr_noshift_3arcmin_fov/xflare_observing_totals_foxsi_20210106_1437.txt',
                 'n100_mcintosh_bloomfield_with_srs_corrections_12hr_noshift_5arcmin_fov/xflare_observing_totals_foxsi_20210106_1454.txt',
                 'n100_mcintosh_bloomfield_with_srs_corrections_12hr_noshift_7arcmin_fov/xflare_observing_totals_foxsi_20210106_1516.txt',
                 'n100_mcintosh_bloomfield_with_srs_corrections_12hr_noshift_9arcmin_fov/xflare_observing_totals_foxsi_20210106_1531.txt',
                 'n100_mcintosh_bloomfield_with_srs_corrections_12hr_noshift_11arcmin_fov/xflare_observing_totals_foxsi_20210106_1548.txt',
                 'n100_mcintosh_bloomfield_with_srs_corrections_12hr_noshift_13arcmin_fov/xflare_observing_totals_foxsi_20210106_1604.txt',
                 'n100_mcintosh_bloomfield_with_srs_corrections_12hr_noshift_15arcmin_fov/xflare_observing_totals_foxsi_20210106_1620.txt',
                 'n100_mcintosh_bloomfield_with_srs_corrections_12hr_noshift_17.5arcmin_fov/xflare_observing_totals_foxsi_20210106_1638.txt',
                 'n100_mcintosh_bloomfield_with_srs_corrections_12hr_noshift_20arcmin_fov/xflare_observing_totals_foxsi_20210106_1652.txt',
                 'n100_mcintosh_bloomfield_with_srs_corrections_12hr_noshift_25arcmin_fov/xflare_observing_totals_foxsi_20210106_1705.txt',
                 'n100_mcintosh_bloomfield_with_srs_corrections_12hr_noshift_30arcmin_fov/xflare_observing_totals_foxsi_20210106_1746.txt',
                 'n100_mcintosh_bloomfield_with_srs_corrections_12hr_noshift_35arcmin_fov/xflare_observing_totals_foxsi_20210107_1047.txt']



    fov_diameters = [1,2,3,5,7,9,11,13,15,17.5,20,25,30,35]

    means = []
    standard_devs = []
    for i, f in enumerate(files):
        mresult = np.loadtxt(os.path.join(base_path,f))
        xresult = np.loadtxt(os.path.join(base_path,xfiles[i]))
        result = (mresult + xresult) / 589.
                                 
        m = np.mean(result)
        std = np.std(result)

        means.append(m)
        standard_devs.append(std)


    plt.figure(1)
    plt.plot(fov_diameters, means, 'bo')
    plt.errorbar(fov_diameters, means, yerr = standard_devs)
    plt.xlabel('FOV side length (arcminutes)')
    plt.ylabel('Fraction of >M1 flares observed')
    plt.axvline(5.0, label='Typical AR size scale', linestyle='dashed',color='black')
    plt.axvline(31.66, label='Solar diameter', linestyle='dashdot', color='black')
    plt.axvline(9.8, label='Baseline simulation', linestyle='dashed', color = 'red')
    plt.legend()
    plt.title('Mcintosh 12hr delay')
    plt.grid()
    plt.savefig('fov_figure.pdf')
    
    

def flares_from_other_ar():

    pointing_commands = pandas.read_csv('/Users/ainglis/physics/mission_simulator/out/number_of_pointings/list_of_pointing_commands_mcintosh.csv',index_col=0)

    base_path = '/Users/ainglis/physics/mission_simulator/out'
    flare_files = ['n100_mcintosh_bloomfield_with_srs_corrections_12hr_noshift_1arcmin_fov/flares_in_view_foxsi_20210107_1231.csv',
                 'n100_mcintosh_bloomfield_with_srs_corrections_12hr_noshift_2arcmin_fov/flares_in_view_foxsi_20210107_1249.csv',
                 'n100_mcintosh_bloomfield_with_srs_corrections_12hr_noshift_3arcmin_fov/flares_in_view_foxsi_20210106_1437.csv',
                 'n100_mcintosh_bloomfield_with_srs_corrections_12hr_noshift_5arcmin_fov/flares_in_view_foxsi_20210106_1454.csv',
                 'n100_mcintosh_bloomfield_with_srs_corrections_12hr_noshift_7arcmin_fov/flares_in_view_foxsi_20210106_1516.csv',
                 'n100_mcintosh_bloomfield_with_srs_corrections_12hr_noshift_9arcmin_fov/flares_in_view_foxsi_20210106_1531.csv',
                 'n100_mcintosh_bloomfield_with_srs_corrections_12hr_noshift_11arcmin_fov/flares_in_view_foxsi_20210106_1548.csv',
                 'n100_mcintosh_bloomfield_with_srs_corrections_12hr_noshift_13arcmin_fov/flares_in_view_foxsi_20210106_1604.csv',
                 'n100_mcintosh_bloomfield_with_srs_corrections_12hr_noshift_15arcmin_fov/flares_in_view_foxsi_20210106_1620.csv',
                 'n100_mcintosh_bloomfield_with_srs_corrections_12hr_noshift_17.5arcmin_fov/flares_in_view_foxsi_20210106_1638.csv',
                 'n100_mcintosh_bloomfield_with_srs_corrections_12hr_noshift_20arcmin_fov/flares_in_view_foxsi_20210106_1652.csv',
                 'n100_mcintosh_bloomfield_with_srs_corrections_12hr_noshift_25arcmin_fov/flares_in_view_foxsi_20210106_1705.csv',
                 'n100_mcintosh_bloomfield_with_srs_corrections_12hr_noshift_30arcmin_fov/flares_in_view_foxsi_20210106_1746.csv',
                 'n100_mcintosh_bloomfield_with_srs_corrections_12hr_noshift_35arcmin_fov/flares_in_view_foxsi_20210107_1047.csv']

    allflare_file = '/Users/ainglis/physics/mission_simulator/out/n100_mcintosh_bloomfield_with_srs_corrections_12hr_noshift_1arcmin_fov/phase_e_flares_foxsi_20210107_1231.csv'
    all_flares = pandas.read_csv(allflare_file, index_col=0)          

    same_target_total = []
    different_target_total = []
    same_target_fraction = []
    different_target_fraction = []
    fov_diameters = [1,2,3,5,7,9,11,13,15,17.5,20,25,30,35]
    
    for f in flare_files:
        
        flares = pandas.read_csv(os.path.join(base_path,f), index_col=0)

        flare_from_same_target = 0
        flare_from_different_target = 0
        
        #get the source ARs for each flare observed
        for index, fl in flares.iterrows():
            l = all_flares.loc[index]
            loc = np.searchsorted(pointing_commands.index, index) - 1
            target_ar = pointing_commands.iloc[loc]['ar_target']

            if l['noaa'] == target_ar:
                flare_from_same_target += 1
            else:
                flare_from_different_target += 1

        same_target_total.append(flare_from_same_target)
        different_target_total.append(flare_from_different_target)

        same_target_fraction.append(flare_from_same_target / len(flares))
        different_target_fraction.append(flare_from_different_target / len(flares))

    print(same_target_total)
    print(different_target_total)
    print(same_target_fraction)
    print(different_target_fraction)
    
    print(len(all_flares))

    plt.figure(1)
    plt.plot(fov_diameters, same_target_fraction, 'ro')#, label='flares from mission AR target')
   # plt.plot(fov_diameters, different_target_fraction, 'ro', label ='flares from different AR')
    plt.xlabel('FOV side length (arcminutes)')
    plt.ylabel('Fraction of flares observed from current AR target')
    plt.axvline(5.0, linestyle='dashed',color='black')
    plt.axvline(31.66, linestyle='dashdot', color='black')
    plt.axvline(9.8, linestyle='dashed',color='red')
   # plt.legend(loc=1)
    plt.title('Mcintosh 12hr delay')
    plt.grid()
    plt.savefig('fov_flares_from_ar_target.pdf')



    



        

            

            
        
        
    
    

