import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

sns.set_context('paper')
sns.set_style('ticks')

def hist_to_probability_weight(d):
    """
    :param d:
    :return:

    Reference
    ---------
    This link
    https: // github.com / matplotlib / matplotlib / issues / 10398  # issuecomment-366021979
    explains how to plot a probability distribution where each of the bars indicates the total
    probability in the bin.  This function returns the weighting needed to do this.
    """
    return np.ones_like(d)/len(d)


def regular_hist_bins(d, bin_width):
    """
    Create same-sized bins that span the values in the input data.

    :param d:
    :param bin_width:
    :return:
    """
    lower = bin_width*(min(d) // bin_width)
    upper = bin_width*(1 + (max(d) // bin_width))
    return np.arange(lower, upper, bin_width)


def plot_performance_figure():

  #  mcintosh_fractions_file_24hr = '/Users/ainglis/physics/mission_simulator/out/new_run_n1000_mcintosh_with_srs_corrections_24hr_noshift/mflare_observing_fractions_foxsi_20200916_1139.txt'
  #  hale_fractions_file_24hr = '/Users/ainglis/physics/mission_simulator/out/new_run_n1000_hale_with_srs_corrections_24hr_noshift/mflare_observing_fractions_foxsi_20200916_1204.txt'
  #  flare_index_fractions_file_24hr = '/Users/ainglis/physics/mission_simulator/out/new_run_n1000_flare_index_with_srs_corrections_24hr_noshift/mflare_observing_fractions_foxsi_20200916_1346.txt'
   # random_fractions_file_24hr = '/Users/ainglis/physics/mission_simulator/out/n1000_random_with_srs_corrections_24hr_noshift_random_fixed/mflare_observing_fractions_foxsi_20201008_1949.txt'
    
   # mcintosh_fractions_file_12hr = '/Users/ainglis/physics/mission_simulator/out/new_run_n1000_mcintosh_with_srs_corrections_12hr_noshift/mflare_observing_fractions_foxsi_20200915_1649.txt'
   # hale_fractions_file_12hr = '/Users/ainglis/physics/mission_simulator/out/new_run_n1000_hale_with_srs_corrections_12hr_noshift/mflare_observing_fractions_foxsi_20200915_1733.txt'
   # flare_index_fractions_file_12hr = '/Users/ainglis/physics/mission_simulator/out/new_run_n1000_flare_index_with_srs_corrections_12hr_noshift/mflare_observing_fractions_foxsi_20200915_1946.txt'
   # random_fractions_file_12hr = '/Users/ainglis/physics/mission_simulator/out/n1000_random_with_srs_corrections_12hr_noshift_random_fixed/mflare_observing_fractions_foxsi_20201008_1654.txt'

   # mcintosh_fractions_file_1hr = '/Users/ainglis/physics/mission_simulator/out/new_run_n1000_mcintosh_with_srs_corrections_1hr_noshift/mflare_observing_fractions_foxsi_20200915_1143.txt'
   # hale_fractions_file_1hr = '/Users/ainglis/physics/mission_simulator/out/new_run_n1000_hale_with_srs_corrections_1hr_noshift/mflare_observing_fractions_foxsi_20200915_1216.txt'
   # flare_index_fractions_file_1hr = '/Users/ainglis/physics/mission_simulator/out/new_run_n1000_flare_index_with_srs_corrections_1hr_noshift/mflare_observing_fractions_foxsi_20200915_1411.txt'
   # random_fractions_file_1hr = '/Users/ainglis/physics/mission_simulator/out/n1000_random_with_srs_corrections_1hr_noshift_random_fixed/mflare_observing_fractions_foxsi_20201009_1127.txt'

   # mcintosh_fractions_file_1hr_geo = '/Users/ainglis/physics/mission_simulator/out/new_run_n1000_mcintosh_with_srs_corrections_1hr_geo_noshift/mflare_observing_fractions_foxsi_20200914_1646.txt'
   # hale_fractions_file_1hr_geo = '/Users/ainglis/physics/mission_simulator/out/new_run_n1000_hale_with_srs_corrections_1hr_geo_noshift/mflare_observing_fractions_foxsi_20200914_1729.txt'
   # flare_index_fractions_file_1hr_geo = '/Users/ainglis/physics/mission_simulator/out/new_run_n1000_flare_index_with_srs_corrections_1hr_geo_noshift/mflare_observing_fractions_foxsi_20200914_1857.txt'
   # random_fractions_file_1hr_geo = '/Users/ainglis/physics/mission_simulator/out/n1000_random_with_srs_corrections_1hr_geo_noshift_random_fixed/mflare_observing_fractions_foxsi_20201009_1610.txt'


  #  mcintosh_mflare_24hr = np.loadtxt('/Users/ainglis/physics/mission_simulator/out/new_run_n1000_mcintosh_with_srs_corrections_24hr_noshift/mflare_observing_totals_foxsi_20200916_1139.txt')
  #  mcintosh_xflare_24hr = np.loadtxt('/Users/ainglis/physics/mission_simulator/out/new_run_n1000_mcintosh_with_srs_corrections_24hr_noshift/xflare_observing_totals_foxsi_20200916_1139.txt')
    mcintosh_mflare_24hr = np.loadtxt('/Users/ainglis/physics/mission_simulator/out/n1000_mcintosh_bloomfield_with_srs_corrections_24hr_noshift/mflare_observing_totals_foxsi_20201221_1300.txt')
    mcintosh_xflare_24hr = np.loadtxt('/Users/ainglis/physics/mission_simulator/out/n1000_mcintosh_bloomfield_with_srs_corrections_24hr_noshift/xflare_observing_totals_foxsi_20201221_1300.txt')
    hale_mflare_24hr = np.loadtxt('/Users/ainglis/physics/mission_simulator/out/new_run_n1000_hale_with_srs_corrections_24hr_noshift/mflare_observing_totals_foxsi_20200916_1204.txt')
    hale_xflare_24hr = np.loadtxt('/Users/ainglis/physics/mission_simulator/out/new_run_n1000_hale_with_srs_corrections_24hr_noshift/xflare_observing_totals_foxsi_20200916_1204.txt')
    flare_index_mflare_24hr = np.loadtxt('/Users/ainglis/physics/mission_simulator/out/new_run_n1000_flare_index_with_srs_corrections_24hr_noshift/mflare_observing_totals_foxsi_20200916_1346.txt')
    flare_index_xflare_24hr = np.loadtxt('/Users/ainglis/physics/mission_simulator/out/new_run_n1000_flare_index_with_srs_corrections_24hr_noshift/xflare_observing_totals_foxsi_20200916_1346.txt')
    random_mflare_24hr = np.loadtxt('/Users/ainglis/physics/mission_simulator/out/n1000_random_with_srs_corrections_24hr_noshift_random_fixed/mflare_observing_totals_foxsi_20201008_1949.txt')
    random_xflare_24hr = np.loadtxt('/Users/ainglis/physics/mission_simulator/out/n1000_random_with_srs_corrections_24hr_noshift_random_fixed/xflare_observing_totals_foxsi_20201008_1949.txt')
    
   # mcintosh_mflare_12hr = np.loadtxt('/Users/ainglis/physics/mission_simulator/out/new_run_n1000_mcintosh_with_srs_corrections_12hr_noshift/mflare_observing_totals_foxsi_20200915_1649.txt')
   # mcintosh_xflare_12hr = np.loadtxt('/Users/ainglis/physics/mission_simulator/out/new_run_n1000_mcintosh_with_srs_corrections_12hr_noshift/xflare_observing_totals_foxsi_20200915_1649.txt')
    mcintosh_mflare_12hr = np.loadtxt('/Users/ainglis/physics/mission_simulator/out/n1000_mcintosh_bloomfield_with_srs_corrections_12hr_noshift/mflare_observing_totals_foxsi_20201210_1624.txt')
    mcintosh_xflare_12hr = np.loadtxt('/Users/ainglis/physics/mission_simulator/out/n1000_mcintosh_bloomfield_with_srs_corrections_12hr_noshift/xflare_observing_totals_foxsi_20201210_1624.txt') 
    hale_mflare_12hr = np.loadtxt('/Users/ainglis/physics/mission_simulator/out/new_run_n1000_hale_with_srs_corrections_12hr_noshift/mflare_observing_totals_foxsi_20200915_1733.txt')
    hale_xflare_12hr = np.loadtxt('/Users/ainglis/physics/mission_simulator/out/new_run_n1000_hale_with_srs_corrections_12hr_noshift/xflare_observing_totals_foxsi_20200915_1733.txt')
    flare_index_mflare_12hr = np.loadtxt('/Users/ainglis/physics/mission_simulator/out/new_run_n1000_flare_index_with_srs_corrections_12hr_noshift/mflare_observing_totals_foxsi_20200915_1946.txt')
    flare_index_xflare_12hr = np.loadtxt('/Users/ainglis/physics/mission_simulator/out/new_run_n1000_flare_index_with_srs_corrections_12hr_noshift/xflare_observing_totals_foxsi_20200915_1946.txt')
    random_mflare_12hr = np.loadtxt('/Users/ainglis/physics/mission_simulator/out/n1000_random_with_srs_corrections_12hr_noshift_random_fixed/mflare_observing_totals_foxsi_20201008_1654.txt')
    random_xflare_12hr = np.loadtxt('/Users/ainglis/physics/mission_simulator/out/n1000_random_with_srs_corrections_12hr_noshift_random_fixed/xflare_observing_totals_foxsi_20201008_1654.txt')

   # mcintosh_mflare_1hr = np.loadtxt('/Users/ainglis/physics/mission_simulator/out/new_run_n1000_mcintosh_with_srs_corrections_1hr_noshift/mflare_observing_totals_foxsi_20200915_1143.txt')
   # mcintosh_xflare_1hr = np.loadtxt('/Users/ainglis/physics/mission_simulator/out/new_run_n1000_mcintosh_with_srs_corrections_1hr_noshift/xflare_observing_totals_foxsi_20200915_1143.txt')
    mcintosh_mflare_1hr = np.loadtxt('/Users/ainglis/physics/mission_simulator/out/n1000_mcintosh_bloomfield_with_srs_corrections_1hr_noshift/mflare_observing_totals_foxsi_20201221_1456.txt')
    mcintosh_xflare_1hr = np.loadtxt('/Users/ainglis/physics/mission_simulator/out/n1000_mcintosh_bloomfield_with_srs_corrections_1hr_noshift/xflare_observing_totals_foxsi_20201221_1456.txt')
    hale_mflare_1hr = np.loadtxt('/Users/ainglis/physics/mission_simulator/out/new_run_n1000_hale_with_srs_corrections_1hr_noshift/mflare_observing_totals_foxsi_20200915_1216.txt')
    hale_xflare_1hr = np.loadtxt('/Users/ainglis/physics/mission_simulator/out/new_run_n1000_hale_with_srs_corrections_1hr_noshift/xflare_observing_totals_foxsi_20200915_1216.txt')
    flare_index_mflare_1hr = np.loadtxt('/Users/ainglis/physics/mission_simulator/out/new_run_n1000_flare_index_with_srs_corrections_1hr_noshift/mflare_observing_totals_foxsi_20200915_1411.txt')
    flare_index_xflare_1hr = np.loadtxt('/Users/ainglis/physics/mission_simulator/out/new_run_n1000_flare_index_with_srs_corrections_1hr_noshift/xflare_observing_totals_foxsi_20200915_1411.txt')
    random_mflare_1hr = np.loadtxt('/Users/ainglis/physics/mission_simulator/out/n1000_random_with_srs_corrections_1hr_noshift_random_fixed/mflare_observing_totals_foxsi_20201009_1127.txt')
    random_xflare_1hr = np.loadtxt('/Users/ainglis/physics/mission_simulator/out/n1000_random_with_srs_corrections_1hr_noshift_random_fixed/xflare_observing_totals_foxsi_20201009_1127.txt')

  #  mcintosh_mflare_1hr_geo = np.loadtxt('/Users/ainglis/physics/mission_simulator/out/new_run_n1000_mcintosh_with_srs_corrections_1hr_geo_noshift/mflare_observing_totals_foxsi_20200914_1646.txt')
  #  mcintosh_xflare_1hr_geo = np.loadtxt('/Users/ainglis/physics/mission_simulator/out/new_run_n1000_mcintosh_with_srs_corrections_1hr_geo_noshift/xflare_observing_totals_foxsi_20200914_1646.txt')
    mcintosh_mflare_1hr_geo = np.loadtxt('/Users/ainglis/physics/mission_simulator/out/n1000_mcintosh_bloomfield_with_srs_corrections_1hr_geo_noshift/mflare_observing_totals_foxsi_20201221_1647.txt')
    mcintosh_xflare_1hr_geo = np.loadtxt('/Users/ainglis/physics/mission_simulator/out/n1000_mcintosh_bloomfield_with_srs_corrections_1hr_geo_noshift/xflare_observing_totals_foxsi_20201221_1647.txt')
    hale_mflare_1hr_geo = np.loadtxt('/Users/ainglis/physics/mission_simulator/out/new_run_n1000_hale_with_srs_corrections_1hr_geo_noshift/mflare_observing_totals_foxsi_20200914_1729.txt')
    hale_xflare_1hr_geo = np.loadtxt('/Users/ainglis/physics/mission_simulator/out/new_run_n1000_hale_with_srs_corrections_1hr_geo_noshift/xflare_observing_totals_foxsi_20200914_1729.txt')
    flare_index_mflare_1hr_geo = np.loadtxt('/Users/ainglis/physics/mission_simulator/out/new_run_n1000_flare_index_with_srs_corrections_1hr_geo_noshift/mflare_observing_totals_foxsi_20200914_1857.txt')
    flare_index_xflare_1hr_geo = np.loadtxt('/Users/ainglis/physics/mission_simulator/out/new_run_n1000_flare_index_with_srs_corrections_1hr_geo_noshift/xflare_observing_totals_foxsi_20200914_1857.txt')
    random_mflare_1hr_geo = np.loadtxt('/Users/ainglis/physics/mission_simulator/out/n1000_random_with_srs_corrections_1hr_geo_noshift_random_fixed/mflare_observing_totals_foxsi_20201009_1610.txt')
    random_xflare_1hr_geo = np.loadtxt('/Users/ainglis/physics/mission_simulator/out/n1000_random_with_srs_corrections_1hr_geo_noshift_random_fixed/xflare_observing_totals_foxsi_20201009_1610.txt')

    mcintosh_mflare_24hr_monwedfri = np.loadtxt('/Users/ainglis/physics/mission_simulator/out/n1000_mcintosh_bloomfield_with_srs_corrections_24hr_noshift_monwedfri/mflare_observing_totals_foxsi_20210128_1518.txt')
    mcintosh_xflare_24hr_monwedfri = np.loadtxt('/Users/ainglis/physics/mission_simulator/out/n1000_mcintosh_bloomfield_with_srs_corrections_24hr_noshift_monwedfri/xflare_observing_totals_foxsi_20210128_1518.txt')
    hale_mflare_24hr_monwedfri = np.loadtxt('/Users/ainglis/physics/mission_simulator/out/n1000_hale_with_srs_corrections_24hr_noshift_monwedfri/mflare_observing_totals_foxsi_20210129_1113.txt')
    hale_xflare_24hr_monwedfri = np.loadtxt('/Users/ainglis/physics/mission_simulator/out/n1000_hale_with_srs_corrections_24hr_noshift_monwedfri/xflare_observing_totals_foxsi_20210129_1113.txt')
    flare_index_mflare_24hr_monwedfri = np.loadtxt('/Users/ainglis/physics/mission_simulator/out/n1000_flare_index_with_srs_corrections_24hr_noshift_monwedfri/mflare_observing_totals_foxsi_20210128_1723.txt')
    flare_index_xflare_24hr_monwedfri = np.loadtxt('/Users/ainglis/physics/mission_simulator/out/n1000_flare_index_with_srs_corrections_24hr_noshift_monwedfri/xflare_observing_totals_foxsi_20210128_1723.txt')
    random_mflare_24hr_monwedfri = np.loadtxt('/Users/ainglis/physics/mission_simulator/out/n1000_random_with_srs_corrections_24hr_noshift_monwedfri/mflare_observing_totals_foxsi_20210204_1456.txt')
    random_xflare_24hr_monwedfri = np.loadtxt('/Users/ainglis/physics/mission_simulator/out/n1000_random_with_srs_corrections_24hr_noshift_monwedfri/xflare_observing_totals_foxsi_20210204_1456.txt')

    mcintosh_mflare_24hr_wed = np.loadtxt('/Users/ainglis/physics/mission_simulator/out/n1000_mcintosh_bloomfield_with_srs_corrections_24hr_noshift_wed/mflare_observing_totals_foxsi_20210201_1519.txt')
    mcintosh_xflare_24hr_wed = np.loadtxt('/Users/ainglis/physics/mission_simulator/out/n1000_mcintosh_bloomfield_with_srs_corrections_24hr_noshift_wed/xflare_observing_totals_foxsi_20210201_1519.txt')
    hale_mflare_24hr_wed = np.loadtxt('/Users/ainglis/physics/mission_simulator/out/n1000_hale_with_srs_corrections_24hr_noshift_wed/mflare_observing_totals_foxsi_20210201_1112.txt')
    hale_xflare_24hr_wed = np.loadtxt('/Users/ainglis/physics/mission_simulator/out/n1000_hale_with_srs_corrections_24hr_noshift_wed/xflare_observing_totals_foxsi_20210201_1112.txt')
    flare_index_mflare_24hr_wed = np.loadtxt('/Users/ainglis/physics/mission_simulator/out/n1000_flare_index_with_srs_corrections_24hr_noshift_wed/mflare_observing_totals_foxsi_20210201_1314.txt')
    flare_index_xflare_24hr_wed = np.loadtxt('/Users/ainglis/physics/mission_simulator/out/n1000_flare_index_with_srs_corrections_24hr_noshift_wed/xflare_observing_totals_foxsi_20210201_1314.txt')
    random_mflare_24hr_wed = np.loadtxt('/Users/ainglis/physics/mission_simulator/out/n1000_random_with_srs_corrections_24hr_noshift_wed/mflare_observing_totals_foxsi_20210204_1757.txt')
    random_xflare_24hr_wed = np.loadtxt('/Users/ainglis/physics/mission_simulator/out/n1000_random_with_srs_corrections_24hr_noshift_wed/xflare_observing_totals_foxsi_20210204_1757.txt')
    
    
   # mcintosh_mflare_24hr = np.loadtxt(mcintosh_fractions_file_24hr)
   # hale_mflare_24hr = np.loadtxt(hale_fractions_file_24hr)
   # flare_index_mflare_24hr = np.loadtxt(flare_index_fractions_file_24hr)
   # random_mflare_24hr = np.loadtxt(random_fractions_file_24hr)
    
   # mcintosh_mflare_12hr = np.loadtxt(mcintosh_fractions_file_12hr)
   # hale_mflare_12hr = np.loadtxt(hale_fractions_file_12hr)
   # flare_index_mflare_12hr = np.loadtxt(flare_index_fractions_file_12hr)
   # random_mflare_12hr = np.loadtxt(random_fractions_file_12hr)
    
   # mcintosh_mflare_1hr = np.loadtxt(mcintosh_fractions_file_1hr)
   # hale_mflare_1hr = np.loadtxt(hale_fractions_file_1hr)
   # flare_index_mflare_1hr = np.loadtxt(flare_index_fractions_file_1hr)
   # random_mflare_1hr = np.loadtxt(random_fractions_file_1hr)

   # mcintosh_mflare_1hr_geo = np.loadtxt(mcintosh_fractions_file_1hr_geo)
   # hale_mflare_1hr_geo = np.loadtxt(hale_fractions_file_1hr_geo)
   # flare_index_mflare_1hr_geo = np.loadtxt(flare_index_fractions_file_1hr_geo)
   # random_mflare_1hr_geo = np.loadtxt(random_fractions_file_1hr_geo)
    

  
    fontsize=8
    labelsize=8
    mean_line = {"color": 'black', "linestyle": 'dashed'}
    percentile_line  = {"color": 'red', "linestyle": 'dotted'}
    grid_style = {"linestyle": "dotted"}
    x0 = 0.1
    x1 = 0.8
    
    
    plt.figure(1,figsize=(12,12))
    plt.subplots_adjust(bottom=0.1,top=0.9,left=0.05,right=0.95)

    # first row

    plt.subplot(4,4,1)
  #  flare_fractions = mcintosh_mflare_24hr
    flare_fractions = (mcintosh_mflare_24hr + mcintosh_xflare_24hr) / 589.
    bins = regular_hist_bins(flare_fractions, 0.01)
    weights = np.ones_like(flare_fractions)/len(flare_fractions)
    plt.hist(flare_fractions, bins=bins, weights=weights,color='blue')
    plt.ylabel('probability', fontsize=fontsize)
   # plt.xlabel('fraction of all M-class flares observed', fontsize=fontsize)
    plt.tick_params(labelsize=labelsize)
    mn = np.mean(flare_fractions)
    plt.axvline(mn, label=f'mean={mn:.3f}', **mean_line)
    percentiles = np.percentile(flare_fractions, [2.5, 97.5])
    plt.axvline(percentiles[0], **percentile_line)
    plt.axvline(percentiles[1], label=f'95% interval', **percentile_line) # is {percentiles[0]:.3f} to {percentiles[1]:.3f}', **percentile_line)
    plt.legend(fontsize=8)
    plt.grid(**grid_style)
    plt.xlim([x0, x1])
    plt.ylim([0,0.3])
    plt.title('24 hr delay',fontsize=8)
    plt.text(0.6,0.025,'Mcintosh')
    
    
    plt.subplot(4,4,2)
   # flare_fractions = mcintosh_mflare_12hr
    flare_fractions = (mcintosh_mflare_12hr + mcintosh_xflare_12hr) / 589.
    bins = regular_hist_bins(flare_fractions, 0.01)
    weights = np.ones_like(flare_fractions)/len(flare_fractions)
    plt.hist(flare_fractions, bins=bins, weights=weights,color='blue')
    plt.ylabel('probability', fontsize=fontsize)
   # plt.xlabel('fraction of all M-class flares observed', fontsize=fontsize)
    plt.tick_params(labelsize=labelsize)
    mn = np.mean(flare_fractions)
    plt.axvline(mn, label=f'mean={mn:.3f}', **mean_line)
    percentiles = np.percentile(flare_fractions, [2.5, 97.5])
    plt.axvline(percentiles[0], **percentile_line)
    plt.axvline(percentiles[1], label=f'95% interval', **percentile_line) # is {percentiles[0]:.3f} to {percentiles[1]:.3f}', **percentile_line)
    plt.legend(fontsize=8)
    plt.grid(**grid_style)
    plt.xlim([x0, x1])
    plt.ylim([0,0.3])
    plt.title('12 hr delay',fontsize=8)
    plt.text(0.6,0.025,'Mcintosh')

    plt.subplot(4,4,3)
    #flare_fractions = mcintosh_mflare_1hr
    flare_fractions = (mcintosh_mflare_1hr + mcintosh_xflare_1hr) / 589.
    bins = regular_hist_bins(flare_fractions, 0.01)
    weights = np.ones_like(flare_fractions)/len(flare_fractions)
    plt.hist(flare_fractions, bins=bins, weights=weights,color='blue')
    plt.ylabel('probability', fontsize=fontsize)
   # plt.xlabel('fraction of all M-class flares observed', fontsize=fontsize)
    plt.tick_params(labelsize=labelsize)
    mn = np.mean(flare_fractions)
    plt.axvline(mn, label=f'mean={mn:.3f}', **mean_line)
    percentiles = np.percentile(flare_fractions, [2.5, 97.5])
    plt.axvline(percentiles[0], **percentile_line)
    plt.axvline(percentiles[1], label=f'95% interval', **percentile_line) # is {percentiles[0]:.3f} to {percentiles[1]:.3f}', **percentile_line)
    plt.legend(fontsize=8)
    plt.grid(**grid_style)
    plt.xlim([x0, x1])
    plt.ylim([0,0.3])
    plt.title('1 hr delay',fontsize=8)
    plt.text(0.6,0.025,'Mcintosh')

    plt.subplot(4,4,4)
   # flare_fractions = mcintosh_mflare_1hr_geo
    flare_fractions = (mcintosh_mflare_1hr_geo + mcintosh_xflare_1hr_geo) / 589.
    bins = regular_hist_bins(flare_fractions, 0.01)
    weights = np.ones_like(flare_fractions)/len(flare_fractions)
    plt.hist(flare_fractions, bins=bins, weights=weights,color='blue')
    plt.ylabel('probability', fontsize=fontsize)
   # plt.xlabel('fraction of all M-class flares observed', fontsize=fontsize)
    plt.tick_params(labelsize=labelsize)
    mn = np.mean(flare_fractions)
    plt.axvline(mn, label=f'mean={mn:.3f}', **mean_line)
    percentiles = np.percentile(flare_fractions, [2.5, 97.5])
    plt.axvline(percentiles[0], **percentile_line)
    plt.axvline(percentiles[1], label=f'95% interval', **percentile_line) # is {percentiles[0]:.3f} to {percentiles[1]:.3f}', **percentile_line)
    plt.legend(fontsize=8)
    plt.grid(**grid_style)
    plt.xlim([x0, x1])
    plt.ylim([0,1.0])
    plt.title('1 hr delay, no eclipse, SAA',fontsize=8)
    plt.text(0.6,0.083,'Mcintosh')


    
    # second row
    plt.subplot(4,4,5)
   # flare_fractions = hale_mflare_24hr
    flare_fractions = (hale_mflare_24hr + hale_xflare_24hr) / 589.0
    bins = regular_hist_bins(flare_fractions, 0.01)
    weights = np.ones_like(flare_fractions)/len(flare_fractions)
    plt.hist(flare_fractions, bins=bins, weights=weights,color='green')
    plt.ylabel('probability', fontsize=fontsize)
  #  plt.xlabel('fraction of all M-class flares observed', fontsize=fontsize)
    plt.tick_params(labelsize=labelsize)
    mn = np.mean(flare_fractions)
    plt.axvline(mn, label=f'mean={mn:.3f}', **mean_line)
    percentiles = np.percentile(flare_fractions, [2.5, 97.5])
    plt.axvline(percentiles[0], **percentile_line)
    plt.axvline(percentiles[1], label=f'95% interval', **percentile_line) # is {percentiles[0]:.3f} to {percentiles[1]:.3f}', **percentile_line)
    plt.legend(fontsize=8)
    plt.grid(**grid_style)
    plt.xlim([x0, x1])
    plt.ylim([0,0.3])
    plt.title('24 hr delay',fontsize=8)
    plt.text(0.55,0.025,'Hale + spots')
    
    plt.subplot(4,4,6)
    #flare_fractions = hale_mflare_12hr
    flare_fractions = (hale_mflare_12hr + hale_xflare_12hr) / 589.0
    bins = regular_hist_bins(flare_fractions, 0.01)
    weights = np.ones_like(flare_fractions)/len(flare_fractions)
    plt.hist(flare_fractions, bins=bins, weights=weights,color='green')
    plt.ylabel('probability', fontsize=fontsize)
  #  plt.xlabel('fraction of all M-class flares observed', fontsize=fontsize)
    plt.tick_params(labelsize=labelsize)
    mn = np.mean(flare_fractions)
    plt.axvline(mn, label=f'mean={mn:.3f}', **mean_line)
    percentiles = np.percentile(flare_fractions, [2.5, 97.5])
    plt.axvline(percentiles[0], **percentile_line)
    plt.axvline(percentiles[1], label=f'95% interval', **percentile_line) # is {percentiles[0]:.3f} to {percentiles[1]:.3f}', **percentile_line)
    plt.legend(fontsize=8)
    plt.grid(**grid_style)
    plt.xlim([x0, x1])
    plt.ylim([0,0.3])
    plt.title('12 hr delay',fontsize=8)
    plt.text(0.55,0.025,'Hale + spots')

    plt.subplot(4,4,7)
    #flare_fractions = hale_mflare_1hr
    flare_fractions = (hale_mflare_1hr + hale_xflare_1hr) / 589.0
    bins = regular_hist_bins(flare_fractions, 0.01)
    weights = np.ones_like(flare_fractions)/len(flare_fractions)
    plt.hist(flare_fractions, bins=bins, weights=weights,color='green')
    plt.ylabel('probability', fontsize=fontsize)
  #  plt.xlabel('fraction of all M-class flares observed', fontsize=fontsize)
    plt.tick_params(labelsize=labelsize)
    mn = np.mean(flare_fractions)
    plt.axvline(mn, label=f'mean={mn:.3f}', **mean_line)
    percentiles = np.percentile(flare_fractions, [2.5, 97.5])
    plt.axvline(percentiles[0], **percentile_line)
    plt.axvline(percentiles[1], label=f'95% interval', **percentile_line) # is {percentiles[0]:.3f} to {percentiles[1]:.3f}', **percentile_line)
    plt.legend(fontsize=8)
    plt.grid(**grid_style)
    plt.xlim([x0, x1])
    plt.ylim([0,0.3])
    plt.title('1 hr delay',fontsize=8)
    plt.text(0.55,0.025,'Hale + spots')

    plt.subplot(4,4,8)
   # flare_fractions = hale_mflare_1hr_geo
    flare_fractions = (hale_mflare_1hr_geo + hale_xflare_1hr_geo) / 589.0
    bins = regular_hist_bins(flare_fractions, 0.01)
    weights = np.ones_like(flare_fractions)/len(flare_fractions)
    plt.hist(flare_fractions, bins=bins, weights=weights,color='green')
    plt.ylabel('probability', fontsize=fontsize)
  #  plt.xlabel('fraction of all M-class flares observed', fontsize=fontsize)
    plt.tick_params(labelsize=labelsize)
    mn = np.mean(flare_fractions)
    plt.axvline(mn, label=f'mean={mn:.3f}', **mean_line)
    percentiles = np.percentile(flare_fractions, [2.5, 97.5])
    plt.axvline(percentiles[0], **percentile_line)
    plt.axvline(percentiles[1], label=f'95% interval', **percentile_line) # is {percentiles[0]:.3f} to {percentiles[1]:.3f}', **percentile_line)
    plt.legend(fontsize=8)
    plt.grid(**grid_style)
    plt.xlim([x0, x1])
    plt.ylim([0,1.0])
    plt.title('1 hr delay, no eclipse, SAA',fontsize=8)
    plt.text(0.55,0.083,'Hale + spots')
    

    #third row
    plt.subplot(4,4,9)
   # flare_fractions = flare_index_mflare_24hr
    flare_fractions = (flare_index_mflare_24hr + flare_index_xflare_24hr) / 589.0
    bins = regular_hist_bins(flare_fractions, 0.01)
    weights = np.ones_like(flare_fractions)/len(flare_fractions)
    plt.hist(flare_fractions, bins=bins, weights=weights,color='purple')
    plt.ylabel('probability', fontsize=fontsize)
  #  plt.xlabel('fraction of all M-class flares observed', fontsize=fontsize)
    plt.tick_params(labelsize=labelsize)
    mn = np.mean(flare_fractions)
    plt.axvline(mn, label=f'mean={mn:.3f}', **mean_line)
    percentiles = np.percentile(flare_fractions, [2.5, 97.5])
    plt.axvline(percentiles[0], **percentile_line)
    plt.axvline(percentiles[1], label=f'95% interval', **percentile_line)# is {percentiles[0]:.3f} to {percentiles[1]:.3f}', **percentile_line)
    plt.legend(fontsize=8)
    plt.grid(**grid_style)
    plt.xlim([x0, x1])
    plt.ylim([0,0.3])
    plt.title('24 hr delay',fontsize=8)
    plt.text(0.57,0.025,'Flare index')

    plt.subplot(4,4,10)
   # flare_fractions = flare_index_mflare_12hr
    flare_fractions = (flare_index_mflare_12hr + flare_index_xflare_12hr) / 589.0
    bins = regular_hist_bins(flare_fractions, 0.01)
    weights = np.ones_like(flare_fractions)/len(flare_fractions)
    plt.hist(flare_fractions, bins=bins, weights=weights,color='purple')
    plt.ylabel('probability', fontsize=fontsize)
  #  plt.xlabel('fraction of all M-class flares observed', fontsize=fontsize)
    plt.tick_params(labelsize=labelsize)
    mn = np.mean(flare_fractions)
    plt.axvline(mn, label=f'mean={mn:.3f}', **mean_line)
    percentiles = np.percentile(flare_fractions, [2.5, 97.5])
    plt.axvline(percentiles[0], **percentile_line)
    plt.axvline(percentiles[1], label=f'95% interval', **percentile_line)# is {percentiles[0]:.3f} to {percentiles[1]:.3f}', **percentile_line)
    plt.legend(fontsize=8)
    plt.grid(**grid_style)
    plt.xlim([x0, x1])
    plt.ylim([0,0.3])
    plt.title('12 hr delay',fontsize=8)
    plt.text(0.57,0.025,'Flare index')

    plt.subplot(4,4,11)
   # flare_fractions = flare_index_mflare_1hr
    flare_fractions = (flare_index_mflare_1hr + flare_index_xflare_1hr) / 589.0
    bins = regular_hist_bins(flare_fractions, 0.01)
    weights = np.ones_like(flare_fractions)/len(flare_fractions)
    plt.hist(flare_fractions, bins=bins, weights=weights,color='purple')
    plt.ylabel('probability', fontsize=fontsize)
  #  plt.xlabel('fraction of all M-class flares observed', fontsize=fontsize)
    plt.tick_params(labelsize=labelsize)
    mn = np.mean(flare_fractions)
    plt.axvline(mn, label=f'mean={mn:.3f}', **mean_line)
    percentiles = np.percentile(flare_fractions, [2.5, 97.5])
    plt.axvline(percentiles[0], **percentile_line)
    plt.axvline(percentiles[1], label=f'95% interval', **percentile_line)# is {percentiles[0]:.3f} to {percentiles[1]:.3f}', **percentile_line)
    plt.legend(fontsize=8)
    plt.grid(**grid_style)
    plt.xlim([x0, x1])
    plt.ylim([0,0.3])
    plt.title('1 hr delay',fontsize=8)
    plt.text(0.57,0.025,'Flare index')

    plt.subplot(4,4,12)
    #flare_fractions = flare_index_mflare_1hr_geo
    flare_fractions = (flare_index_mflare_1hr_geo + flare_index_xflare_1hr_geo) / 589.0
    bins = regular_hist_bins(flare_fractions, 0.01)
    weights = np.ones_like(flare_fractions)/len(flare_fractions)
    plt.hist(flare_fractions, bins=bins, weights=weights,color='purple')
    plt.ylabel('probability', fontsize=fontsize)
  #  plt.xlabel('fraction of all M-class flares observed', fontsize=fontsize)
    plt.tick_params(labelsize=labelsize)
    mn = np.mean(flare_fractions)
    plt.axvline(mn, label=f'mean={mn:.3f}', **mean_line)
    percentiles = np.percentile(flare_fractions, [2.5, 97.5])
    plt.axvline(percentiles[0], **percentile_line)
    plt.axvline(percentiles[1], label=f'95% interval', **percentile_line)# is {percentiles[0]:.3f} to {percentiles[1]:.3f}', **percentile_line)
    plt.legend(fontsize=8)
    plt.grid(**grid_style)
    plt.xlim([x0, x1])
    plt.ylim([0,1.0])
    plt.title('1 hr delay, no eclipse, SAA',fontsize=8)
    plt.text(0.15,0.083,'Flare index')

    
    # fourth row
    plt.subplot(4,4,13)
    #flare_fractions = random_mflare_24hr
    flare_fractions = (random_mflare_24hr + random_xflare_24hr) / 589.0
    bins = regular_hist_bins(flare_fractions, 0.01)
    weights = np.ones_like(flare_fractions)/len(flare_fractions)
    plt.hist(flare_fractions, bins=bins, weights=weights,color='red')
    plt.ylabel('probability', fontsize=fontsize)
    plt.xlabel('fraction of all M-class flares observed', fontsize=fontsize)
    plt.tick_params(labelsize=labelsize)
    mn = np.mean(flare_fractions)
    plt.axvline(mn, label=f'mean={mn:.3f}', **mean_line)
    percentiles = np.percentile(flare_fractions, [2.5, 97.5])
    plt.axvline(percentiles[0], **percentile_line)
    plt.axvline(percentiles[1], label=f'95% interval', **percentile_line)# is {percentiles[0]:.3f} to {percentiles[1]:.3f}', **percentile_line)
    plt.legend(fontsize=8)
    plt.grid(**grid_style)
    plt.xlim([x0, x1])
    plt.ylim([0,0.3])
    plt.title('24 hr delay',fontsize=8)
    plt.text(0.6,0.025,'Random')
    
    plt.subplot(4,4,14)
   # flare_fractions = random_mflare_12hr
    flare_fractions = (random_mflare_12hr + random_xflare_12hr) / 589.0
    bins = regular_hist_bins(flare_fractions, 0.01)
    weights = np.ones_like(flare_fractions)/len(flare_fractions)
    plt.hist(flare_fractions, bins=bins, weights=weights,color='red')
    plt.ylabel('probability', fontsize=fontsize)
    plt.xlabel('fraction of all M-class flares observed', fontsize=fontsize)
    plt.tick_params(labelsize=labelsize)
    mn = np.mean(flare_fractions)
    plt.axvline(mn, label=f'mean={mn:.3f}', **mean_line)
    percentiles = np.percentile(flare_fractions, [2.5, 97.5])
    plt.axvline(percentiles[0], **percentile_line)
    plt.axvline(percentiles[1], label=f'95% interval', **percentile_line)# is {percentiles[0]:.3f} to {percentiles[1]:.3f}', **percentile_line)
    plt.legend(fontsize=8)
    plt.grid(**grid_style)
    plt.xlim([x0, x1])
    plt.ylim([0,0.3])
    plt.title('12 hr delay',fontsize=8)
    plt.text(0.6,0.025,'Random')
    
    plt.subplot(4,4,15)
 #   flare_fractions = random_mflare_1hr
    flare_fractions = (random_mflare_1hr + random_xflare_1hr) / 589.0
    bins = regular_hist_bins(flare_fractions, 0.01)
    weights = np.ones_like(flare_fractions)/len(flare_fractions)
    plt.hist(flare_fractions, bins=bins, weights=weights,color='red')
    plt.ylabel('probability', fontsize=fontsize)
    plt.xlabel('fraction of all M-class flares observed', fontsize=fontsize)
    plt.tick_params(labelsize=labelsize)
    mn = np.mean(flare_fractions)
    plt.axvline(mn, label=f'mean={mn:.3f}', **mean_line)
    percentiles = np.percentile(flare_fractions, [2.5, 97.5])
    plt.axvline(percentiles[0], **percentile_line)
    plt.axvline(percentiles[1], label=f'95% interval', **percentile_line)# is {percentiles[0]:.3f} to {percentiles[1]:.3f}', **percentile_line)
    plt.legend(fontsize=8)
    plt.grid(**grid_style)
    plt.xlim([x0, x1])
    plt.ylim([0,0.3])
    plt.title('1 hr delay',fontsize=8)
    plt.text(0.6,0.025,'Random')

    plt.subplot(4,4,16)
   # flare_fractions = random_mflare_1hr_geo
    flare_fractions = (random_mflare_1hr_geo + random_xflare_1hr_geo) / 589.0
    bins = regular_hist_bins(flare_fractions, 0.01)
    weights = np.ones_like(flare_fractions)/len(flare_fractions)
    plt.hist(flare_fractions, bins=bins, weights=weights,color='red')
    plt.ylabel('probability', fontsize=fontsize)
    plt.xlabel('fraction of all M-class flares observed', fontsize=fontsize)
    plt.tick_params(labelsize=labelsize)
    mn = np.mean(flare_fractions)
    plt.axvline(mn, label=f'mean={mn:.3f}', **mean_line)
    percentiles = np.percentile(flare_fractions, [2.5, 97.5])
    plt.axvline(percentiles[0], **percentile_line)
    plt.axvline(percentiles[1], label=f'95% interval', **percentile_line)# is {percentiles[0]:.3f} to {percentiles[1]:.3f}', **percentile_line)
    plt.legend(fontsize=8)
    plt.grid(**grid_style)
    plt.xlim([x0, x1])
    plt.ylim([0,1.0])
    plt.title('1 hr delay, no eclipse, SAA',fontsize=8)
    plt.text(0.6,0.083,'Random')
    
    plt.savefig('performance_figure.pdf')
    


    plt.figure(2,figsize=(8,5))
    plt.subplots_adjust(bottom=0.1,top=0.9,left=0.1,right=0.90)

    plt.subplot(1,1,1)

    flare_fractions = (mcintosh_mflare_12hr + mcintosh_xflare_12hr) / 589.
    bins = regular_hist_bins(flare_fractions, 0.01)
    weights = np.ones_like(flare_fractions)/len(flare_fractions)
    plt.hist(flare_fractions, bins=bins, weights=weights)
    plt.ylabel('probability', fontsize=14)
    plt.xlabel('Fraction of >M1 flares observed', fontsize=14)
   # plt.xlabel('fraction of all M-class flares observed', fontsize=fontsize)
    plt.tick_params(labelsize=labelsize)
    mn = np.mean(flare_fractions)
    std = np.std(flare_fractions)
    plt.axvline(mn, label=f'mean={mn:.3f}', **mean_line)
   # percentiles = np.percentile(flare_fractions, [16.0, 84.0])
    plt.axvline(mn - std, **percentile_line, linewidth=2)
    plt.axvline(mn + std, label=f'1$\sigma$ interval', **percentile_line, linewidth=2) # is {percentiles[0]:.3f} to {percentiles[1]:.3f}', **percentile_line)
    plt.legend(fontsize=12)
    plt.grid(**grid_style)
    plt.xlim([0.3, 0.50])
    plt.ylim([0,0.3])
    plt.title('Mcintosh daily repoint with 12 hr delay (scenario M4)',fontsize=12)
    #plt.text(0.5,0.025,'Mcintosh')
    plt.tick_params(labelsize=12)

    plt.savefig('performance_example.pdf')

    plt.figure(3,figsize=(10,5))
    plt.subplots_adjust(bottom=0.1,top=0.9,left=0.1,right=0.90)
    
    plt.subplot(1,1,1)

    flare_fractions = (random_mflare_24hr_wed + random_xflare_24hr_wed) / 589.0
    mn = np.mean(flare_fractions)
    std = np.std(flare_fractions)
    plt.errorbar(1.0, mn, xerr = None, yerr = std, color='red', marker = 'o')
    
    flare_fractions = (random_mflare_24hr_monwedfri + random_xflare_24hr_monwedfri) / 589.0
    mn = np.mean(flare_fractions)
    std = np.std(flare_fractions)
    plt.errorbar(3.0, mn, xerr = None, yerr = std, color='red', marker = 'o')
    
    flare_fractions = (random_mflare_24hr + random_xflare_24hr) / 589.0
    mn = np.mean(flare_fractions)
    std = np.std(flare_fractions)
    plt.errorbar(5.0, mn, xerr = None, yerr = std, color='red', marker = 'o')
    
    flare_fractions = (random_mflare_12hr + random_xflare_12hr) / 589.0
    mn = np.mean(flare_fractions)
    std = np.std(flare_fractions)
    plt.errorbar(7.0, mn, xerr = None, yerr = std, color='red', marker = 'o')

    flare_fractions = (random_mflare_1hr + random_xflare_1hr) / 589.0
    mn = np.mean(flare_fractions)
    std = np.std(flare_fractions)
    plt.errorbar(9.0, mn, xerr = None, yerr = std, color = 'red', marker = 'o')

    flare_fractions = (random_mflare_1hr_geo + random_xflare_1hr_geo) / 589.0
    mn = np.mean(flare_fractions)
    std = np.std(flare_fractions)
    plt.errorbar(11.0, mn, xerr = None, yerr = std, color = 'red', marker = 'o', label = 'Random')


    flare_fractions = (mcintosh_mflare_24hr_wed + mcintosh_xflare_24hr_wed) / 589.0
    mn = np.mean(flare_fractions)
    std = np.std(flare_fractions)
    plt.errorbar(1.0, mn, xerr = None, yerr = std, color='blue', marker = 'o')
    
    flare_fractions = (mcintosh_mflare_24hr_monwedfri + mcintosh_xflare_24hr_monwedfri) / 589.0
    mn = np.mean(flare_fractions)
    std = np.std(flare_fractions)
    plt.errorbar(3.0, mn, xerr = None, yerr = std, color='blue', marker = 'o')

    flare_fractions = (mcintosh_mflare_24hr + mcintosh_xflare_24hr) / 589.0
    mn = np.mean(flare_fractions)
    std = np.std(flare_fractions)
    plt.errorbar(5.0, mn, xerr = None, yerr = std, color = 'blue', marker = 'o', label = 'Mcintosh')
    
    flare_fractions = (mcintosh_mflare_12hr + mcintosh_xflare_12hr) / 589.0
    mn = np.mean(flare_fractions)
    std = np.std(flare_fractions)
    plt.errorbar(7.0, mn, xerr = None, yerr = std, color = 'blue', marker = 'o')

    flare_fractions = (mcintosh_mflare_1hr + mcintosh_xflare_1hr) / 589.0
    mn = np.mean(flare_fractions)
    std = np.std(flare_fractions)
    plt.errorbar(9.0, mn, xerr = None, yerr = std, color = 'blue', marker = 'o')

    flare_fractions = (mcintosh_mflare_1hr_geo + mcintosh_xflare_1hr_geo) / 589.0
    mn = np.mean(flare_fractions)
    std = np.std(flare_fractions)
    plt.errorbar(11.0, mn, xerr = None, yerr = std, color = 'blue', marker = 'o')


    
    flare_fractions = (hale_mflare_24hr_wed + hale_xflare_24hr_wed) / 589.0
    mn = np.mean(flare_fractions)
    std = np.std(flare_fractions)
    plt.errorbar(1.0, mn, xerr = None, yerr = std, color='green', marker = 'o')
    
    flare_fractions = (hale_mflare_24hr_monwedfri + hale_xflare_24hr_monwedfri) / 589.0
    mn = np.mean(flare_fractions)
    std = np.std(flare_fractions)
    plt.errorbar(3.0, mn, xerr = None, yerr = std, color='green', marker = 'o')

    flare_fractions = (hale_mflare_24hr + hale_xflare_24hr) / 589.0
    mn = np.mean(flare_fractions)
    std = np.std(flare_fractions)
    plt.errorbar(5.0, mn, xerr = None, yerr = std, color = 'green', marker = 'o')
    
    flare_fractions = (hale_mflare_12hr + hale_xflare_12hr) / 589.0
    mn = np.mean(flare_fractions)
    std = np.std(flare_fractions)
    plt.errorbar(7.0, mn, xerr = None, yerr = std, color = 'green', marker = 'o')

    flare_fractions = (hale_mflare_1hr + hale_xflare_1hr) / 589.0
    mn = np.mean(flare_fractions)
    std = np.std(flare_fractions)
    plt.errorbar(9.0, mn, xerr = None, yerr = std, color = 'green', marker = 'o')

    flare_fractions = (hale_mflare_1hr_geo + hale_xflare_1hr_geo) / 589.0
    mn = np.mean(flare_fractions)
    std = np.std(flare_fractions)
    plt.errorbar(11.0, mn, xerr = None, yerr = std, color = 'green', marker = 'o', label = 'Hale + spots')



    flare_fractions = (flare_index_mflare_24hr_wed + flare_index_xflare_24hr_wed) / 589.0
    mn = np.mean(flare_fractions)
    std = np.std(flare_fractions)
    plt.errorbar(1.0, mn, xerr = None, yerr = std, color='purple', marker = 'o')
    
    flare_fractions = (flare_index_mflare_24hr_monwedfri + flare_index_xflare_24hr_monwedfri) / 589.0
    mn = np.mean(flare_fractions)
    std = np.std(flare_fractions)
    plt.errorbar(3.0, mn, xerr = None, yerr = std, color='purple', marker = 'o')
    
    flare_fractions = (flare_index_mflare_24hr + flare_index_xflare_24hr) / 589.0
    mn = np.mean(flare_fractions)
    std = np.std(flare_fractions)
    plt.errorbar(5.0, mn, xerr = None, yerr = std, color = 'purple', marker = 'o', label = 'Flare index')
    
    flare_fractions = (flare_index_mflare_12hr + flare_index_xflare_12hr) / 589.0
    mn = np.mean(flare_fractions)
    std = np.std(flare_fractions)
    plt.errorbar(7.0, mn, xerr = None, yerr = std, color = 'purple', marker = 'o')

    flare_fractions = (flare_index_mflare_1hr + flare_index_xflare_1hr) / 589.0
    mn = np.mean(flare_fractions)
    std = np.std(flare_fractions)
    plt.errorbar(9.0, mn, xerr = None, yerr = std, color = 'purple', marker = 'o')

    flare_fractions = (flare_index_mflare_1hr_geo + flare_index_xflare_1hr_geo) / 589.0
    mn = np.mean(flare_fractions)
    std = np.std(flare_fractions)
    plt.errorbar(11.0, mn, xerr = None, yerr = std, color = 'purple', marker = 'o')

    plt.ylabel('Fraction of >M1 flares observed', fontsize=14)
    plt.legend(fontsize=12)
    plt.xlim([0,12])
    plt.ylim([0,0.7])
    plt.grid(axis = 'y')

    xlabels = ['','24hr \n weekly', '', '24hr \n 3/week', '', '24 hr \n daily', '', '12 hr \n daily','', ' 1 hr \n daily ', '', '1hr geo \n daily', '']
    plt.xticks(np.linspace(0,12,13), xlabels)
    plt.tick_params(labelsize=12)

    plt.axvline(2.0, linestyle='dashed', color='black')
    plt.axvline(4.0, linestyle='dashed', color='black')
    plt.axvline(6.0, linestyle='dashed', color = 'black')
    plt.axvline(8.0, linestyle='dashed', color = 'black')
    plt.axvline(10.0, linestyle='dashed', color = 'black')

    plt.axvspan(10.0, 12.0, color='darkorange', alpha=0.5, lw=0)

    plt.savefig('performance_figure_new.pdf')
    plt.show()


def compare_mcintosh_strategies():

    mcintosh_mflare_12hr = np.loadtxt('/Users/ainglis/physics/mission_simulator/out/new_run_n1000_mcintosh_with_srs_corrections_12hr_noshift/mflare_observing_totals_foxsi_20200915_1649.txt')
    mcintosh_xflare_12hr = np.loadtxt('/Users/ainglis/physics/mission_simulator/out/new_run_n1000_mcintosh_with_srs_corrections_12hr_noshift/xflare_observing_totals_foxsi_20200915_1649.txt')

    bloomfield_mcintosh_mflare_12hr = np.loadtxt('/Users/ainglis/physics/mission_simulator/out/n1000_mcintosh_bloomfield_with_srs_corrections_12hr_noshift/mflare_observing_totals_foxsi_20201210_1624.txt')
    bloomfield_mcintosh_xflare_12hr = np.loadtxt('/Users/ainglis/physics/mission_simulator/out/n1000_mcintosh_bloomfield_with_srs_corrections_12hr_noshift/xflare_observing_totals_foxsi_20201210_1624.txt')

    bornmann_shaw_mcintosh_mflare_12hr = np.loadtxt('/Users/ainglis/physics/mission_simulator/out/n1000_mcintosh_bornmann_shaw_with_srs_corrections_12hr_noshift/mflare_observing_totals_foxsi_20201211_1315.txt')
    bornmann_shaw_mcintosh_xflare_12hr = np.loadtxt('/Users/ainglis/physics/mission_simulator/out/n1000_mcintosh_bornmann_shaw_with_srs_corrections_12hr_noshift/xflare_observing_totals_foxsi_20201211_1315.txt')
    

    fontsize=12
    labelsize=10
    mean_line = {"color": 'black', "linestyle": 'dashed'}
    percentile_line  = {"color": 'red', "linestyle": 'dotted'}
    grid_style = {"linestyle": "dotted"}
    x0 = 0.1
    x1 = 0.8
    
    
    plt.figure(1,figsize=(12,4))
    plt.subplots_adjust(bottom=0.12,top=0.9,left=0.05,right=0.95)
    
    plt.subplot(1,3,3)
   # flare_fractions = mcintosh_mflare_12hr
    flare_fractions = (mcintosh_mflare_12hr + mcintosh_xflare_12hr) / 589.
    bins = regular_hist_bins(flare_fractions, 0.01)
    weights = np.ones_like(flare_fractions)/len(flare_fractions)
    plt.hist(flare_fractions, bins=bins, weights=weights,color='blue')
    plt.ylabel('probability', fontsize=fontsize)
    plt.xlabel('fraction of flares >M1 observed', fontsize=fontsize)
    plt.tick_params(labelsize=labelsize)
    mn = np.mean(flare_fractions)
    plt.axvline(mn, label=f'mean={mn:.3f}', **mean_line)
    percentiles = np.percentile(flare_fractions, [2.5, 97.5])
    plt.axvline(percentiles[0], **percentile_line)
    plt.axvline(percentiles[1], label=f'95% interval', **percentile_line) # is {percentiles[0]:.3f} to {percentiles[1]:.3f}', **percentile_line)
    plt.legend(fontsize=8)
    plt.grid(**grid_style)
    plt.xlim([0.3, 0.6])
    plt.ylim([0,0.3])
    plt.title('\' Averaging\' Method, 12 hr delay',fontsize=8)
   # plt.text(0.6,0.025,'Mcintosh')

    
    plt.subplot(1,3,1)
   # flare_fractions = mcintosh_mflare_12hr
    flare_fractions = (bloomfield_mcintosh_mflare_12hr + bloomfield_mcintosh_xflare_12hr) / 589.
    bins = regular_hist_bins(flare_fractions, 0.01)
    weights = np.ones_like(flare_fractions)/len(flare_fractions)
    plt.hist(flare_fractions, bins=bins, weights=weights,color='blue')
    plt.ylabel('probability', fontsize=fontsize)
    plt.xlabel('fraction of flares >M1 observed', fontsize=fontsize)
    plt.tick_params(labelsize=labelsize)
    mn = np.mean(flare_fractions)
    plt.axvline(mn, label=f'mean={mn:.3f}', **mean_line)
    percentiles = np.percentile(flare_fractions, [2.5, 97.5])
    plt.axvline(percentiles[0], **percentile_line)
    plt.axvline(percentiles[1], label=f'95% interval', **percentile_line) # is {percentiles[0]:.3f} to {percentiles[1]:.3f}', **percentile_line)
    plt.legend(fontsize=8)
    plt.grid(**grid_style)
    plt.xlim([0.3, 0.6])
    plt.ylim([0,0.3])
    plt.title('Bloomfield et al. (2012), 12 hr delay',fontsize=8)
    #plt.text(0.6,0.025,'h')

    plt.subplot(1,3,2)
   # flare_fractions = mcintosh_mflare_12hr
    flare_fractions = (bornmann_shaw_mcintosh_mflare_12hr + bornmann_shaw_mcintosh_xflare_12hr) / 589.
    bins = regular_hist_bins(flare_fractions, 0.01)
    weights = np.ones_like(flare_fractions)/len(flare_fractions)
    plt.hist(flare_fractions, bins=bins, weights=weights,color='blue')
    plt.ylabel('probability', fontsize=fontsize)
    plt.xlabel('fraction of flares >M1 observed', fontsize=fontsize)
    plt.tick_params(labelsize=labelsize)
    mn = np.mean(flare_fractions)
    plt.axvline(mn, label=f'mean={mn:.3f}', **mean_line)
    percentiles = np.percentile(flare_fractions, [2.5, 97.5])
    plt.axvline(percentiles[0], **percentile_line)
    plt.axvline(percentiles[1], label=f'95% interval', **percentile_line) # is {percentiles[0]:.3f} to {percentiles[1]:.3f}', **percentile_line)
    plt.legend(fontsize=8)
    plt.grid(**grid_style)
    plt.xlim([0.3, 0.6])
    plt.ylim([0,0.3])
    plt.title('Bornmann and Shaw (1994), 12 hr delay',fontsize=8)
    #plt.text(0.6,0.025,'h')

    plt.savefig('mcintosh_bloomfield_comparison.pdf')
