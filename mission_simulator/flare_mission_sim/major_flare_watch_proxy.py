from flare_mission_sim.util import mcintosh_class_to_rank, flareclass_to_flux
import numpy as np
import datetime
import pandas as pd


def generate_mfw_proxy(flare_data, ar_data):
    mfw_starts = []
    mfw_ends = []
    noaas = []
    x_coords = []
    y_coords = []
    for id, flare in flare_data.iterrows():
        #  if a > M5 flare occurs, begin a major flare watch
        if flareclass_to_flux(flare['goes_class']) > flareclass_to_flux('M5'):
            mfw_start = id
            noaanum = flare['noaa']
            hpc_x = flare['hpc_x']
            hpc_y = flare['hpc_y']
            time_indexes, classifications, rank_value = extract_ar_class(str(noaanum), ar_data)
            has_enough_complexity = np.asarray(rank_value) < 3.0
            print(has_enough_complexity)
            #find the last date after the flare where the AR has enough complexity
            #this is the end of the MFW
            if (has_enough_complexity.any() == False): #or (np.asarray(time_indexes)[has_enough_complexity][-1] < id):
                mfw_end = id + datetime.timedelta(1)
            elif (np.asarray(time_indexes)[has_enough_complexity][-1] < id):
                mfw_end = id + datetime.timedelta(1)
            else:
                mfw_end = np.asarray(time_indexes)[has_enough_complexity][-1]
            
            mfw_starts.append(mfw_start)
            mfw_ends.append(mfw_end)
            noaas.append(noaanum)
            x_coords.append(hpc_x)
            y_coords.append(hpc_y)
    #last step is to merge overlapping entries
    mfw_series = pd.DataFrame({'mfw_ends': mfw_ends, 'noaa': noaas, 'hpc_x': x_coords, 'hpc_y':y_coords}, index = mfw_starts)
  
    mfw_series_no_duplicates = drop_mfw_duplicates2(mfw_series)
    mfw_series_no_duplicates.to_csv('mfw_proxy.csv')
    return mfw_series_no_duplicates #  _no_duplicates


def drop_mfw_duplicates2(mfw_series):
    # now, search through new series and identify overlapping intervals
    drop_indices = []
    for i in range(0,len(mfw_series)):
        print(i)
        if (i > 0):
                #  check if interval i is wholly contained in interval i-1
                if (mfw_series.index[i] <= mfw_series['mfw_ends'][i - 1]) and (
                mfw_series['mfw_ends'][i] <= mfw_series['mfw_ends'][i - 1]):
                        drop_indices.append(i)
                # check if interval i is partially contained in interval i-1. Adjust interval i-1 and drop interval i
                elif (mfw_series.index[i] <= mfw_series['mfw_ends'][i - 1]) and (
                mfw_series['mfw_ends'][i] >= mfw_series['mfw_ends'][i - 1]):
                        mfw_series['mfw_ends'][i - 1] = mfw_series['mfw_ends'][i]
                        drop_indices.append(i)

    mfw_series2 = mfw_series.drop(mfw_series.index[drop_indices], inplace=False)
    return mfw_series2

    
def drop_mfw_duplicates(mfw_series):
    # now, search through new series and identify overlapping intervals
    drop_indices = []
    for i, t in enumerate(mfw_series):
        # check if interval i is wholly contained in interval i-1
        if (i > 0) and (mfw_series.index[i] <= mfw_series[i - 1]) and (
            mfw_series[i] <= mfw_series[i - 1]):
            drop_indices.append(i)
        # check if interval i is partially contained in interval i-1. Adjust interval i-1 and drop interval i
        elif (i > 0) and (mfw_series.index[i] <= mfw_series[i - 1]) and (
            mfw_series[i] >= mfw_series[i - 1]):
            mfw_series[i - 1] = mfw_series[i]
            drop_indices.append(i)

    mfw_series.drop(mfw_series.index[drop_indices], inplace=True)
    return mfw_series


def extract_ar_class(ar_noaanum, ar_data):

    
    new_ar_data = ar_data[ar_data['noaa'].str.contains(ar_noaanum)]

    time_indexes = []
    classifications = []
    rank_value = []
    
    for id, ar_entry in new_ar_data.iterrows():
        noaa_nums = str.split(ar_entry['noaa'],',')
        ind = noaa_nums.index(ar_noaanum)
        ar_class = str.split(ar_entry['classification'],',')[ind]
        classifications.append(ar_class)
        time_indexes.append(id)
        rank_value.append(mcintosh_class_to_rank(ar_class))


    return time_indexes, classifications, rank_value
        

    



