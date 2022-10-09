"""The purpose of this script is to create a flare database file"""
from sunpy.net import hek, attrs
import time
import datetime
import numpy as np
import h5py

client = hek.HEKClient()
current_year = datetime.datetime.now().year

start_year = 1996

number_of_years = current_year - start_year
years = np.arange(start_year, current_year+1, 1)

max_i = 0
i = 0

# create a local file
this_filename = f'goes_flares.hdf5'
print(f'Creating file {this_filename}.')
f = h5py.File(this_filename, 'w')


def get_hek_flares_for_year(year):
    start_time = datetime.datetime(year, 1, 1)
    end_time = datetime.datetime(year + 1, 1, 1)
    if start_time.year == current_year:
        end_time = datetime.datetime.now()

    flares_hek = client.search(attrs.Time(start_time, end_time),
                               hek.attrs.FL, hek.attrs.FRM.Name == 'SWPC')

    return flares_hek


def add_year_to_datafile(group_name, flares_hek, hdf5_handler):
    grp = hdf5_handler.create_group(group_name)
    num_flares = len(flares_hek)
    for k in flares_hek.keys():
        try:
            grp.create_dataset(k, (num_flares,), dtype="S20", data=np.array(flares_hek[k]).astype('S20'))
        except ValueError:
            print(f'skipping key {k}')
            pass


# create a file for each year
for this_year in years:

    print(f'Processing {this_year} ({i} of {len(years)})')
    if max_i != 0:
        if i >= max_i:
            break

    flares_hek = get_hek_flares_for_year(this_year)

    print(f'Adding {this_year} to the file.')
    add_year_to_datafile(str(this_year), flares_hek, f)

    i += 1
    # wait for a few seconds to not overload the hek server
    time.sleep(5)

f.close()

print('Remember to move file into the data directory after inspection.')
