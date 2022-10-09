"""
The purpose of this script is to create an AR database file.  The queried repository of these data
is the HEK.  The data obtained are the NOAA SWPC Observer records held by the HEK.
"""
from sunpy.net import hek, attrs
import time
import datetime
import numpy as np
import h5py
#import flare_mission_sim as fm
import __init__ as fm

client = hek.HEKClient()
current_year = datetime.datetime.now().year

start_year = 1996

number_of_years = current_year - start_year
years = np.arange(start_year, current_year+1, 1)

max_i = 50
i = 0

log = fm.setup_logging()

# create a local file
this_filename = f'swpc_ars.hdf5'

log.info(f'Creating file {this_filename}.')
f = h5py.File(this_filename, 'w')


def get_hek_ar_for_year(year):
    start_time = datetime.datetime(year, 1, 1)
    end_time = datetime.datetime(year + 1, 1, 1)

    if start_time.year == current_year:
        end_time = datetime.datetime.now()

    ar_hek = client.search(attrs.Time(start_time, end_time), hek.attrs.AR,
                           hek.attrs.FRM.Name == 'NOAA SWPC Observer')

    return ar_hek


def add_year_to_datafile(group_name, ar_hek, hdf5_handler):
    grp = hdf5_handler.create_group(group_name)
    num_ars = len(ar_hek)
    for k in ar_hek.keys():
        try:
            grp.create_dataset(k, (num_ars,), dtype="S20",
                               data=np.array(ar_hek[k]).astype('S20'))
        except ValueError:
            log.info(f'skipping key {k}')
            pass


for this_year in years:

    log.info(f'Processing {years.min(), years.max()}, current year {this_year} ({i} of {len(years)})')
    if max_i != 0:
        if i >= max_i:
            break

    ar_hek = get_hek_ar_for_year(this_year)

    log.info(f'Adding year({this_year}) to the file.')
    add_year_to_datafile(str(this_year), ar_hek, f)

    i += 1
    # wait for a few seconds to not overload the hek server
    time.sleep(1)

f.close()
log.info('Remember to move this file into the data directory after inspection.')
