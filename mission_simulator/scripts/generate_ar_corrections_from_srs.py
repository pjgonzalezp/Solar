"""
The purpose of this script is to compare ARs recorded in the HEK with the ARs recorded in NOAA's
SRS files.  It is found that the HEK is missing a few records compared to the SRS file.  This
script finds where the differences are and writes a CSV file with corrections to the HEK record.

The corrections to the HEK-based AR information can be incorporated in to the AR information used
in this project using the existing tools that read the HEK-based files.

"""
import os
import datetime
import matplotlib.pyplot as plt
from astropy.io.ascii import InconsistentTableError
from astropy.coordinates import SkyCoord

from sunpy.io.special import read_srs
from sunpy.coordinates import frames
import numpy as np
import pandas as pd
#import flare_mission_sim as fm
import __init__ as fm
#from flare_mission_sim import flare_mission_sim as ms
import flare_mission_sim as ms

log = fm.setup_logging()

# The location of NOAA's SRS files.
# A SRS file is a daily report on solar active regions.  They are stored in one year tar's gzipped
# files.  The are available for download from ftp://ftp.swpc.noaa.gov/pub/warehouse/
# Note that the rest of this script expects that the SRS files are extracted from the tar'd
# qzipped files and are stored in year-based subdirectories, for example
# <srs_directory>/2011/
srs_directory = os.path.expanduser('~/Data/SRS/')

# Filepath for the output CSV file containing the AR corrections
srs_save_filepath = os.path.join(fm.data_dir, 'srs_derived_ar_corrections.csv')

# Start and end dates to examine
start_date = datetime.datetime(2011, 1, 1)
end_date = datetime.datetime(2018, 12, 31)


def days(start_date, end_date):
    """A little practice routine, trying out coroutines in Python"""
    this_day = start_date
    while this_day <= end_date:
        yield this_day
        this_day = this_day + datetime.timedelta(days=1)


def ar_data_table_format(column_data):
    """Format data in the form used by the ar_data table"""
    n = len(column_data)
    if n == 0:
        output = 'NaN'
    elif n == 1:
        output = str(column_data[0])
    else:
        output = str(column_data[0])
        for i in range(1, n):
            output += ',' + str(column_data[i])
    return output.upper()


# Load the HEK data
ar_data = ms.load_ar_data()

# Get the dates
these_days = []
for day in days(start_date, end_date):
    # Store the days
    these_days.append(day)

# Get the number of active regions and the NOAA active region numbers as recorded by the HEK
hek_info = dict()
hek_info['n_I'] = []
hek_info['NOAA_I'] = []
for this_day in these_days:
    hek_record = ms.get_ars_for_day(str(this_day.date()), ar_data)
    hek_info['n_I'].append(hek_record[0])
    if hek_record[1] is None:
        hek_info['NOAA_I'].append(np.asarray([]))
    else:
        hek_info['NOAA_I'].append(np.asarray(sorted(hek_record[1])))


# Storage for the SRS information
srs_info = dict()
srs_info['n_I'] = []
srs_info['n_IA'] = []
srs_info['n_II'] = []
srs_info['NOAA_I'] = []
srs_info['NOAA_IA'] = []
srs_info['NOAA_II'] = []
srs_info['bad'] = []
srs_info['missing'] = []


# SRS keys that will be used in the
srs_keys = ['Number', 'Z', 'Mag Type', 'Number of Sunspots']
srs_data_storage_keys = ['Number', 'Z', 'Mag Type', 'Number of Sunspots', 'hpc_x', 'hpc_y']
srs_data = {key: [] for key in srs_data_storage_keys}

# Load in the AR information for this day from a SRS file, and keep the relevant data if there is
# a difference between the number of ARs recorded by the HEK and number recorded in the SRS files.
current_year = -1
srs_days = []
srs_n_I = []
for iday, day in enumerate(days(start_date, end_date)):
    # Get the SRS file for this day
    this_year = day.year
    directory = os.path.join(srs_directory, str(this_year))
    filename = day.strftime("%Y%m%d") + 'SRS.txt'
    filepath = os.path.join(directory, filename)

    if os.path.exists(filepath):
        try:
            srs = read_srs(filepath)
            # Get the AR types recorded by NOAA, count the number per type and get their AR numbers
            for ar_type in ('I', 'IA', 'II'):
                ar_type_index = srs['ID'] == ar_type
                srs_info['n_' + ar_type].append(np.sum(ar_type_index))
                srs_info['NOAA_' + ar_type].append(np.asarray(sorted([c for c in srs['Number'][ar_type_index]])))
        except InconsistentTableError or ValueError:
            srs_info['bad'].append(day)
            print('Bad = ' + filepath)
            for ar_type in ('I', 'IA', 'II'):
                srs_info['n_' + ar_type].append(-1)
                srs_info['NOAA_' + ar_type].append([])
    else:
        print('Not found = ' + filepath)
        srs_info['missing'].append(day)
        srs_info['n_I'].append(None)
        srs_info['n_IA'].append(None)
        srs_info['n_II'].append(None)

    n_difference_1 = srs_info['n_I'][-1] - hek_info['n_I'][iday]
    if n_difference_1 != 0:
        srs_days.append(day)
        n_I = srs_info['n_I'][-1]  # Number of ARs of the type we want to track
        srs_n_I.append(n_I)
        if n_I is not None:
            if n_I > 0:
                type_index_I = srs['ID'] == 'I'
                for key in srs_keys:
                    d = srs[key][type_index_I]
                    formatted_data = ar_data_table_format(d)
                    srs_data[key].append(formatted_data)

                # Convert the latitude and longitude to helioprojective Cartesian
                lat = srs["Latitude"][type_index_I]
                lon = srs["Longitude"][type_index_I]
                ar_location = SkyCoord(lon=lon,
                                       lat=lat,
                                       obstime=day,
                                       observer='earth',
                                       frame=frames.HeliographicStonyhurst).transform_to(frames.Helioprojective)

                srs_data['hpc_x'].append(ar_data_table_format(ar_location.Tx.value))
                srs_data['hpc_y'].append(ar_data_table_format(ar_location.Ty.value))


# Convert to a pandas dataframe - enables easy storage and inspection
srs_ar_data = pd.DataFrame(index=srs_days,
                           data={'noaa': srs_data['Number'],
                                 'number': srs_n_I,
                                 'hpc_x': srs_data['hpc_x'],
                                 'hpc_y': srs_data['hpc_y'],
                                 'classification': srs_data['Z'],
                                 'hale_classification': srs_data['Mag Type'],
                                 'numspots': srs_data['Number of Sunspots']})

# Save as a CSV file
srs_ar_data.to_csv(srs_save_filepath)


# Do the comparison
hek_consistency = np.zeros(shape=(len(these_days),))
srs_consistency = np.zeros_like(hek_consistency)
n_difference_1 = np.zeros_like(hek_consistency)
n_difference_2 = np.zeros_like(hek_consistency)
same_ar_numbers = np.zeros_like(hek_consistency)

for i, this_day in enumerate(these_days):
    # Is the HEK record consistent?
    hek_consistency[i] = hek_info['n_I'][i] - hek_info['NOAA_I'][i].size

    # Is the SRS record consistent?
    srs_consistency[i] = srs_info['n_I'][i] - srs_info['NOAA_I'][i].size

    # Are the number of active regions the same?
    n_difference_1[i] = hek_info['n_I'][i] - srs_info['n_I'][i]
    n_difference_2[i] = hek_info['NOAA_I'][i].size - srs_info['NOAA_I'][i].size

    # Look at the actual active region numbers
    hek_noaa = hek_info['NOAA_I'][i]
    srs_noaa = srs_info['NOAA_I'][i]

    # Different numbers of active regions on a day?
    if len(hek_noaa) != len(srs_noaa):
        same_ar_numbers[i] = -np.abs(len(hek_noaa) - len(srs_noaa))
    # Number of active regions is zero, so report zero
    elif len(hek_noaa) == 0:
        same_ar_numbers[i] = 0
    # Number of active regions is non zero and identical
    else:
        same_ar_numbers[i] = np.sum(hek_noaa == srs_noaa)

# Make some plots
plt.ion()

fig, ax = plt.subplots()
ax.plot(these_days, hek_info['n_I'], label='HEK')
ax.plot(these_days, srs_info['n_I'], label='SRS')
ax.plot(these_days, np.asarray(hek_info['n_I']) - np.asarray(srs_info['n_I']), label='HEK - SRS')
ax.set_ylabel('number of active regions recorded')
ax.set_title('daily comparison of number of ARs noted by SRS and HEK')
ax.legend()
plt.grid('on', linestyle=":")


fig, ax = plt.subplots()
ax.plot(these_days, same_ar_numbers, label='NOAA AR similarity')
ax.set_ylabel('NOAA AR number similarity for the HEK and SRS')
ax.set_title('similarity S')
ax.axhline(0, color='red', label='no ARs (both HEK and SRS agree)')
ax.legend()
plt.grid('on', linestyle=":")
