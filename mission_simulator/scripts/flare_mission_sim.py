import csv
import pickle
import collections
import os.path
import h5py
import pandas

import matplotlib.pyplot as plt
from matplotlib import dates

from util import *
import __init__ as fm

log = fm.setup_logging()
common_datetime_fmt = "%Y/%m/%d %H:%M"

# The length of a solar cycle
# equals floor(11.5 * 365.25), using a non int here breaks everything!
days_per_year = 365.25
years_per_solar_cycle = 11.5
SOLAR_CYCLE_DT = datetime.timedelta(days=int(np.floor(years_per_solar_cycle*days_per_year)))


def load_major_flare_watches(filename='major_flare_watches.csv', tshift=True):
    """
    Load the major flare watch file.

    Parameter
    =========
    filename: str
        The filename to load the data from.

    Returns
    =======
    data: pandas.DataFrame
    """
    url = os.path.join(fm.data_dir, filename)
    data = pd.read_csv(url, index_col=0, parse_dates=True,
                       names=['end_time', 'noaa', 'hpc_x', 'hpc_y'],
                       skiprows=1, dtype={'noaa': np.uint16})
    data['noaa'] = data['noaa'].astype(str).astype(int)
    data['end_time'] = pd.to_datetime(data['end_time'])
    # shift into the future
    if tshift:
        data.index = data.index + SOLAR_CYCLE_DT
        data['end_time'] = data['end_time'] + SOLAR_CYCLE_DT
    return data


def load_ar_data_from_hdf5(ar_data_hdf5_filename='swpc_ars.hdf5', srs_corrections_filename='srs_derived_ar_corrections.csv'):

    date_rows = []
    noaa_rows = []
    hpc_x_rows = []
    hpc_y_rows = []
    class_rows = []
    number_rows = []
    haleclass_rows = []
    numspot_rows = []

    datafile_url = os.path.join(fm.data_dir, ar_data_hdf5_filename)
    log.info('Loading AR data from ' + datafile_url)

    result = h5py.File(datafile_url, 'r')
    years = list(result.keys())

    for y in years:
        log.info('Processing AR data from ' + str(y))
        results_this_year = result[y]
        unique_times = np.unique(results_this_year['event_starttime'])

        for i, t in enumerate(unique_times):
            key = datetime.datetime.strptime(t.astype('str'),
                                             "%Y-%m-%dT%H:%M:%S")
            ind = np.where(results_this_year['event_starttime'][:] == t)[0]

            if len(ind) > 0:
                noaa_str = ''
                classification_str = ''
                hale_classification_str = ''
                hpc_x_str = ''
                hpc_y_str = ''
                comma_str = ''
                numspot_str = ''
                noaa_list = []

                for loc, j in enumerate(ind):
                    noaa_ar = results_this_year['ar_noaanum'][j].astype('str')
                    hpc_x = results_this_year['hpc_x'][j].astype('str')
                    hpc_y = results_this_year['hpc_y'][j].astype('str')
                    arcls = results_this_year['ar_mcintoshcls'][j].astype('str')
                    halecls = results_this_year['ar_mtwilsoncls'][j].astype('str')
                    numspot = results_this_year['ar_numspots'][j].astype('str')
                    if loc > 0:
                        comma_str = ','
                    if noaa_ar not in noaa_list:
                        noaa_list.append(noaa_ar)
                        noaa_str += comma_str + str(noaa_ar)
                        hpc_x_str += comma_str + str(hpc_x)
                        hpc_y_str += comma_str + str(hpc_y)
                        classification_str += comma_str + str(arcls)
                        hale_classification_str += comma_str + str(halecls)
                        numspot_str += comma_str + str(numspot)
                number_rows.append(len(noaa_str.split(',')))
                date_rows.append(key)
                noaa_rows.append(noaa_str)
                hpc_x_rows.append(hpc_x_str)
                hpc_y_rows.append(hpc_y_str)
                class_rows.append(classification_str)
                haleclass_rows.append(hale_classification_str)
                numspot_rows.append(numspot_str)

    data = pd.DataFrame(index=date_rows,
                        data={'noaa': noaa_rows, 'number': number_rows,
                              'hpc_x': hpc_x_rows, 'hpc_y': hpc_y_rows,
                              'classification': class_rows, 'hale_classification': haleclass_rows,
                              'numspots': numspot_rows})

    # Read and apply the SRS-based AR corrections to the HEK
    if isinstance(srs_corrections_filename, str):
        corrections_url = os.path.join(fm.data_dir, srs_corrections_filename)
        log.info('Loading SRS-derived AR corrections from ' + corrections_url)
        srs_ar_corrections = pd.read_csv(corrections_url)
        data.update(srs_ar_corrections)

    return data


def load_ar_data(filename='swpc_ars3.csv'):
    datafile_url = os.path.join(fm.data_dir, filename)
    log.info('Loading AR data from ' + datafile_url)
    data = pd.read_csv(datafile_url, index_col=0, parse_dates=True)
    data.sort_index(inplace=True)
    #  drop duplicate indices
    data = data[~data.index.duplicated(keep='first')]
    return data


def load_flare_data_from_hdf5(filename='goes_flares.hdf5'):
    """
    Load the raw HEK flare data.

    Parameter
    =========
    None

    Return
    ======
    data: pandas.DataFrame
    """

    datafile_url = os.path.join(fm.data_dir, filename)
    log.info('Loading flare data from ' + datafile_url)

    f = h5py.File(datafile_url, 'r')

    hpc_x = []
    hpc_y = []
    goes_class = []
    start_time = []
    peak_time = []
    end_time = []
    noaa = []

    for this_year in f:
        grp = f[this_year]

        hpc_x.extend(list(f[this_year]['hpc_x'][...].astype(float)))
        hpc_y.extend(list(f[this_year]['hpc_y'][...].astype(float)))

        goes_class.extend(list(f[this_year]['fl_goescls'][...].astype(str)))
        noaa.extend(list(f[this_year]['ar_noaanum'][...].astype(int)))
        peak_time.extend(list(f[this_year]['event_peaktime'][...].astype(str)))
        end_time.extend(list(f[this_year]['event_endtime'][...].astype(str)))
        start_time.extend(list(f[this_year]['event_starttime'][...].astype(str)))

    start_time = [pd.to_datetime(this_time.replace('T', ' ')) for this_time in start_time]
    end_time = [pd.to_datetime(this_time.replace('T', ' ')) for this_time in end_time]
    peak_time = [pd.to_datetime(this_time.replace('T', ' ')) for this_time in peak_time]

    data = pd.DataFrame(index=start_time,
                        data={'peak_time': peak_time, 'end_time': end_time,
                              'goes_class': goes_class,
                              'hpc_x': hpc_x, 'hpc_y': hpc_y, 'noaa': noaa})
    f.close()
    return data


def load_ar_flare_data():
    """
    Load the raw data from the pickle files.

    Returns
    =======
        ar_data, flare_data, mclass_flare_data, num_ars, dates
    """
    # load in the solar cycle data from files
    data_files = ['solar_cycle_data_2010_py3.pickle', 'solar_cycle_data_2011_py3.pickle',
                  'solar_cycle_data_2012_py3.pickle', 'solar_cycle_data_2013_py3.pickle']
    fhs = []
    for this_file in data_files:
        fhs.append(open(this_file, 'rb'))

    num_ar_2010, date_2010, data_2010, flare_data_2010, mclass_flare_data_2010 = pickle.load(fhs[0])
    num_ar_2011, date_2011, data_2011, flare_data_2011, mclass_flare_data_2011 = pickle.load(fhs[1])
    num_ar_2012, date_2012, data_2012, flare_data_2012, mclass_flare_data_2012 = pickle.load(fhs[2])
    num_ar_2013, date_2013, data_2013, flare_data_2013, mclass_flare_data_2013 = pickle.load(fhs[3])

    #  recombine all the solar cycle and flare data
    dates = date_2010 + date_2011 + date_2012 + date_2013
    num_ars = num_ar_2010 + num_ar_2011 + num_ar_2012 + num_ar_2013
    ar_data = dict(data_2010)
    ar_data.update(data_2011)
    ar_data.update(data_2012)
    ar_data.update(data_2013)
    flare_data = dict(flare_data_2010)
    flare_data.update(flare_data_2011)
    flare_data.update(flare_data_2012)
    flare_data.update(flare_data_2013)
    mclass_flare_data = dict(mclass_flare_data_2010)
    mclass_flare_data.update(mclass_flare_data_2011)
    mclass_flare_data.update(mclass_flare_data_2012)
    mclass_flare_data.update(mclass_flare_data_2013)

    #  ar_data is not sorted by date! Need to sort first!!
    ar_data = collections.OrderedDict(sorted(ar_data.items()))

    # please remember to close the file!
    for this_fh in fhs:
        this_fh.close()

    return ar_data, flare_data, mclass_flare_data, num_ars, dates


def _process_ar_data(ar_data, tshift=True):
    """
    OLD!
    Process the active region data into pandas data frame.

    Parameters
    ==========
    ar_data: dict
        As produced by load_ar_flare_data()

    tshift: bool (default True)
        If set to true than shift the flares forward by SOLAR_CYCLE_DT.

    Returns
    =======
    ar_series: pandas.DataFrame

    """
    date_rows = []
    noaa_rows = []
    hpc_x_rows = []
    hpc_y_rows = []
    class_rows = []
    number_rows = []
    for index, row in ar_data.iterrows():
        this_day = datetime.datetime.strptime(index, "%Y-%m-%d")
        ars_list = row['noaa']

        if row['number'] > 0:
            noaa_str = ''
            classification_str = ''
            hpc_x_str = ''
            hpc_y_str = ''
            comma_str = ''
            noaa_list = []
            for i, this_ar in enumerate(ars_list):
                if i > 0:
                    comma_str = ','
                if not this_ar['ar_noaanum'] in noaa_list:  # do not add duplicates
                    noaa_list.append(this_ar['ar_noaanum'])
                    noaa_str += comma_str + str(this_ar['ar_noaanum'])
                    hpc_x_str += comma_str + str(this_ar['hpc_x'])
                    hpc_y_str += comma_str + str(this_ar['hpc_y'])
                    classification_str += comma_str + this_ar['ar_mcintoshcls']
            number_rows.append(len(noaa_str.split(',')))
            date_rows.append(this_day)
            noaa_rows.append(noaa_str)
            hpc_x_rows.append(hpc_x_str)
            hpc_y_rows.append(hpc_y_str)
            class_rows.append(classification_str)

    data = pd.DataFrame(index=date_rows,
                        data={'noaa': noaa_rows, 'number': number_rows,
                              'hpc_x': hpc_x_rows, 'hpc_y': hpc_y_rows,
                              'classification': class_rows})
    # shift into the future
    if tshift:
        data.index = data.index + SOLAR_CYCLE_DT
    return data


def shift_by_solar_cycle(data, n=1):
    """Shift the data by n solar cycles.

    Parameters
    ==========
    n: int (default 1)
    """
    data.index = data.index + n * SOLAR_CYCLE_DT
    if isinstance(data, pd.DataFrame):
        if 'peak_time' in data.columns:
            data['peak_time'] += n * SOLAR_CYCLE_DT
        if 'end_time' in data.columns:
            data['end_time'] += n * SOLAR_CYCLE_DT
    elif (isinstance(data, pd.Series)) and (isinstance(data.values[0], np.datetime64)):
        data += n * SOLAR_CYCLE_DT
    return data


def _process_flare_data(flare_data, mclass_flare_data, tshift=True):
    """Process the flare list data into pandas data frame.

    Parameters
    ==========
    tshift: bool (default True)
        If set to true than shift the flares forward by SOLAR_CYCLE_DT.
    """
    start_times = []
    peak_times = []
    end_times = []
    noaa = []
    goes_class = []
    hpc_x = []
    hpc_y = []

    if tshift:
        my_time_shift = SOLAR_CYCLE_DT
    else:
        my_time_shift = datetime.timedelta(days=0)

    # process x class flares
    for flare_list in [flare_data, mclass_flare_data]:
        for key in flare_list:
            if len(flare_list[key]) != 0:
                this_day = datetime.datetime.strptime(key, "%Y-%m-%d %H:%M:%S")
                for this_flare in flare_list[key]:
                    start_times.append(datetime.datetime.strptime(this_flare['event_starttime'], "%Y-%m-%dT%H:%M:%S") + my_time_shift)
                    peak_times.append(datetime.datetime.strptime(this_flare['event_peaktime'], "%Y-%m-%dT%H:%M:%S") + my_time_shift)
                    end_times.append(datetime.datetime.strptime(this_flare['event_endtime'], "%Y-%m-%dT%H:%M:%S") + my_time_shift)
                    hpc_x.append(this_flare['hpc_x'])
                    hpc_y.append(this_flare['hpc_y'])
                    goes_class.append(this_flare['fl_goescls'])
                    noaa.append(this_flare['ar_noaanum'])

    data = pd.DataFrame(index=start_times,
                        data={'peak_time': peak_times, 'end_time': end_times,
                              'goes_class': goes_class, 'hpc_x': hpc_x,
                              'hpc_y': hpc_y, 'noaa': noaa})
    data.sort_index(inplace=True)
    data.drop_duplicates(inplace=True)

    return data


def fix_flare_data(flare_series, ar_series):
    """Replace locations where the flare positions are likely bad and replace
    with location of the source active region.

    Parameters
    ==========
    """

    bad_flare_series = flare_series[flare_series['hpc_x'] == 0]
    log.info(f"Found {len(bad_flare_series)} likely bad positions. Fixing.")
    for this_time, this_row in bad_flare_series.iterrows():
        source_noaa = this_row['noaa']
        if source_noaa != 0:
            ar_pos = get_ar_positions_over_time(source_noaa, ar_series)
            this_date = this_time.to_pydatetime().date()
            if this_date in ar_pos.index:
                flare_series.at[this_time, 'hpc_x'] = ar_pos['hpc_x'][this_date]
                flare_series.at[this_time, 'hpc_y'] = ar_pos['hpc_y'][this_date]
            else:
                flare_series.at[this_time, 'hpc_x'] = 0.0
                flare_series.at[this_time, 'hpc_y'] = 0.0

    #  now do some hand fixing, the following positions were found through the
    #  RHESSI browser plots except for
    # 2012/03/07 00:02:00 which was found in Patsourakos et al 2016 as N18, E31
    time = ['2011/02/15 01:44:00', '2011/09/22 10:29:00', '2012/03/07 00:02:00',
            '2012/07/06 23:01:00', '2012/10/23 03:13:00', '2013/05/13 01:53:00',
            '2013/05/14 00:00:00', '2013/10/25 07:53:00', '2013/10/25 14:51:00',
            '2013/10/29 21:42:00', '2013/11/19 10:14:00']
    position = [[196, -231], [-967, 192], [-417, 297],
                [708, -307], [-798, -273], [-946, 177],
                [-928, 204], [-917, -161], [-896, -167],
                [998, 59], [906, -252]]
    for this_time, this_position in zip(time, position):
        shift_time = datetime.datetime.strptime(this_time,
                                                common_datetime_fmt + ':%S')
        flare_series.at[shift_time, 'hpc_x'] = this_position[0]
        flare_series.at[shift_time, 'hpc_y'] = this_position[1]

    bad_flare_series = flare_series[flare_series['hpc_x'] == 0]
    log.info("{0} likely bad positions left ({1} X class flares {2} M class flares.".format(len(bad_flare_series),
                                                                                            np.sum(bad_flare_series['goes_class'].str.contains('X')),
                                                                                            np.sum(bad_flare_series['goes_class'].str.contains('M'))))
    return flare_series


def fix_flare_positions_sff(flare_data):

    sff_reference_flare_data = pandas.read_csv(os.path.join(fm.data_dir, 'ssw_sff_list.txt'),header=None)
    #need to convert the SFF reference data time column format to do matching
    sff_starttimes = []
    for t in sff_reference_flare_data[1].values:
        sff_starttimes.append(datetime.datetime.fromisoformat(t))

    sff_starttimes = np.asarray(sff_starttimes)
    sff_reference_flare_data[1] = sff_starttimes


    for ind, flare in flare_data.iterrows():
        #find the matching flare in the SFF database
        reference_flare = sff_reference_flare_data[sff_reference_flare_data[1] == flare.name]
        #print(len(reference_flare))
        if (not reference_flare.empty) and (len(reference_flare) == 1):
            #replace with the new flare position from SFF database
            flare_data.at[ind,'hpc_x'] = float(reference_flare[5])
            flare_data.at[ind,'hpc_y'] = float(reference_flare[6])



    return flare_data

def fix_flare_positions_manual(flare_data):
    """
    Use a separate file (`flare_manually_corrected_positions.csv`) to correct
    the positions in the flare list with manually defined flare positions.

    As a special case, if the manually defined flare position is (999999,
    999999), the flare will be *removed* from the flare list.
    """
    manual_flare_data = pandas.read_csv(os.path.join(fm.data_dir, 'flare_manually_corrected_positions.csv'))
    #need to convert the manual reference data time column format to do matching
    manual_starttimes = []
    for t in manual_flare_data['start_time'].values: 
        manual_starttimes.append(datetime.datetime.fromisoformat(t))

    manual_starttimes = np.asarray(manual_starttimes)
    manual_flare_data['start_time'] = manual_starttimes

    flares_to_remove = []

    for ind, flare in flare_data.iterrows():
        #find the matching flare in the manual database
        reference_flare = manual_flare_data[manual_flare_data['start_time'] == flare.name]
        if (not reference_flare.empty) and (len(reference_flare) == 1):
            if (int(reference_flare['hpc_x']) == 999999) and (int(reference_flare['hpc_y']) == 999999):
                #remove the flare if HPC (x, y) are (999999, 999999)
                flares_to_remove.append(ind)
            else:
                #replace with the new flare position from manual database
                flare_data.at[ind,'hpc_x'] = float(reference_flare['hpc_x'])
                flare_data.at[ind,'hpc_y'] = float(reference_flare['hpc_y'])

    if len(flares_to_remove) > 0:
        flare_data = flare_data.drop(flares_to_remove)

    return flare_data
    

def _ardata_fromhdf5_to_csv(output_filename_csv, ar_data_hdf5_filename='swpc_ars.hdf5',
                            srs_corrections_filename='srs_derived_ar_corrections.csv'):
    """
    Reads the HDF5 file of AR data and a AR data correction file, and writes out the AR data in
    CSV format.
    """
    ar_data = fm.load_ar_data_from_hdf5(ar_data_hdf5_filename=ar_data_hdf5_filename,
                                        srs_corrections_filename=srs_corrections_filename)
    out_path = os.path.join(fm.data_dir, output_filename_csv)
    log.info("Writing AR data to CSV format at " + out_path)
    ar_data.to_csv(out_path)


def read_orbit_events(filename=None, pad_saa=0, pad_eclipse=0, pad_polar=0):
    """Read the orbit events file and return Eclipse, SAA, and contact times.

    Parameters
    ==========
    filename: str (optional)
        the file to be read. If no file is given then defaults to one defined
        by mission_parameters.csv
    pad_polar: int (number of seconds)
    pad_eclipse: int (number of seconds)
    pad_saa: int (number of seconds)
    """
    if filename is None:
        url = os.path.join(fm.data_dir, fm.FILE_ORBIT_EVENTS)

    else:
        url = filename

    fh = open(url, encoding='utf-8')
    result = csv.DictReader(fh, fieldnames=['time', 'event_type'])
    saa_entries = []
    saa_exits = []
    eclipse_entries = []
    eclipse_exits = []
    #  observing_status = []
    ground_contact = []
    ground_contact_station = []
    polar_entries = []
    polar_exits = []

    eclipse_entry_strings = ('Eclipse Entry - Earth', 'Penumbra entry')
    eclipse_exit_strings = ('Eclipse Exit - Earth', 'Penumbra exit')

    polar_entries_strings = ('South Polar Region Entry', 'North Polar Region Entry')
    polar_exits_strings = ('South Polar Region Exit', 'North Polar Region Exit')

    for i, row in enumerate(result):
        # read in each row from the orbit file and format appropriately
        # ignore the first few entries for simplicity
        if i > 3:
            dat, tim = row['time'].split(' ')
            dat2 = dat.split('/')
            tim2 = tim.split(':')
            # reformat the time string into a datetime object
            actual_time = pd.to_datetime(row['time'])
            #actual_time = datetime.datetime(int('20' + dat2[2]), int(dat2[0]),
            #                                int(dat2[1]), int(tim2[0]),
            #                                int(tim2[1]), int(tim2[2]))
            # actual_time = datetime.datetime.strptime(row['time'], '%-m/%-d/%y %-H:%M:%S')
            row['time'] = actual_time

            if row['event_type'] == 'SAA Entry':
                saa_entries.append(actual_time)
            elif (row['event_type'] == 'SAA Exit') and len(saa_entries) > 0:  # don't record an exit unless there is an entry
                saa_exits.append(actual_time + datetime.timedelta(0, pad_saa))
            elif any(x in row['event_type'] for x in eclipse_entry_strings):
                eclipse_entries.append(actual_time - datetime.timedelta(0, pad_eclipse))
            elif any(x in row['event_type'] for x in eclipse_exit_strings):
                if len(eclipse_entries) > 0:   # don't record an exit unless there is an entry
                    eclipse_exits.append(actual_time + datetime.timedelta(0, pad_eclipse))
            elif any(x in row['event_type'] for x in ('Ground', 'AOS')):
                ground_contact.append(actual_time)
                if 'Ground' in row['event_type']:
                    ground_contact_station.append(row['event_type'].split('-')[1][1:])
                if 'AOS' in row['event_type']:
                    ground_contact_station.append(row['event_type'].replace('AOS', '').replace(' ', ''))
            elif any(x in row['event_type'] for x in polar_entries_strings):
                polar_entries.append(actual_time - datetime.timedelta(0, pad_polar))
            elif any(x in row['event_type'] for x in polar_exits_strings):
                if len(polar_entries) > 0:  # don't record an exit unless there is an entry
                    polar_exits.append(actual_time + datetime.timedelta(0, pad_polar))

    # there must be as many entries as exits
    saa = pd.Series(data=saa_exits, index=saa_entries[0:len(saa_exits)])

    # there must be as many entries as exits
    eclipse = pd.Series(data=eclipse_exits, index=eclipse_entries[0:len(eclipse_exits)])
    contacts = pd.Series(data=ground_contact_station, index=ground_contact)

    if len(polar_entries) > 0:
        polar = pd.Series(data=polar_exits, index=polar_entries[0:len(polar_exits)])
    else:
        polar = None

    fh.close()
    return eclipse, saa, contacts, polar


def shuffle_orbit_events(eclipse, saa, contacts):
    """Shuffle the default input eclipse, saa and contact times for monte-carlo
    purposes"""
    offset = datetime.timedelta(0, int(np.random.uniform(low=0, high=86400)))
    eclipse2 = eclipse + offset
    eclipse2.index = eclipse.index + offset
    saa2 = saa + offset
    saa2.index = saa.index + offset
    contacts2 = copy.deepcopy(contacts)
    contacts2.index = contacts.index + offset

    return eclipse2, saa2, contacts2


def load_goes_fitsfile_index(filename='flare_goes_fits_filenames.csv'):
    """Load the index files which associates goes fits files with flares."""
    url = os.path.join(fm.data_dir, filename)
    return pd.read_csv(url, parse_dates=True, index_col=0)


def get_flare_timeseries(flare_series, skip_index=None):
    """Download goes fits files for each flare."""
    import time
    from sunpy.net import Fido, attrs as a
    from sunpy.time import TimeRange
    files = []
    times = []
    goes_class = []
    # download the time profile for each x class flare
    i = 0
    for start_time, this_row in flare_series.iterrows():
        tr = TimeRange([start_time, this_row['end_time']])
        tr = TimeRange(tr.start - SOLAR_CYCLE_DT, tr.end - SOLAR_CYCLE_DT)
        log.info(f'{i} {this_row["goes_class"]}: {tr.start_time}-{tr.end_time}')
        if not (i in skip_index):
            results = Fido.search(a.Time(tr), a.Instrument('XRS'))
            these_files = Fido.fetch(results, wait=True, progress=True)
            if these_files is not None:
                times.append(start_time)
                goes_class.append(this_row['goes_class'])
                files.append(os.path.basename(these_files[0]))
        else:
            log.info('Skipping')
        i += 1
        time.sleep(1)
    data = pd.DataFrame(index=times, data={'goes_class': goes_class,
                                           'filename': files})
    return data


def get_ars_for_day(time, ar_series):
    """Given a date return all of the active regions present on that date."""
    this_day = str(pd.to_datetime(time).to_pydatetime().date())
    # check if data exists first
    try:
        these_ars = [int(this_noaa) for this_noaa in ar_series['noaa'][this_day].split(',')]
        these_hpc_xs = [float(this_x) for this_x in ar_series['hpc_x'][this_day].split(',')]
        these_hpc_ys = [float(this_y) for this_y in ar_series['hpc_y'][this_day].split(',')]
        these_classs = ar_series['classification'][this_day].split(',')
        these_hale_classes = ar_series['hale_classification'][this_day].split(',')
        these_numspots = ar_series['numspots'][this_day].split(',')
        the_number = ar_series['number'][this_day]
    except KeyError:
        #print(f'{this_day} no ARs found.')
        return 0, None, None, None, None, None, None
    except AttributeError:
        return 0, None, None, None, None, None, None
    return the_number, np.array(these_ars), np.array(these_hpc_xs), np.array(these_hpc_ys), these_classs, these_hale_classes, these_numspots


def get_ar_positions_over_time(noaa, ar_series):
    """Given a NOAA index return all positions for that active region as a
    function of time."""
    index = ar_series['noaa'].str.contains(str(noaa))
    indices = []
    hpc_y = []
    hpc_x = []
    if index.any():
        dates = ar_series['noaa'][index].index
        # get the indices

        for this_row in ar_series['noaa'][index]:
            temp = np.array(this_row.split(','))
            indices.append(np.where(temp == str(noaa))[0][0])
        for i, this_row in enumerate(ar_series['hpc_x'][index]):
            temp = np.array(this_row.split(','))
            hpc_x.append(float(temp[indices[i]]))
        for i, this_row in enumerate(ar_series['hpc_y'][index]):
            temp = np.array(this_row.split(','))
            hpc_y.append(float(temp[indices[i]]))
    else:
        dates = [ar_series['noaa'].index[0]]
        hpc_x = hpc_y = None
    return pd.DataFrame(index=dates, data={'hpc_x': hpc_x, 'hpc_y': hpc_y})


def get_flares(time, flare_series):
    """Given a date return all of the flares that occurred on that date."""
    try:
        these_goes_classs = list(flare_series['goes_class'][time].values)
        these_hpc_xs = [float(this_x) for this_x in flare_series['hpc_x'][time]]
        these_hpc_ys = [float(this_y) for this_y in flare_series['hpc_y'][time]]
        the_number = len(flare_series['hpc_x'][time])
    except KeyError:
        the_number = 0
        these_hpc_xs = these_hpc_ys = these_goes_classs = None
    return the_number, these_goes_classs, these_hpc_xs, these_hpc_ys


def is_sun_observed(time, dilate=False, verbose=False, eclipse=None, saa=None):
    """Given a time return True if the Sun is observable meaning is the
    Observatory outside of both eclipse and SAA.

    Parameters
    ----------
    dilate: the time to dilate the ranges.
    """
    time_dt = datetime.datetime.strptime(time.replace('-', '/')[0:16], common_datetime_fmt)

    if eclipse is not None:
        is_in_eclipse, eclipse_its_in = is_in_ranges(time_dt, eclipse)
        if is_in_eclipse and verbose:
            log.info(f"{time} is {time_dt - eclipse_its_in[0]} in eclipse [{eclipse_its_in[0], eclipse_its_in[1]}]")
    else:
        is_in_eclipse = False
    if saa is not None:
        is_in_saa, saa_its_in = is_in_ranges(time_dt, saa)
        if is_in_saa and verbose:
            log.info(
                "{0} is {dt} in saa [{1}, {2}]".format(time, saa_its_in[0], saa_its_in[1], dt=time_dt - saa_its_in[0]))
    else:
        is_in_saa = False
    return not (is_in_eclipse or is_in_saa)


def percent_is_sun_observed(start_time, end_time, eclipse_and_saa):
    """
    Given a start and end time provide the amount of that time range which is observable.

    If the provided eclipse/SAA times is ``None``, the entire range is considered observed.
    """
    # TODO: this function double counts non observables times when both in eclipse and saa at the same time.
    # SAA times are not very long so ignoring this therefore provides a lower bound on percent observable time
    if eclipse_and_saa is None:
        return 1

    if (end_time - start_time).total_seconds() == 0:  # can't calculate percent in this case
        return 0

    start_times = []
    end_times = []
    # print("Start time is {0} end time is {1}".format(start_time, end_time))
    # first check if start is in eclipse time
    is_in_eclipse_or_saa, this_range = is_in_ranges(start_time, eclipse_and_saa)

    if is_in_eclipse_or_saa:
        # print("start is in eclipse or saa {0}".format(this_range))
        start_times.append(start_time)
        end_times.append(this_range[1])

    # print(eclipse_and_saa[start_time:end_time])
    for this_start_time, this_end_time in eclipse_and_saa[start_time:end_time].iteritems():
        start_times.append(this_start_time)
        end_times.append(this_end_time)

    # check if end is in eclipse or saa time
    is_in_eclipse_or_saa, this_range = is_in_ranges(end_time, eclipse_and_saa)
    if is_in_eclipse_or_saa:
        # print("end is in eclipse or saa {0}".format(this_range))
        end_times[-1] = end_time

    non_observed_time_sec = 0.0

    #  now add up all of the times that could not observe
    if len(start_times) != 0:
        data = pd.Series(index=start_times, data=end_times)
        # print(data)
        for this_start_time, this_end_time in data.iteritems():
            non_observed_time_sec += (this_end_time - this_start_time).total_seconds()
    if (end_time - start_time).total_seconds() == 0:
        print(start_time)
        print(end_time)
    return 1 - non_observed_time_sec / (end_time - start_time).total_seconds()


def generate_ar_targets(ar_series, mfw=None, allowed_days_of_week = [1,2,3,4,5,6,7],
                        target_method='mcintosh',
                        rank_threshold=0, fallback_target=0):
    """
    Generate the daily time series of AR targets to simulate SOC decisions.

    There are several special AR numbers that are not "real":
    * 0: point to Sun center
    * -1: continue pointing to previous AR (valid only for `fallback_target`, will not appear in output)
    * -1xx: point to cover the edge of the solar disk
      * xx corresponds to numbers on a clock, so it goes from 01 to 12
      * For example, -110 points to the eastern limb of the northern latitude band (~10 o'clock)
    """

    pointing_location_decision = []
    decision_date = []
    ar_targets = []
    descriptions = []

    current_noaa_target = None
    in_major_flare_watch = False

    # generate some extra information for specific targeting methods only if needed
    if target_method == 'hale':
        allowed_labels = find_unique_hale_classes()
    elif target_method == 'flare_index':
        flare_index = generate_flare_index()

    for key, this_row in ar_series.iterrows():
        major_flare_watch_start = False
        # this_day = datetime.datetime.strptime(key, "%Y-%m-%d %H:%M:%S")
        this_day = key
        dayofweek = this_day.isoweekday()

        # if a major flare watch list was provided than use it
        if mfw is not None:
            # check if a major flare watch was initiated today
            try:
                mfw_today = mfw[this_day.strftime('%Y-%m-%d')]
                if len(mfw_today) > 0:
                    current_noaa_target = mfw_today['noaa'].values[0]
                    decision_date.append(mfw_today.index[0].to_pydatetime())
                    ar_targets.append(current_noaa_target)
                    descriptions.append('DECISION: pointing change, major flare watch')
                    major_flare_watch_start = True
            except KeyError:
                pass

            # check if a major flare watch is still on going, if so than stay on target
            in_major_flare_watch, flare_watch_range = is_in_ranges(this_day, mfw['end_time'], consider_range_days=20)
            if in_major_flare_watch:
                decision_date.append(this_day.date() + datetime.timedelta(hours=10))
                ar_targets.append(current_noaa_target)
                descriptions.append('DECISION: no change, major flare watch stay')

        if (not in_major_flare_watch) and (not major_flare_watch_start):  # no major flare watch so we are on our own to decide where to point
            this_date = this_day.date()
            # SOC decision time is 00:30 UT, coinciding with when new SWPC reports are published
            this_decision_time = datetime.datetime(this_date.year, this_day.month, this_date.day, 0, 30, 0)
            decision_date.append(this_decision_time)
            if dayofweek in allowed_days_of_week:  # make sure it is a day where repoints are allowed
                # get active region from yesterday
                # using ar_series.shift(1) to use yesterday's information on active region classification. - NO! AR report already reflects yesterday's data
                yesterday_num, yesterday_noaa, yesterday_hpc_x, yesterday_hpc_y, yesterday_mcintosh, yesterday_hale, yesterday_numspots = get_ars_for_day(this_day, ar_series.shift(0))
                if yesterday_num != 0:
                    if target_method == 'mcintosh':
                        ar_rank = [mcintosh_class_productivity(this_mcintosh) for this_mcintosh in yesterday_mcintosh]
                    elif target_method == 'mcintosh_bloomfield':
                        ar_rank = [mcintosh_class_productivity_bloomfield_2012(this_mcintosh) for this_mcintosh in yesterday_mcintosh]
                    elif target_method == 'mcintosh_bornmann_shaw':
                        ar_rank = [mcintosh_class_productivity_bloomfield_2012(this_mcintosh, use_bornmann_shaw = True) for this_mcintosh in yesterday_mcintosh]
                    elif target_method == 'old_mcintosh':
                        ar_rank = [mcintosh_class_to_rank(this_mcintosh) for this_mcintosh in yesterday_mcintosh]
                    elif target_method == 'hale':
                        ar_rank = [hale_class_productivity(this_hale, allowed_labels) for this_hale in yesterday_hale]
                        # break ties in Hale system by using number of sunspots as a tiebreaker
                        ar_rank = list(np.array(ar_rank) + np.array(yesterday_numspots, dtype='float'))
                    elif target_method == 'random':
                        ar_rank = [mcintosh_class_to_random_rank(this_mcintosh) for this_mcintosh in yesterday_mcintosh]
                    elif target_method == 'flare_index':
                        try:
                            # need to get flare index data from YESTERDAY'S flares, so -datetime.timedelta(1)
                            flare_index_entry = flare_index.loc[ (this_day.date() - datetime.timedelta(1)).isoformat()]
                            ar_rank = [noaanum_to_flare_index_rank(this_noaa, flare_index_entry) for this_noaa in yesterday_noaa]
                        except:
                            ar_rank = [0] * len(yesterday_noaa)
                        
                    else:
                        # raise an error if an invalid method is chosen
                        raise ValueError('Unrecognized target determination method. Options are mcintosh, old_mcintosh, hale, random.')

                    # pick the highest 'ranked' AR as the target
                    if np.max(ar_rank) >= rank_threshold:
                        best_noaa_target = yesterday_noaa[np.argmax(ar_rank)]
                    else:
                        best_noaa_target = fallback_target if fallback_target != -1 else current_noaa_target

                    # if the fallback target is "previous AR", but there is no previous AR, set to Sun center
                    if best_noaa_target is None:
                        best_noaa_target = 0

                    # make sure that yesterdays active region still exists today!
                    if best_noaa_target != current_noaa_target:
                        current_noaa_target = best_noaa_target
                        if current_noaa_target > 0:
                            descriptions.append(f'DECISION: pointing change, NEW best {target_method}')
                        else:
                            descriptions.append('DECISION: pointing change, reverting to fallback')
                    else:
                        descriptions.append('DECISION: no change')
                    ar_targets.append(current_noaa_target)
                else:
                    ar_targets.append(0)
                    descriptions.append('DECISION: no change, no active regions')
            else:  # if its a forbidden day to repoint, make no changes
                ar_targets.append(current_noaa_target)
                descriptions.append('DECISION: no change, repoints not allowed this day')
    result = pd.DataFrame(index=decision_date,
                          data={'ar_target': ar_targets, 'reason': descriptions})
    return result


def generate_pointing_commands(ar_targets, contacts, monte_carlo=False, delay_params=[4.0, 4.0, 2.0], use_random_delay = False):
    """Given a set of AR targets, generate the corresponding list of pointing commands for
    *changes* in the AR target, accounting for contact times/delays.

    Tracking of a target as it moves across the Sun is not part of this set of
    commands.
    """

    pointing_change_index = ar_targets['reason'].str.contains('pointing change')
    ar_target_changes = ar_targets[pointing_change_index]

    repoint_time = []
    ar_target = []
    station = []
    decision_time = []
    reason = []
    station_used = []

    for this_time, this_row in ar_target_changes.iterrows():
        check_mfw = 'major flare watch' in this_row['reason']
        if monte_carlo and check_mfw == False:
            shall_we_repoint = np.random.uniform(low=0,high=1)
        else:
            shall_we_repoint = 1.0

        # if desired, use a delay drawn from a distribution rather than specific ground pass information
        if use_random_delay:
            # delay_hours is the amount of time it takes for command to the spacecraft to be generated and validated by the MOC and SOC.
            # baseline constant delay PLUS a random extra delay drawn from a Gaussian. Default base 4 hours, mean 4 hours, std 2 hours.
            # i.e. typical delay is ~ 8 +/- 2 hours
            delay_hours = delay_params[0] + max(0.0, np.random.normal(loc = delay_params[1], scale = delay_params[2]))
            to_cmd = datetime.timedelta(hours=delay_hours)
            repoint_time.append(this_time + to_cmd)
            station_used.append('None')
            station.append('None')
        # otherwise, use ground station pass information to determine repoint delay
        else:
            # find the next scheduled contact which is between 3 and 5 pm ET, (~8-10pm UT) after the command is ready.
            today_str = str(this_time.date())
            available_contacts = contacts[today_str + ' 20:00': today_str + ' 22:00']
            if len(available_contacts) == 0: # or (shall_we_repoint < 0.0):
                # there are no contacts available today, get the first contact during working hours tomorrow.
                tomorrow_str = str((this_time + pd.Timedelta(days=1)).date())
                available_contacts = contacts[tomorrow_str + ' 14:00': tomorrow_str + ' 22:00']
            # TODO this assumes that there necessarily is a pass the next day which may not be true!
            contact_time, contact_station = [available_contacts.index[0], available_contacts[0]]
            # b = (contacts.index > cmd_ready_time).argmax()  # index of first entry after target
            # contact_time, contact_station = [contacts.index[b], contacts[b]]
            station_used.append(contact_station)
            repoint_time.append(contact_time)
            station.append(contact_station)
            
        ar_target.append(this_row['ar_target'])
        reason.append(this_row['reason'])
        decision_time.append(this_time)

    data = {'ar_target': ar_target, 'reason': reason,
            'decision_time': decision_time, 'station': station_used}
    result = pd.DataFrame(index=repoint_time, data=data)
    return result


def get_pointing(time, pointing_commands, ar_series):
    """
    Given a time and the pointing command list, return the position of pointing
    in HPC.

    Currently performs a linear interpolation between the daily locations of
    the active region to simulate tracking.  Ideally should not know the
    future and instead track based on known latitude-dependent rotation rates.

    There are several special AR numbers that are not "real":
    * 0: point to Sun center
    * -1xx: point to cover the edge of the solar disk
      * xx corresponds to numbers on a clock, so it goes from 01 to 12
      * For example, -110 points to the eastern limb of the northern latitude band (~10 o'clock)
    """
    # this_time = str(datetime.datetime.strptime(time, common_datetime_fmt))
    this_time = time
    ar_target = int(extract_series_at_time(pointing_commands['ar_target'], this_time))

    if ar_target > 0:  # real NOAA AR number
        ar_pos = get_ar_positions_over_time(ar_target, ar_series)
        # If the active region is no longer on the disk, the interpolated value is the last known value
        target_hpc_x = extract_series_at_time(ar_pos['hpc_x'], this_time, interpolate=True)
        target_hpc_y = extract_series_at_time(ar_pos['hpc_y'], this_time, interpolate=True)
    elif ar_target == 0:  # point to Sun center
        target_hpc_x = 0
        target_hpc_y = 0
    elif (ar_target >= -112) and (ar_target <= -101):
        target_hpc_x = 750 * np.sin((-ar_target - 100) * 30*u.deg)
        target_hpc_y = 750 * np.cos((-ar_target - 100) * 30*u.deg)
    else:
        raise ValueError("Invalid AR number")

    return float(target_hpc_x), float(target_hpc_y), ar_target


def get_flare_observations(flare_series, pointings, ar_series, eclipse_and_saa):
    distance_to_flare_start_arcmin = []
    offset_xy = []
    goes_class = []
    flare_time = []
    peak_observable = []
    start_observable = []
    end_observable = []
    percent_observable = []
    end_time = []
    ar_match = []
    peak_time = []
    pointing_hpc_x_all = []
    pointing_hpc_y_all = []
    flare_hpc_x_all = []
    flare_hpc_y_all = []
    this_flare_series = flare_series[fm.PHASE_E[0]: fm.PHASE_E[1]]

    for start_time, this_row in this_flare_series.iterrows():

        #  consider the distance between the flare and the pointing axis
        flare_hpc_x, flare_hpc_y = [this_row['hpc_x'], this_row['hpc_y']]
        flare_time.append(start_time)
        goes_class.append(this_row['goes_class'])
        end_time.append(this_row['end_time'])
        # find out where was pointing for flare start
        # this_date = start_time.date().strftime(common_datetime_fmt)
        this_date = start_time.strftime(common_datetime_fmt)
        pointing_hpc_x, pointing_hpc_y, target_noaa = get_pointing(this_date, pointings, ar_series)
        if pointing_hpc_y is not None and pointing_hpc_y is not None:
            offset_xy.append(u.Quantity([flare_hpc_x - pointing_hpc_x,
                                         flare_hpc_y - pointing_hpc_y], 'arcsec'))
            r = np.sqrt((offset_xy[-1]**2).sum())
            distance_to_flare_start_arcmin.append(r.to_value('arcmin'))
        else:
            distance_to_flare_start_arcmin.append(np.nan)

        if this_row['noaa'] == target_noaa:
            ar_match.append(True)
        else:
            ar_match.append(False)

        #  check if the start, peak, and end of the flare are observable
        peak_observable.append(is_sun_observed(str(this_row['peak_time']), eclipse=eclipse_and_saa))
        start_observable.append(is_sun_observed(str(start_time), eclipse=eclipse_and_saa))
        # end_observable.append(is_sun_observed(str(this_row['end_time']), eclipse=eclipse_and_saa))
        end_observable.append(is_sun_observed(str(this_row['peak_time']), eclipse=eclipse_and_saa))
        percent_observable.append(percent_is_sun_observed(start_time, this_row['peak_time'], eclipse_and_saa))
        # percent_observable.append(percent_is_sun_observed(start_time, this_row['end_time'], eclipse_and_saa))
        peak_time.append(this_row['peak_time'])
        pointing_hpc_x_all.append(pointing_hpc_x)
        pointing_hpc_y_all.append(pointing_hpc_y)

        # Get the location of the flare
        flare_hpc_x_all.append(flare_hpc_x)
        flare_hpc_y_all.append(flare_hpc_y)

    result = pd.DataFrame(index=flare_time, data={'end_time': end_time,
                                                  'peak_time': peak_time,
                                                  'goes_class': goes_class,
                                                  'r_hpc_arcmin': distance_to_flare_start_arcmin,
                                                  'offset_xy': offset_xy,
                                                  'start_observable': start_observable,
                                                  'peak_observable': peak_observable,
                                                  'end_observable': end_observable,
                                                  'percent_observable': percent_observable,
                                                  'flare_ar_matches_target_ar': ar_match,
                                                  'pointing_hpc_x': pointing_hpc_x_all,
                                                  'pointing_hpc_y': pointing_hpc_y_all,
                                                  'flare_hpc_x': flare_hpc_x_all,
                                                  'flare_hpc_y': flare_hpc_y_all,
                                                  })
    return result


def _part_of_flare_duration_was_observed(flare_list, percent_observable_lower_limit):
    """helper function to return a whether or not the percentage of the flare that was observed exceeded a lower
    limit."""
    return flare_list['percent_observable'] > percent_observable_lower_limit


def count_flares_in_view(flare_to_pointing, percent_observable_lower_limit=0.0, class_list=('X', 'M')):
    """
    Takes the output of get_flare_observations and provides a count of the flares in the view of
    the observatory.

    Parameters
    ----------
    flare_to_pointing : `~pandas.DataFrame`


    percent_observable_lower_limit : ~float
        The minimum percentage of the flare duration that the flare was observable.

    class_list : ~tuple

    Returns
    -------
    ~tuple
        A two element tuple.  The first element is a dictionary that collects the number of flares
        of each type in the class list that satisfy the observation criteria.  The second element
        is a `~pandas.DataFrame` of the observatory pointing when the observatory observed a
        flare in the full (detector) FOV.
    """
    on_detector = []
    on_detector_and_observed = []
    in_clean_fov = []
    in_clean_fov_and_observed = []

    for i, this_class in enumerate(class_list):
        flare_list = flare_to_pointing[flare_to_pointing['goes_class'].str.contains(this_class)]

        # Was the flare observed for at least some time?
        observed = _part_of_flare_duration_was_observed(flare_list, percent_observable_lower_limit)

        # Was the flare in the detector FOV?
        flare_in_detector_fov = [offset_xy in fm.full_fov for offset_xy in flare_list['offset_xy']]
        on_detector.append(np.sum(flare_in_detector_fov))
        if i == 0:
            on_detector_position = positions_and_flare_class(flare_list, this_class, flare_in_detector_fov)
        else:
            on_detector_position = pd.concat([positions_and_flare_class(flare_list, this_class, flare_in_detector_fov),
                                             on_detector_position])

        # Was the flare in the detector FOV and did the instrument observer for a percentage of the
        # flare duration that was above the lower limit?
        condition = flare_in_detector_fov & observed
        on_detector_and_observed.append(np.sum(condition))
        if i == 0:
            on_detector_and_observed_position = positions_and_flare_class(flare_list, this_class, condition)
        else:
            on_detector_and_observed_position = pd.concat([positions_and_flare_class(flare_list, this_class, condition),
                                                           on_detector_and_observed_position])

        # Was the flare in the clean FOV and did the instrument observer for a percentage of the
        # flare duration that was above the lower limit?

        flare_in_clean_fov = [offset_xy in fm.central_fov for offset_xy in flare_list['offset_xy']]
        in_clean_fov.append(np.sum(flare_in_clean_fov))
        condition = flare_in_clean_fov & observed
        in_clean_fov_and_observed.append(np.sum(condition))
        if i == 0:
            in_clean_fov_and_observed_position = positions_and_flare_class(flare_list, this_class, condition)
        else:
            in_clean_fov_and_observed_position = pd.concat([positions_and_flare_class(flare_list, this_class, condition),
                                                            in_clean_fov_and_observed_position])

    return pd.DataFrame(index=class_list, data={'on detector': on_detector,
                                                'on detector and observed': on_detector_and_observed,
                                                'in clean fov': in_clean_fov,
                                                'in clean fov and observed': in_clean_fov_and_observed}),\
           on_detector_and_observed_position, in_clean_fov_and_observed_position


def positions_and_flare_class(flare_list, flare_class, condition):
    """Return the HPC location given a selection condition.  Also make sure the
    returned dataframe includes the flare class."""
    f = flare_list[['pointing_hpc_x', 'pointing_hpc_y', 'flare_hpc_x', 'flare_hpc_y']][condition]
    f['class'] = flare_class
    return f


def count_flare_ar_matches_target_ar(flare_to_pointing, percent_observable_lower_limit=0.0, class_list=['X', 'M']):
    """Counts the number of times the flare AR matches the target AR as a function of GOES flare class"""
    number_flare_ar_matches_target_ar = []
    number_flare_ar_matches_target_ar_and_observed = []
    for this_class in class_list:
        flare_list = flare_to_pointing[flare_to_pointing['goes_class'].str.contains(this_class)]
        flare_ar_matches_target_ar = flare_list['flare_ar_matches_target_ar']
        number_flare_ar_matches_target_ar.append(np.sum(flare_ar_matches_target_ar))

        observed = _part_of_flare_duration_was_observed(flare_list, percent_observable_lower_limit)
        number_flare_ar_matches_target_ar_and_observed.append(np.sum(flare_ar_matches_target_ar & observed))

    return pd.DataFrame(index=class_list, data={'flare AR matches target AR': number_flare_ar_matches_target_ar,
                                                'flare AR matches target AR and observed': number_flare_ar_matches_target_ar_and_observed})


def plot_orbit_events(start_time, end_time, eclipse, saa, contacts=None, ax=None):
    if ax is None:
        fig, ax = plt.subplots()

    plt.xlim(datetime.datetime.strptime(start_time, common_datetime_fmt),
             datetime.datetime.strptime(end_time, common_datetime_fmt))
    these_eclipses = eclipse[start_time:end_time]
    these_saa = saa[start_time:end_time]

    i = 0
    for index, row in these_eclipses.items():
        plt.axvspan(index.strftime('%x %X'), row.strftime('%x %X'), label='_' * i + 'Eclipse', facecolor='b', alpha=0.5)
        i += 1

    i = 0
    for index, row in these_saa.items():
        plt.axvspan(index.strftime('%x %X'), row.strftime('%x %X'), label='_' * i + 'SAA', facecolor='g', alpha=0.5)
        i += 1

    if contacts is not None:
        these_contacts = contacts[start_time:end_time]
        i = 0
        for index, row in these_contacts.items():
            plt.axvline(index.strftime('%x %X'), label='_' * i + 'Contact', color='black', alpha=1.0)
            i += 1
    # format of the labels
    hfmt = dates.DateFormatter('%x %H:%M')
    ax.xaxis.set_major_formatter(hfmt)


def plot_solar_disk():
    fig, ax = plt.subplots()
    ax.set_aspect('equal')
    plt.xlim(-1100, 1100)
    plt.ylim(-1100, 1100)
    sun = plt.Circle((0, 0), 980, facecolor='lightyellow', edgecolor='orange', zorder=1)
    ax.add_artist(sun)
    return fig, ax


def plot_ars_on_disk(time, ar_series, ax=None):
    if ax is None:
        fig, ax = plot_solar_disk()

    # get color list to keep color of regions the same depending on their NOAA number
    cmap = plt.cm.Paired
    cmaplist = [cmap(i) for i in range(cmap.N)]
    ar_num, noaa, hpc_x, hpc_y, classification = get_ars_for_day(time, ar_series)
    if ar_num != 0:
        colors = [cmaplist[this_noaa % len(cmaplist)] for this_noaa in noaa]
        scatter = ax.scatter(hpc_x, hpc_y, label='AR', color=colors, zorder=2)
    else:
        scatter = None
    return scatter


def plot_flares_on_disk(time, flare_series, ax=None):
    if ax is None:
        fig, ax = plot_solar_disk()
    the_number, goes_class, hpc_x, hpc_y = get_flares(time[0:10], flare_series)
    if the_number != 0:
        colors = ['red' if this_class.count('X') > 0 else 'blue' for this_class in goes_class]
        scatter = ax.scatter(hpc_x, hpc_y, marker='x', zorder=2, color=colors,
                             label='Flare')
    else:
        scatter = None
    return scatter


def detector_fov_coverage(pointing):
    """
    Returns the lower-left corner of a rectangular detector FOV, and its width and height.   Note that coordinates are
    assumed to be helioprojective Cartesian and distances/locations are measured in arcseconds.

    Parameters
    ----------
    pointing : ~tuple
        A tuple describing the center of the detector FOV in arcseconds.

    Returns
    -------
    xy, width, height : ~tuple
        A tuple describing the lower left hand corner of the detector FOV, its width and its height in arcseconds.

    """
    xy = np.array(pointing[0:2]) - fm.FULL_FOV_SIDE.to('arcsec').value / 2.
    width = fm.FULL_FOV_SIDE.to('arcsec').value
    height = fm.FULL_FOV_SIDE.to('arcsec').value
    return xy, width, height


def fov_coverage_mask(xy, fov, shape):
    """
    Create a FOV coverage mask, given the position of FOV, and the shape of the array the FOV is
    painted in to.

    Parameters
    ----------
    xy : ~tuple
        A tuple describing the position of the FOV in arcseconds.
    fov : `~.fov.FOV`
        A FOV object
    shape : ~tuple
        A length-2 tuple that describes the shape of the coverage mask.
    """

    # Pointing relative to Sun center in pixels, assuming that the center of the Sun is at the
    # center of the shape
    rx = xy[1] + shape[1]/2
    ry = xy[0] + shape[0]/2
    if fov.shape in ['square', 'rectangle']:
        mask = np.zeros(shape)
        width = fov.width.to(u.arcsec).value
        height = fov.height.to(u.arcsec).value
        bottom_left_x = rx - height/2
        bottom_left_y = ry - width/2
        mask[int(bottom_left_x):int(bottom_left_x + width), int(bottom_left_y):int(bottom_left_y + height)] = 1
    elif fov.shape == 'circle':
        # Adapted from https://stackoverflow.com/questions/10031580/how-to-write-simple-geometric-shapes-into-numpy-arrays
        xx, yy = np.mgrid[:shape[0], :shape[1]]
        circle = (xx - rx) ** 2 + (yy - ry) ** 2
        mask = 1.0*(circle < (fov.diameter.to(u.arcsec).value/2) ** 2)
    else:
        raise ValueError('FOV shape is not recognized.')
    return mask


def get_mission_pointings(pointing_commands, ar_series, cadence=pandas.Timedelta(24, 'hour')):
    """
    Get where the mission was pointed between the first pointing and the last, sampling at a
    user-defined cadence

    Parameters
    ----------
    pointing_commands : `~pandas.Dataframe`
        A data frame of pointing commands.

    ar_series : `~pandas.Dataframe`
        A data frame

    cadence : `~pandas.Timedelta`
        The pointing sampling cadence.

    Return
    ------
    pandas.DataFrame
        A dataframe that contains the pointing and target at each sample time.
    """
    start_time = pointing_commands.index.min()
    end_time = pointing_commands.index.max()
    t = start_time
    hpc_x = []
    hpc_y = []
    times = []
    targets = []
    while t < end_time:
        times.append(t)
        pointing = get_pointing(t, pointing_commands, ar_series)
        hpc_x.append(pointing[0])
        hpc_y.append(pointing[1])
        targets.append(pointing[2])
        t += cadence
    return pd.DataFrame({'hpc_x': hpc_x, 'hpc_y': hpc_y, 'targets': targets},
                        columns=['hpc_x', 'hpc_y', 'targets'],
                        index=times)


def make_heatmap(pointings, fov, pointing_heatmap_shape):
    """
    Create a heatmap given a list of pointings, the FOV at each pointing, and the size of the
    heatmap the FOVs will be pasted into.

    :param pointings:
    :param fov:
    :param pointing_heatmap_shape:
    :return:
    """
    heatmap = np.zeros(pointing_heatmap_shape)
    ps = list(zip(pointings['hpc_x'].to_list(), pointings['hpc_y'].to_list()))
    for p in ps:
        mask = fov_coverage_mask(p, fov, pointing_heatmap_shape)
        heatmap += mask
    return heatmap


def plot_pointing(time, pointing_series, ar_series, ax=None):
    if ax is None:
        fig, ax = plot_solar_disk()
    pointing = get_pointing(time, pointing_series, ar_series)
    if pointing is not None:
        clean_fov = plt.Circle((pointing[0], pointing[1]), fm.CENTRAL_FOV_DIAMETER.to('arcsec').value/2., edgecolor='red',
                               zorder=1, fill=False)
        ax.add_artist(clean_fov)

        # Get the coverage of this pointing
        xy, width, height = detector_fov_coverage(pointing)

        # Create a rectangle
        detector_fov = plt.Rectangle(xy, width=width, height=height, edgecolor='black', fill=False)
        ax.add_artist(detector_fov)
        # plt.plot(pointing[0], pointing[1], 'x', label='Pointing Vector', markersize=12)


def plot_flare_timeseries(start_time, end_time, goes_series, eclipse, saa):
    fig, ax = plt.subplots()

    plot_orbit_events(start_time, end_time, eclipse, saa, ax=ax)
    goes_series.plot(ax=ax, color='red', label='GOES 1-8 $\AA$')

    ax.set_ylim(1e-9, 1e-2)
    ax.set_yscale("log")
    ax.set_ylabel('Watts m$^{-2}$')
    ax.yaxis.grid(True, 'major')
    ax.xaxis.grid(False, 'major')

    ax2 = ax.twinx()
    ax2.set_yscale("log")
    ax2.set_ylim(1e-9, 1e-2)
    ax2.set_yticks((1e-9, 1e-8, 1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2))
    ax2.set_yticklabels((' ', 'A', 'B', 'C', 'M', 'X', ' '))
    return fig, ax


def plot_all_flare_timeseries(flare_files, flare_series, eclipse_series,
                              saa_series, plot_index=None):
    from sunpy.time import parse_time
    from sunpy.timeseries import TimeSeries
    for start_time, this_row in flare_files.iterrows():
        goes_fits = this_row['filename']
        this_flare = flare_series[flare_series.index == start_time]

        plot_start_time = (parse_time(start_time) - datetime.timedelta(minutes=60)).strftime(common_datetime_fmt)
        plot_end_time = (parse_time(this_flare['end_time'][0]) + datetime.timedelta(minutes=60)).strftime(common_datetime_fmt)

        goes = TimeSeries(goes_fits)
        goes_shifted_series = shift_by_solar_cycle(goes.data)
        fig, ax = plot_flare_timeseries(plot_start_time, plot_end_time,
                                        goes_shifted_series['xrsb'],
                                        eclipse_series, saa_series)
        ax.axvspan(parse_time(start_time), parse_time(this_flare['end_time'][0]), alpha=0.2, color='red',
                   label='NOAA Flare')
        ax.axvline(parse_time(this_flare['peak_time'][0]), label="NOAA Peak time", color='grey')
        ax.set_title(f'GOES {this_flare["goes_class"]} {start_time} ({start_time - SOLAR_CYCLE_DT})')
        ax.legend()

        fig.autofmt_xdate()
        time_string = str(start_time).replace('-', '').replace(' ', '_').replace(':', '')
        pdf_filename = f'flare_timeseries_{this_row["goes_class"][0:2]}_{time_string}.pdf'
        log.info(f'Saving to {pdf_filename}')
        plt.savefig(pdf_filename)

