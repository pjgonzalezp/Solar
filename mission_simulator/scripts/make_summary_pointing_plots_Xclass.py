import pandas as pd
import matplotlib.pyplot as plt

#import flare_mission_sim as fm
import __init__ as fm
#from flare_mission_sim import flare_mission_sim as ms
import flare_mission_sim as ms

pointing_commands = pd.read_csv(fm.in_output_dir('pointing_changes', 'csv'),
                                parse_dates=True, index_col=0)
flare_to_pointing = pd.read_csv(fm.in_output_dir('flare_to_pointing', 'csv'),
                                parse_dates=True, index_col=0)

# load flare and active region data
ar_data = ms.load_ar_data()
flare_data = ms.load_flare_data_from_hdf5()
flare_data = ms.shift_by_solar_cycle(flare_data)
ar_data = ms.shift_by_solar_cycle(ar_data)


# create summary plots for all X class flares
x_class_flare_list = flare_to_pointing[flare_to_pointing['goes_class'].str.contains('X')]

for start_time, this_row in x_class_flare_list.iterrows():
    this_date = start_time.date().strftime(ms.common_datetime_fmt)
    fig, ax = ms.plot_solar_disk()
    ms.plot_ars_on_disk(this_date, ar_data, ax=ax)
    ms.plot_flares_on_disk(this_date, flare_data, ax=ax)
    ms.plot_pointing(this_date, pointing_commands, ar_data, ax=ax)
    ax.set_title('{date} {goes_class} r={distance:0.2f} arcmin'.format(date=this_date,
                                                                       goes_class=this_row['goes_class'],
                                                                       distance=this_row['r_hpc_arcmin']))
    time_string = str(start_time).replace('-', '').replace(' ', '_').replace(':', '')
    pdf_filename = 'flare_pointing_{0}_{1}'.format(this_row['goes_class'][0:2], time_string)

    plt.savefig(fm.in_output_dir(pdf_filename, 'pdf'))
