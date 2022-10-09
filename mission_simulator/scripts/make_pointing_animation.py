# warning, needs imagemagick installed. Probably use homebrew to install.

# create an animation of the pointing throughout the mission.

import matplotlib.animation as animation
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import datetime
import numpy as np

#import flare_mission_sim as fm
import __init__ as fm
#from flare_mission_sim import flare_mission_sim as ms
import flare_mission_sim as ms

# Load in the necessary flare, AR and orbit data
# ----------------------------------------------
# ----------------------------------------------

start_time = fm.PHASE_E[0]
num_days = (fm.PHASE_E[1] - fm.PHASE_E[0]).days

mfw = ms.load_major_flare_watches()
# load flare and active region data
ar_data = ms.load_ar_data()
flare_data = ms.load_flare_data_from_hdf5()
flare_data = flare_data[flare_data['goes_class'].str.contains('X|M')]
flare_data = ms.fix_flare_positions_sff(flare_data)
flare_data = ms.fix_flare_positions_manual(flare_data)

# read in orbit events
eclipse, saa, contacts, polar = ms.read_orbit_events()
#shift orbit data back to solar cycle 24, where we have real flare and AR data
eclipse = ms.shift_by_solar_cycle(eclipse, n=-1)  
saa = ms.shift_by_solar_cycle(saa, n=-1)
contacts = ms.shift_by_solar_cycle(contacts, n=-1)  
eclipse_and_saa = ms.concatenate_range_series(eclipse, saa)

# Run the pointing decisions and pointing commands codes
# ------------------------------------------------------
# ------------------------------------------------------

soc_decisions = ms.generate_pointing_decisions(ar_data[fm.PHASE_E[0]: fm.PHASE_E[1]], mfw=mfw)
pointings = ms.generate_pointing_commands(soc_decisions, contacts,
                                          monte_carlo=False, use_random_delay = True)

# Make the animation
# ------------------
# ------------------


# create the basic plot
fig, ax = ms.plot_solar_disk()
flare_scatter = ax.scatter([], [], marker='x', color='red', zorder=2)
ar_scatter = ax.scatter([], [], zorder=2)

clean_fov = plt.Circle((0, 0), fm.CENTRAL_FOV_DIAMETER.value / 2., edgecolor='red', zorder=2, fill=False)

xy = np.array([0, 0]) - fm.FULL_FOV_SIDE.value / 2.
detector_fov = plt.Rectangle(xy, width=fm.FULL_FOV_SIDE.value,
                             height=fm.FULL_FOV_SIDE.value, edgecolor='black', fill=False, zorder=2)
ax.add_artist(detector_fov)
ax.add_patch(clean_fov)

# create a list of colors to be used by the NOAA active regions so it depends on their NOAA number and does not change.
cmap = plt.cm.tab10
cmaplist = [cmap(i) for i in range(cmap.N)]
# color_dict = {'F':3,'E':1,'D':0,'C':4,'H':2,'B':7,'A':8}


def update_plot(frame_number):
    new_time = start_time + datetime.timedelta(days=int(frame_number))
    new_time_str = new_time.strftime(ms.common_datetime_fmt)

    ar_num, noaa, hpc_x, hpc_y, classification = ms.get_ars(new_time_str[0:10],
                                                            ar_data)
    if ar_num != 0:
        colors = np.array([cmaplist[this_noaa % len(cmaplist)] for this_noaa in noaa])
        # colorlist = []
        # for c in classification:
        #   colorlist.append(color_dict.get(c[0]))
        # colors = np.array([cmaplist[c] for c in colorlist])
        
        ar_scatter.set_offsets(np.vstack([hpc_x, hpc_y]).transpose())
        ar_scatter.set_facecolors(colors)
    else:
        ar_scatter.set_offsets(np.vstack([[], []]).transpose())

    the_number, goes_class, hpc_x, hpc_y = ms.get_flares(new_time_str[0:10], flare_data)
    if the_number != 0:
        flare_scatter.set_offsets(np.vstack([hpc_x, hpc_y]).transpose())
    else:
        flare_scatter.set_offsets(np.vstack([[], []]).transpose())

    hpc_x, hpc_y, target_noaa = ms.get_pointing(new_time_str, pointings,
                                                ar_data)
    pointing = [hpc_x, hpc_y]
    if pointing is not None:
        clean_fov.center = pointing
        xy = np.array(pointing) - fm.FULL_FOV_SIDE.value / 2.
        detector_fov.xy = xy
        # scatter = ax.scatter(hpc_x, hpc_y, marker='x', zorder=2, color='red', label='Flare')
    ax.set_title(new_time, fontsize=16)


anim = FuncAnimation(fig, update_plot, frames=np.arange(0, 100), interval=200, repeat=True)
output_filename = f'pointing_on_disk_vis_{datetime.datetime.now().strftime("%Y%m%d_%H%M")}.mp4'
anim.save(output_filename, dpi=300, writer='imagemagick')
# plt.show()
