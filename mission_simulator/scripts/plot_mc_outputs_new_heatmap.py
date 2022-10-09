import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Rectangle
import seaborn as sns
import numpy as np
from matplotlib import cm
from matplotlib.ticker import NullFormatter
import astropy.units as u
#import flare_mission_sim as fm
import __init__ as fm
import flare_mission_sim as ms

sns.set_context('paper')
sns.set_style('ticks')
plt.rcParams['figure.figsize'] = (10, 6)
plt.rcParams['text.usetex'] = True


figsize = (12, 4)
fontsize = 14
labelsize = 12
subplots_adjust = {"left": 0.1, "bottom": 0.2, "right": 0.95, "top": 0.95}
mean_line = {"color": 'black', "linestyle": 'dashed'}
percentile_line = {"color": 'red', "linestyle": 'dotted'}
grid_style = {"linestyle": "dotted"}
output_format = 'pdf'

fov_side = fm.FULL_FOV_SIDE.to(u.arcsec).value
diameter = fm.CENTRAL_FOV_DIAMETER.to(u.arcsec).value

heatmap_map = '$H(x, y)$'
heatmap_lat = '$H_{\mbox{lat}}(y)$'
heatmap_lon = '$H_{\mbox{lon}}(x)$'
heatmap_north = '$H_{\mbox{north}}$'
heatmap_south = '$H_{\mbox{south}}$'
heatmap_west = '$H_{\mbox{west}}$'
heatmap_east = '$H_{\mbox{east}}$'

meridian_linestyle = 'dotted'
equator_linestyle = 'dashed'
pole_linestyle = (0, (3, 10, 1, 10))
ew_limb_linestyle = (0, (5, 10))

nullfmt = NullFormatter()  # no labels


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


def plot_flare_observing_total_histograms(xflare_observing_totals, mflare_observing_totals):
    """
    What does "observe" mean here?  It means that the instrument observed the flare while it was in progress.  If the
    instrument saw ANY part of the flare, this counts as a successful observation.
    """
    plt.figure(1, figsize=figsize)
    plt.subplots_adjust(**subplots_adjust)
    ylim = [0, 1]

    # X-class flares
    bins = regular_hist_bins(xflare_observing_totals, 1)
    weights = hist_to_probability_weight(xflare_observing_totals)
    plt.subplot(1, 2, 1)
    plt.hist(xflare_observing_totals, bins=bins, weights=weights, density=False, align='left')
    plt.ylabel('probability', fontsize=fontsize)
    plt.xlabel('Number of X-class flares observed', fontsize=fontsize)
    plt.tick_params(labelsize=labelsize)
    mn = np.mean(xflare_observing_totals)
    percentiles = np.percentile(xflare_observing_totals, [2.5, 97.5])
    plt.axvline(mn, label=f'mean={mn:.1f}', **mean_line)
    plt.axvline(percentiles[0], **percentile_line)
    plt.axvline(percentiles[1], label=f'95% interval is {percentiles[0]:.1f} to {percentiles[1]:.1f}', **percentile_line)
    plt.legend()
    plt.grid(**grid_style)
    plt.ylim(ylim)

    # M-class flares
    bins = regular_hist_bins(mflare_observing_totals, 5)
    weights = hist_to_probability_weight(mflare_observing_totals)
    plt.subplot(1, 2, 2)
    plt.hist(mflare_observing_totals, bins=bins, weights=weights, density=False)
    plt.ylabel('probability', fontsize=fontsize)
    plt.xlabel('Number of M-class flares observed', fontsize=fontsize)
    plt.tick_params(labelsize=labelsize)
    mn = np.mean(mflare_observing_totals)
    plt.axvline(mn, label=f'mean={mn:.1f}', **mean_line)
    percentiles = np.percentile(mflare_observing_totals, [2.5, 97.5])
    plt.axvline(percentiles[0], **percentile_line)
    plt.axvline(percentiles[1], label=f'95% interval is {percentiles[0]:.1f} to {percentiles[1]:.1f}', **percentile_line)
    plt.legend()
    plt.grid(**grid_style)
    plt.ylim(ylim)

    # Save the complete output
    plt.savefig(fm.in_output_dir('flare_observing_total_histograms', output_format))
    plt.close()


def plot_observing_efficiency_histogram(xflare_observing_fractions, mflare_observing_fractions):
    """
    What does "observe" mean here?  It means that the instrument observed the flare while it was in progress.  If the
    instrument saw ANY part of the flare, this counts as a successful observation.
    """
    plt.figure(1, figsize=figsize)
    plt.subplots_adjust(**subplots_adjust)
    ylim = [0, 1]

    # X-class flares
    bins = regular_hist_bins(xflare_observing_fractions, 0.025)
    weights = hist_to_probability_weight(xflare_observing_fractions)
    plt.subplot(1, 2, 1)
    plt.hist(xflare_observing_fractions, bins=bins, weights=weights)
    plt.ylabel('probability', fontsize=fontsize)
    plt.xlabel('fraction of all X-class flares observed', fontsize=fontsize)
    plt.tick_params(labelsize=labelsize)
    mn = np.mean(xflare_observing_fractions)
    plt.axvline(mn, label=f'mean={mn:.3f}', **mean_line)
    percentiles = np.percentile(xflare_observing_fractions, [2.5, 97.5])
    plt.axvline(percentiles[0], **percentile_line)
    plt.axvline(percentiles[1], label=f'95% interval is {percentiles[0]:.3f} to {percentiles[1]:.3f}', **percentile_line)
    plt.legend()
    plt.grid(**grid_style)
    plt.xlim([0.0, 1.0])
    plt.ylim(ylim)

    bins = regular_hist_bins(mflare_observing_fractions, 0.025)
    weights = np.ones_like(mflare_observing_fractions)/len(mflare_observing_fractions)
    plt.subplot(1, 2, 2)
    plt.hist(mflare_observing_fractions, bins=bins, weights=weights)
    plt.ylabel('probability', fontsize=fontsize)
    plt.xlabel('fraction of all M-class flares observed', fontsize=fontsize)
    plt.tick_params(labelsize=labelsize)
    mn = np.mean(mflare_observing_fractions)
    plt.axvline(mn, label=f'mean={mn:.3f}', **mean_line)
    percentiles = np.percentile(mflare_observing_fractions, [2.5, 97.5])
    plt.axvline(percentiles[0], **percentile_line)
    plt.axvline(percentiles[1], label=f'95% interval is {percentiles[0]:.3f} to {percentiles[1]:.3f}', **percentile_line)
    plt.legend()
    plt.grid(**grid_style)
    plt.xlim([0.0, 1.0])
    plt.ylim(ylim)

    # Save the final plot
    plt.savefig(fm.in_output_dir('flare_observing_fraction_histograms', output_format))
    plt.close()


def plot_flare_percentages(flare_percentages):
    """
    What does "observe" mean here?  It means that the instrument observed the flare while it was in progress.  If the
    instrument saw ANY part of the flare, this counts as a successful observation.
    """
    plt.figure(1, figsize=figsize)
    plt.subplots_adjust(**subplots_adjust)

    bins = regular_hist_bins(flare_percentages, 0.05)
    weights = hist_to_probability_weight(flare_percentages)
    plt.subplot(1, 2, 1)
    plt.hist(flare_percentages, bins=bins, weights=weights)
    plt.ylabel('p', fontsize=fontsize)
    plt.xlabel('Fraction of flare observed', fontsize=fontsize)
    plt.tick_params(labelsize=labelsize)
    plt.grid(**grid_style)

    plt.subplot(1, 2, 2)
    plt.hist(flare_percentages * 100.0, bins=xbins*100.0, cumulative=True, normed=True)
    plt.ylabel('Cumulative fraction', fontsize=fontsize)
    plt.xlabel('Percentage of flare observed', fontsize=fontsize)
    plt.tick_params(labelsize=labelsize)
    plt.grid(**grid_style)

    # Save the final plot
    plt.savefig(fm.in_output_dir('flare_percentage_observed_histogram', output_format))
    plt.close()


def plot_pointing_efficiency_histogram(x, m):
    plt.figure(1, figsize=figsize)
    plt.subplots_adjust(**subplots_adjust)

    # X-class flares
    bins = regular_hist_bins(x, 1)

    plt.subplot(1, 2, 1)
    plt.hist(x, bins=bins, density=True, align='left')
    plt.ylabel('probability density', fontsize=fontsize)
    plt.xlabel('pointing efficiency (X-class)', fontsize=fontsize)
    plt.tick_params(labelsize=labelsize)
    mn = np.mean(x)
    plt.axvline(mn, label=f'mean={mn:.1f}', **mean_line)
    plt.legend()
    plt.grid(**grid_style)

    # M-class flares
    bins = regular_hist_bins(m, 1)

    plt.subplot(1, 2, 2)
    plt.hist(m, bins=bins, density=True, align='left')
    plt.ylabel('probability density', fontsize=fontsize)
    plt.xlabel('pointing efficiency (M-class)', fontsize=fontsize)
    plt.tick_params(labelsize=labelsize)
    mn = np.mean(m)
    plt.axvline(mn, label=f'mean={mn:.1f}', **mean_line)
    plt.legend()
    plt.grid(**grid_style)

    # Save the final plot
    plt.savefig(fm.in_output_dir('flare_observing_pointing_efficiency', output_format))
    plt.close()


def _ratio_west_to_east(heatmap):
    """Calculate the heatmap ratio of the western half to the eastern half"""
    image_half_length = _image_half_length(heatmap)
    eastern = heatmap[:, 0:image_half_length]
    western = heatmap[:, image_half_length:]
    return western.sum()/eastern.sum()


def _ratio_north_to_south(heatmap):
    """Calculate the heatmap ratio of the western half to the eastern half"""
    image_half_length = _image_half_length(heatmap)
    south = heatmap[0:image_half_length, :]
    north = heatmap[image_half_length+1:, :]
    return north.sum()/south.sum()


def _image_half_length(heatmap):
    return int(heatmap.shape[0]/2)


def _y_position_calculator(y):
    y_min = y.min()
    y_max = y.max()
    y0 = y_min + 0.25*(y_max - y_min)
    y1 = y_min + 0.30*(y_max - y_min)
    return y0, y1


def _longitude_heat(heatmap, filepath='longitude_heat', title=None):
    # Calculate the mean heat
    longitude_heat = np.mean(heatmap, axis=0)

    # Get the longitude range
    image_half_length = _image_half_length(heatmap)
    longitude = np.arange(-image_half_length, image_half_length)

    # Calculate the ratio of the total heat: western hemisphere / eastern hemos[here
    ratio = np.sum(longitude_heat[image_half_length+1:]) / np.sum(longitude_heat[0:image_half_length])

    # Y positions of the scale lengths
    y0, y1 = _y_position_calculator(longitude_heat)

    # Make the plot
    fig, ax = plt.subplots(1)
    ax.plot(longitude, longitude_heat, label='mean heat over all latitudes')
    ax.set_xlabel('solar X')
    ax.set_ylabel('heat')
    ax.set_title(f"heat as a function of solar X\n{title}")
    ax.axvline(0, label='disk center', linestyle='--', color='r')
    ax.axvline(-960, label='eastern and western limbs', linestyle=':', color='r')
    ax.axvline(960, linestyle=':', color='r')
    ax.plot([-fov_side/2, fov_side/2],
            [y0, y0], linewidth=2, label='full lengthscale', color='k')
    ax.plot([-diameter/2, diameter/2],
            [y1, y1], linewidth=2, label='clean diameter', color='k', linestyle=":")
    ax.text(100, 0, 'western\nhemisphere', color='k')
    ax.text(-100, 0, 'eastern\nhemisphere', ha='right', color='k')
    ax.text(0, 0.8*np.max(longitude_heat), 'ratio west/east = {:.1f}'.format(ratio), color='k')
    ax.set_xlim(-image_half_length, image_half_length)
    ax.grid('on', linestyle=':')
    ax.legend()

    # Save the final plot
    plt.savefig(fm.in_output_dir(filepath, output_format))
    plt.close()


def _latitude_heat(heatmap, filepath='longitude_heat', title=None):
    # Calculate the mean heat
    latitude_heat = np.mean(heatmap, axis=1)

    # Get the latitude range
    image_half_length = _image_half_length(heatmap)
    latitude = np.arange(-image_half_length, image_half_length)

    # Calculate the ratio of the total heat: northern hemisphere / southern hemisphere
    ratio = np.sum(latitude_heat[image_half_length+1:]) / np.sum(latitude_heat[0:image_half_length])

    # Y positions of the scale lengths
    y0, y1 = _y_position_calculator(latitude_heat)

    # Make the plot
    fig, ax = plt.subplots(1)
    ax.plot(latitude, latitude_heat, label='mean heat over all longitudes')
    ax.set_xlabel('solar Y')
    ax.set_ylabel('heat')
    ax.set_title(f"heat as a function of solar Y\n{title}")
    ax.axvline(0, label='equator', linestyle='--', color='r')
    ax.axvline(-960, label='northern and southern poles', linestyle=':', color='r')
    ax.axvline(960, linestyle=':', color='r')
    ax.plot([-fov_side/2, fov_side/2],
            [y0, y0], linewidth=2, label='full lengthscale', color='k')
    ax.plot([-diameter/2, diameter/2],
            [y1, y1], linewidth=2, label='clean diameter', color='k', linestyle=":")
    ax.text(100, 0, 'northern\nhemisphere', color='k')
    ax.text(-100, 0, 'southern\nhemisphere', ha='right', color='k')
    ax.text(0, 0.8*np.max(latitude_heat), 'ratio north/south = {:.1f}'.format(ratio), color='k')
    ax.set_xlim(-image_half_length, image_half_length)
    ax.grid('on', linestyle=':')
    ax.legend()

    # Save the final plot
    plt.savefig(fm.in_output_dir(filepath, output_format))
    plt.close()


def plot_fov_heatmap(heatmap, fov_type, normalization=np.max, clabel=''):
    """
    Plots a spatial histogram of where the FOV was over the course of the mission
    """
    description = 'mission pointing'
    plot_type = 'heatmap'
    fov_usage = f'using {fov_type} FOV'
    sampling = 'sample rate: 1 pointing per day'
    title = f"{description} {plot_type}\n{fov_usage}\n{sampling}"
    l_title = f"{description}\n{fov_usage}\n{sampling}"

    # Properties of the heatmap
    image_half_length = _image_half_length(heatmap)

    # Normalization of the heatmap
    if isinstance(normalization, int):
        n = normalization
    else:
        n = normalization(heatmap)
    heatmap = heatmap/n

    # Eastern and western hemisphere heatmap ratio
    ratio = _ratio_west_to_east(heatmap)

    # Show the image
    fig, ax = plt.subplots(1)
    ax.set_aspect('equal')
    extent = [-image_half_length, image_half_length, -image_half_length, image_half_length]
    im = ax.imshow(heatmap, cmap=cm.inferno, origin='lower', extent=extent)
    limb = Circle((0, 0), radius=960, label='solar disk', fill=False, color='w')
    central = Circle((0, 0), radius=diameter / 2, label='central', fill=False, color='w',
                     linestyle=":")
    full = Rectangle((-fov_side / 2, -fov_side / 2), fov_side, fov_side, fill=False, color='w',
                     linewidth=2, label='FOV')
    ax.add_patch(limb)
    ax.add_patch(full)
    ax.add_patch(central)
    ax.set_xlabel('solar X')
    ax.set_ylabel('solar Y')
    ax.set_title(title)
    ax.axvline(0, color='w', linestyle=':')
    ax.text(100, 500, 'western\nhemisphere', color='w')
    ax.text(-100, -500, 'eastern\nhemisphere', ha='right', color='w')
    ax.text(50, 1000, 'ratio west/east = {:.1f}'.format(ratio), color='w')
    fig.colorbar(im, ax=ax, label=clabel)

    # Save the final plot
    base_filename = f'heatmap_one_fov_per_pointing_fov={fov_type}'
    plt.savefig(fm.in_output_dir(f'{base_filename}__image', output_format))
    plt.close()

    # Make a longitude only plot
    _longitude_heat(heatmap, filepath=f'{base_filename}__longitude', title=l_title)

    # Make a latitude only plot
    _latitude_heat(heatmap, filepath=f'{base_filename}__latitude', title=l_title)


def plot_fov_when_flare_observed(flare_positions, observation_type, shape=(2200, 2200),
                                 class_list=('M', 'X'),
                                 normalization=np.max,
                                 clabel=''):
    """
    Plots the FOVs when a flare was observed
    """
    if observation_type == 'on detector and observed':
        fov = fm.full_fov
        fov_type = 'full'
    else:
        raise ValueError('other observation types not yet supported')

    description = 'mission pointing when flare was observed by mission'
    plot_type = 'heatmap'
    fov_usage = f'using {fov_type} FOV'
    sampling = 'sample rate: 1 FOV per flare observed by mission'
    title = f"{description} ({plot_type})\n{fov_usage}\n{sampling}"
    l_title = f"{description}\n{fov_usage}\n{sampling}"

    # Create the spatial histogram by adding up all the FOV at each pointing
    heatmap = np.zeros(shape)
    image_half_length = _image_half_length(heatmap)
    for flare_class in class_list:
        condition = flare_positions['class'] == flare_class
        x = flare_positions[condition]['pointing_hpc_x'].to_numpy()
        y = flare_positions[condition]['pointing_hpc_y'].to_numpy()

        for i in range(0, len(x)):
            fov_mask = ms.fov_coverage_mask([x[i], y[i]], fov, shape)
            heatmap += fov_mask

    # Normalization
    if isinstance(normalization, int):
        n = normalization
    else:
        n = normalization(heatmap)
    heatmap = heatmap/n

    # Eastern and western hemispheres ratio
    ratio = _ratio_west_to_east(heatmap)

    # Show the image
    fig, ax = plt.subplots(1)
    ax.set_aspect('equal')
    extent = [-image_half_length, image_half_length, -image_half_length, image_half_length]
    im = ax.imshow(heatmap, cmap=cm.inferno, origin='lower', extent=extent)
    ax.set_xlabel('solar X')
    ax.set_ylabel('solar Y')
    ax.set_title(title)
    ax.axvline(0, color='w', linestyle=':')
    ax.text(100, 500, 'western\nhemisphere', color='w')
    ax.text(-100, -500, 'eastern\nhemisphere', ha='right', color='w')
    ax.text(50, 1000, 'ratio west/east = {:.1f}'.format(ratio), color='w')
    limb = Circle((0, 0), radius=960, label='solar disk', fill=False, color='w')
    central = Circle((0, 0), radius=diameter / 2, label='central', fill=False, color='w',
                     linestyle=":")
    full = Rectangle((-fov_side / 2, -fov_side / 2), fov_side, fov_side, fill=False, color='w',
                     linewidth=2, label='FOV')
    ax.add_patch(limb)
    ax.add_patch(full)
    ax.add_patch(central)
    fig.colorbar(im, ax=ax, label=clabel)

    # Save the final plot
    base_filename = f'heatmap_one_fov_per_observed_flare_fov={fov_type}'
    plt.savefig(fm.in_output_dir(f'{base_filename}__image', output_format))
    plt.close()

    # Make a longitude only plot
    _longitude_heat(heatmap, filepath=f'{base_filename}__longitude', title=l_title)

    # Make a latitude only plot
    _latitude_heat(heatmap, filepath=f'{base_filename}__latitude', title=l_title)


def plot_fov_for_all_flares(flare_positions, fov_type, shape=(2200, 2200), normalization=np.max,
                            clabel=''):
    """
    Heatmap built using one FOV per flare regardless of whether the mission saw it or not.
    This is an estimate of the heatmap we would get if we have foresight of all flares.
    """
    # Which FOV
    if fov_type == 'full':
        fov = fm.full_fov
    elif fov_type == 'central':
        fov = fm.central_fov
    else:
        raise ValueError('Unknown FOV specification.')

    description = 'flare occurrence'
    plot_type = 'heatmap'
    fov_usage = f'using {fov_type} FOV'
    sampling = 'sample rate: 1 FOV per flare occurrence'
    title = f"{description} {plot_type}\n{fov_usage}\n{sampling}"
    l_title = f"{description}\n{fov_usage}\n{sampling}"

    # Create the spatial histogram by adding up all the FOV at each pointing
    heatmap = np.zeros(shape)
    image_half_length = _image_half_length(heatmap)
    x = flare_positions['hpc_x'].to_numpy()
    y = flare_positions['hpc_y'].to_numpy()

    # Build the heatmap
    for i in range(0, len(x)):
        fov_mask = ms.fov_coverage_mask([x[i], y[i]], fov, shape)
        heatmap += fov_mask

    # Normalization
    if isinstance(normalization, int):
        n = normalization
    else:
        n = normalization(heatmap)
    heatmap = heatmap/n

    # Eastern and western hemispheres ratio
    ratio = _ratio_west_to_east(heatmap)

    # Show the image
    fig, ax = plt.subplots(1)
    ax.set_aspect('equal')
    extent = [-image_half_length, image_half_length, -image_half_length, image_half_length]
    im = ax.imshow(heatmap, cmap=cm.inferno, origin='lower', extent=extent)
    ax.set_xlabel('solar X')
    ax.set_ylabel('solar Y')
    ax.set_title(title)
    ax.axvline(0, color='w', linestyle=':')
    ax.text(100, 500, 'western\nhemisphere', color='w')
    ax.text(-100, -500, 'eastern\nhemisphere', ha='right', color='w')
    ax.text(50, 1000, 'ratio west/east = {:.1f}'.format(ratio), color='w')
    limb = Circle((0, 0), radius=960, label='solar disk', fill=False, color='w')
    central = Circle((0, 0), radius=diameter / 2, label='central', fill=False, color='w',
                     linestyle=":")
    full = Rectangle((-fov_side / 2, -fov_side / 2), fov_side, fov_side, fill=False, color='w',
                     linewidth=2, label='FOV')
    ax.add_patch(limb)
    ax.add_patch(full)
    ax.add_patch(central)
    fig.colorbar(im, ax=ax, label=clabel)

    # Save the final plot
    base_filename = f'heatmap_one_fov_per_flare_fov={fov_type}'
    plt.savefig(fm.in_output_dir(f'{base_filename}__image', output_format))
    plt.close()

    # Make a longitude only plot
    _longitude_heat(heatmap, filepath=f'{base_filename}__longitude', title=l_title)

    # Make a latitude only plot
    _latitude_heat(heatmap, filepath=f'{base_filename}__latitude', title=l_title)


def plot_heatmaps_and_summaries(heatmap, filename, title, figsize=(8, 8), left=0.15, width=0.6,
                                bottom=0.1, height=0.6, offset=0.03, l_size=0.17):
    """
    Plot a two-dimensional heatmap and its summaries in a single plot.

    figsize

    left

    width

    bottom

    height

    offset
    """

    ratio_we = '{:.1f}'.format(_ratio_west_to_east(heatmap))
    ratio_ns = '{:.1f}'.format(_ratio_north_to_south(heatmap))

    # definitions for the axes
    bottom_h = left + width + offset
    left_h = left + width + offset

    rect_heatmap = [left, bottom, width, height]
    rect_longitude = [left, bottom_h, width, l_size]
    rect_latitude = [left_h, bottom, l_size, height]
  #  rect_colorbar = [left_h, bottom_h, 0.1, l_size]
    rect_colorbar = [left_h, bottom, 0.05, height]

    # start with a rectangular Figure
    fig = plt.figure(9, figsize=figsize)

 #   plt.suptitle(title)

    ax_heatmap = plt.axes(rect_heatmap)
    ax_longitude = plt.axes(rect_longitude)
 #   ax_latitude = plt.axes(rect_latitude)
    ax_colorbar = plt.axes(rect_colorbar)

    # no labels
    ax_longitude.xaxis.set_major_formatter(nullfmt)
  #  ax_latitude.yaxis.set_major_formatter(nullfmt)

    image_half_length = _image_half_length(heatmap)

    img = heatmap/heatmap.max()
    img_masked =np.ma.masked_where((img == 0), img)
    cmap = cm.Oranges
    cmap.set_bad('tab:purple',1.0)
    
    # Two-dimensional heatmap plot
    extent = [-image_half_length, image_half_length, -image_half_length, image_half_length]
    im = ax_heatmap.imshow(img_masked, cmap=cmap, origin='lower', extent=extent, interpolation = None)
    limb = Circle((0, 0), radius=960, label='solar disk', fill=False, color='black')
    central = Circle((0, 0), radius=diameter / 2, label='central', fill=False, color='black',
                     linestyle=":")
    full = Rectangle((-fov_side / 2, -fov_side / 2), fov_side, fov_side, fill=False, color='black',
                     linewidth=2, label='FOV')
    ax_heatmap.add_patch(limb)
    ax_heatmap.add_patch(central)
    ax_heatmap.add_patch(full)
    ax_heatmap.set_xlabel('solar X (arcsec)', fontsize=20)
    ax_heatmap.set_ylabel('solar Y (arcsec)', fontsize=20)
    ax_heatmap.axhline(0, color='black', linestyle=equator_linestyle)
    ax_heatmap.axvline(0, color='black', linestyle=meridian_linestyle)
    ax_heatmap.grid('on', linestyle=':', color='grey', linewidth=0.1)
   # ax_heatmap.set_title('normalized heatmap')
    ax_heatmap.tick_params(labelsize=20)
    cb = fig.colorbar(im, cax=ax_colorbar)
    cb.ax.tick_params(labelsize=20)
    cb.ax.set_ylabel('{:s}/max({:s})'.format(heatmap_map, heatmap_map), fontsize=18)

    # Heat as a function of longitude
    equal_extent = np.arange(-image_half_length, image_half_length)
    longitude_heat = np.mean(heatmap, axis=0)
    ratio = '{:s}/{:s}={:s}'.format(heatmap_west, heatmap_east, ratio_we)
    ylabel = '{:s}/max({:s})'.format(heatmap_lon, heatmap_lon)
    ax_longitude.plot(equal_extent, longitude_heat/longitude_heat.max())
    ax_longitude.tick_params(labelsize=20)
    ax_longitude.grid('on', linestyle=':')
    ax_longitude.set_ylabel(ylabel, fontsize=18)
    ax_longitude.axvline(-960, label='eastern/western limb', linestyle=ew_limb_linestyle, color='r')
    ax_longitude.axvline(960, linestyle=ew_limb_linestyle, color='r')
    ax_longitude.axvline(0, label='meridian', linestyle=meridian_linestyle, color='r')
    ax_longitude.set_xlim(-image_half_length, image_half_length)
    ax_longitude.text(50, 0.1, ratio, fontsize=20)
    ax_longitude.set_title('normalized heat vs. x (summed over all y)', fontsize=18)

    # Heat as a function of longitude
    latitude_heat = np.mean(heatmap, axis=1)
    ones = np.ones_like(latitude_heat)
    l_h_max = latitude_heat/latitude_heat.max()
    ratio = '{:s}/{:s}={:s}'.format(heatmap_north, heatmap_south, ratio_ns)
  #  ax_latitude.plot(l_h_max, equal_extent)
  #  ax_latitude.grid('on', linestyle=':')
  #  xlabel = '{:s}/max({:s})'.format(heatmap_lat, heatmap_lat)
  #  ax_latitude.set_xlabel(xlabel)
  #  ax_latitude.plot(l_h_max, -960*ones, label='northern/southern poles', linestyle=pole_linestyle, color='r')
  #  ax_latitude.plot(l_h_max, 960*ones, linestyle=pole_linestyle, color='r')
  #  ax_latitude.plot(l_h_max, 0*ones, label='equator', linestyle=equator_linestyle, color='r')
   # ax_latitude.set_ylim(-image_half_length, image_half_length)
    ax_heatmap.text(-1000, 900, ratio, fontsize=20, color = 'white')
    ax_heatmap.text(450, 950, 'Scenario H4', fontsize=20, color='white')
     
   # ax_latitude.text(0.1, 50, ratio)
   # ax_latitude.text(1.1, 0.5, 'normalized heat as a function of latitude (summed over all longitudes)',
    #                 transform=ax_latitude.transAxes, horizontalalignment='center',
    #                 verticalalignment='center', rotation=-90)
    #ax_latitude.yaxis.set_label_position("right")

    # Save the final plot
    plt.tight_layout()
    plt.savefig(fm.in_output_dir(filename, output_format))
    plt.close()


def plot_map_of_observed_flares(phase_e_flares, flares_in_view, heatmap_shape, filename,
                                figsize=(8, 8), left=0.1, width=0.6,
                                bottom=0.1, height=0.6, offset=0.03, l_size=0.2):

    def _ratio_ns_ew(x, y, prefix=''):
        r_ns = '$r_{\mbox{N/S}}$'
        r_ew = '$r_{\mbox{E/W}}$'
        ew = '{:.1f}'.format(np.sum(x < 0) / np.sum(x > 0))
        ns = '{:.1f}'.format(np.sum(y > 0) / np.sum(y < 0))
        return '{:s}{:s}={:s}'.format(prefix, r_ew, ew), '{:s}{:s}={:s}'.format(prefix, r_ns, ns)


    fig = plt.figure(9, figsize=figsize)

    # definitions for the axes
    bottom_h = left + width + offset
    left_h = left + width + offset

    rect_heatmap = [left, bottom, width, height]
    rect_longitude = [left, bottom_h, width, l_size]
    rect_latitude = [left_h, bottom, l_size, height]
    #rect_colorbar = [left_h, bottom_h, 0.1, l_size]

    ax_scatter = plt.axes(rect_heatmap)
    ax_longitude = plt.axes(rect_longitude)
    ax_latitude = plt.axes(rect_latitude)
    #ax_colorbar = plt.axes(rect_colorbar)

    ax_longitude.xaxis.set_major_formatter(nullfmt)
    ax_latitude.yaxis.set_major_formatter(nullfmt)

    pef_x = phase_e_flares['hpc_x'].to_numpy()
    pef_y = phase_e_flares['hpc_y'].to_numpy()
    pef_label = 'unobserved M or X'
    pef_color = 'lightsteelblue'
    pef_x_r, pef_y_r = _ratio_ns_ew(pef_x, pef_y, prefix='All: ')

    pwo_x = flares_in_view['hpc_x'].to_numpy()
    pwo_y = flares_in_view['hpc_y'].to_numpy()
    pwo_label = 'pointing when M or X observed'
    pwo_label_v = 'pointing when\nM or X observed'
    pwo_color = 'r'
    pwo_x_r, pwo_y_r = _ratio_ns_ew(pwo_x, pwo_y, prefix='PWO: ')

    fwo_x = flares_in_view['flare_hpc_x'].to_numpy()
    fwo_y = flares_in_view['flare_hpc_y'].to_numpy()
    fwo_label = 'observed M or X'
    fwo_color = 'blue'

    ax_scatter.scatter(pef_x, pef_y, label=pef_label, color=pef_color)
    ax_scatter.scatter(fwo_x, fwo_y, label=fwo_label, color=fwo_color, s=40)
    ax_scatter.scatter(pwo_x, pwo_y, label=pwo_label, color=pwo_color, s=40, marker='x')
    limb = Circle((0, 0), radius=960, fill=False, color='k')
    ax_scatter.add_patch(limb)
    ax_scatter.set_xlabel('solar X (longitude, arcsec)')
    ax_scatter.set_ylabel('solar Y (latitude, arcsec)')
    ax_scatter.axhline(0, color='k', linestyle=equator_linestyle)
    ax_scatter.axvline(0, color='k', linestyle=meridian_linestyle)
    ax_scatter.grid('on', linestyle=':')
    ax_scatter.set_xlim(-heatmap_shape[0]/2, heatmap_shape[0]/2)
    ax_scatter.set_ylim(-heatmap_shape[1]/2, heatmap_shape[1]/2)
    ax_scatter.legend()

    xbins = 30
    ax_longitude.hist(pef_x, bins=xbins, label='all M and X')
    ax_longitude.hist(pwo_x, bins=xbins, range=(pef_x.min(), pef_x.max()), color='r', label=pwo_label)
    ax_longitude.axvline(0, color='k', linestyle=meridian_linestyle)
    ax_longitude.text(-750, 40, pef_x_r + '\n' + pwo_x_r)
    ax_longitude.set_ylabel('number')
    ax_longitude.grid('on', linestyle=":")
    ax_longitude.set_xlim(ax_scatter.get_xlim())
    ax_longitude.legend()

    ax_latitude.hist(pef_y, bins=xbins, orientation='horizontal', label='all M and X')
    ax_latitude.hist(pwo_y, bins=xbins, orientation='horizontal', range=(pef_y.min(), pef_y.max()), color='r', label=pwo_label_v)
    ax_latitude.axhline(0, color='k', linestyle=equator_linestyle)
    ax_latitude.text(5, -750, pef_y_r + '\n' + pwo_y_r)
    ax_latitude.set_xlabel('number')
    ax_latitude.grid('on', linestyle=":")
    ax_latitude.set_ylim(ax_scatter.get_ylim())
    ax_latitude.legend()

  #  plt.tight_layout()
    plt.savefig(fm.in_output_dir(filename, output_format))
    plt.close()


