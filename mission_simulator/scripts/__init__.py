name = 'mission_sim'

from pathlib import Path
import os.path
import logging
import datetime
import sys

import pandas as pd

import astropy.units as u

from fov import FOV

data_dir = Path(os.path.join(Path(__file__).parents[1], "data"))
out_dir = Path(os.path.join(Path(__file__).parents[1], "out"))

# string format to timestamp output files
now_str_format = "%Y%m%d_%H%M"


def setup_logging(filepath=None):
    """
    Setup up logging to a file as well as to the screen
    """
    log = logging.getLogger()
    log.setLevel("INFO")

    formatter = logging.Formatter('[%(asctime)s - %(filename)s:%(lineno)s]:%(levelname)s - %(message)s')

    if filepath is not None:
        # create file handler which logs even debug messages
        fh = logging.FileHandler(filepath)
        fh.setFormatter(formatter)
        log.addHandler(fh)

    # Add the stdout handler if the logger doesn't already have a stdout handler
    streamhandlers = filter(lambda h: isinstance(h, logging.StreamHandler), log.handlers)
    if sys.stdout not in [sh.stream for sh in streamhandlers]:
        ch = logging.StreamHandler(sys.stdout)
        ch.setFormatter(formatter)
        log.addHandler(ch)

    return log


def _load_mission_parameters(filename):
    """Load the mission configuration file (mission_parameters.csv)."""
    params = pd.read_csv(filename, index_col=0)

    result = {}
    for row in params.iterrows():
        if row[0].count('date') == 1:
            result.update({row[0]: datetime.datetime.strptime(row[1]['value'],
                                                              row[1]['units'])})
        elif any(x in row[0] for x in ('file', 'mission_name')):
            result.update({row[0]: row[1]['value']})
        else:
            result.update({row[0]: u.Quantity(row[1]['value'],
                                              row[1]['units'])})

    return result


# Define the location of the mission parameter file
mission_parameter_file = os.path.join(data_dir, "mission_parameters.csv")

# Load the mission parameters
params = _load_mission_parameters(mission_parameter_file)

# Extract some important parameters from the parameters array
CENTRAL_FOV_DIAMETER = params['central_fov_diameter']
FULL_FOV_SIDE = params['full_fov_diameter']
PHASE_E = [params['date_mission_start'], params['date_mission_end']]
FILE_ORBIT_EVENTS = params['file_orbit_events']
MISSION_NAME = params['mission_name']

# Define FOV instances
central_fov = FOV('circle', diameter=CENTRAL_FOV_DIAMETER)
full_fov = FOV('square', width=FULL_FOV_SIDE)

def in_output_dir(basename, extension, timestamp=None):
    """
    Function that creates a standardized filepath.
    """
    if timestamp is not None:
        output_filename = f'{basename}_{MISSION_NAME}_{timestamp.strftime(now_str_format)}.{extension}'
    else:
        output_filename = f'{basename}_{MISSION_NAME}.{extension}'
    return os.path.join(out_dir, output_filename)
