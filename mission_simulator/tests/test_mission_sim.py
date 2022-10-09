
from flare_mission_sim import flare_mission_sim
import flare_mission_sim as fm
from flare_mission_sim import util
from numpy.testing import assert_almost_equal
import numpy as np

def test_load_ar_data():

    result = flare_mission_sim.load_ar_data()
    assert len(result.columns) == 7
    assert result.columns[0] == 'noaa'
    assert result.columns[1] == 'number'
    assert result.columns[2] == 'hpc_x'
    assert result.columns[3] == 'hpc_y'
    assert result.columns[4] == 'classification'
    assert result.columns[5] == 'hale_classification'
    assert result.columns[6] == 'numspots'

    # check a known entry is correct
    entry = result.loc['2011-02-15']
    assert entry['noaa'] == '11157,11158,11161,11159' 
    assert entry['classification'] == 'DSO,EKC,BXO,BXO'
    assert entry['hpc_x'] == '628.524,220.4454,-753.708,223.191'

    
def test_load_and_fix_flare_data():
    flare_data = flare_mission_sim.load_flare_data_from_hdf5()
    entries = flare_data.loc['2011-02-15']
    assert len(entries) == 8

    # remove all C class flares
    flare_data = flare_data[flare_data['goes_class'].str.contains('X|M')]

    #test two known entries
    
    #entry1 has a bad position that is fixed by cross-reference with SFF database
    entry1 = flare_data.loc['2011-02-15']
    assert len(entry1) == 1
    assert entry1['hpc_x'].values == 0.0
    assert entry1['hpc_y'].values == 115.764

    # entry2 has a bad position that is fixed later by manual correction
    entry2 = flare_data.loc['2011-09-24 16:36:00']
    assert entry2['hpc_x'] == 0.0
    assert entry2['hpc_y'] == -116.8656
    
    # fix bad flare positions by cross-referencing with SFF databases
    # re-check entry1
    flare_data = flare_mission_sim.fix_flare_positions_sff(flare_data)
    entry1 = flare_data.loc['2011-02-15']
    assert len(entry1) == 1
    assert entry1['hpc_x'].values == 159.166
    assert entry1['hpc_y'].values == -224.146

    # re-check entry2
    entry2 = flare_data.loc['2011-09-24 16:36:00']
    assert entry2['hpc_x'] == 0.0
    assert entry2['hpc_y'] == -116.8656

    # fix remaining bad flare positions using manual correction
    # re-check entry1
    flare_data = flare_mission_sim.fix_flare_positions_manual(flare_data)
    entry1 = flare_data.loc['2011-02-15']
    assert len(entry1) == 1
    assert entry1['hpc_x'].values == 159.166
    assert entry1['hpc_y'].values == -224.146

    # re-check entry2
    entry2 = flare_data.loc['2011-09-24 16:36:00']
    assert entry2['hpc_x'] == 883.0  
    assert entry2['hpc_y'] == 345.0


def test_get_ar_positions_over_time():

    ar_data = flare_mission_sim.load_ar_data()
    entries = flare_mission_sim.get_ar_positions_over_time(11158,ar_data)
    assert len(entries) == 9
    assert entries.columns[0] == 'hpc_x'
    assert entries.columns[1] == 'hpc_y'
    assert entries.iloc[0]['hpc_x'] == -390.1458
    assert entries.iloc[-1]['hpc_x'] == 905.82


def test_get_ars_for_day():
    ar_data = flare_mission_sim.load_ar_data()
    num, noaa, hpc_x, hpc_y, mcintosh, hale, numspots = flare_mission_sim.get_ars_for_day('2011-02-15', ar_data.shift(1))
    assert num == 4
    assert noaa[0] == 11159
    assert mcintosh[0] == 'CSO'
    assert numspots[0] == '3'

    num, noaa, hpc_x, hpc_y, mcintosh, hale, numspots = flare_mission_sim.get_ars_for_day('2011-02-15', ar_data.shift(0))
    assert num == 4
    assert noaa[0] == 11157
    assert mcintosh[0] == 'DSO'
    assert numspots[0] == '9'


def test_get_pointing():

    eclipse, saa, contacts, polar = flare_mission_sim.read_orbit_events()
    eclipse = flare_mission_sim.shift_by_solar_cycle(eclipse, n=-1)
    saa = flare_mission_sim.shift_by_solar_cycle(saa, n=-1)
    contacts = flare_mission_sim.shift_by_solar_cycle(contacts, n=-1)
    eclipse_and_saa = flare_mission_sim.concatenate_range_series(eclipse, saa)
    
    ar_data = flare_mission_sim.load_ar_data()
    ar_targets = flare_mission_sim.generate_ar_targets(ar_data[fm.PHASE_E[0]: fm.PHASE_E[1]],
                                    target_method = 'mcintosh')

    assert ar_targets.loc['2011-02-03']['ar_target'].values == 11150
    assert ar_targets.loc['2011-02-04']['ar_target'].values == 11152
    assert ar_targets.loc['2011-02-05']['ar_target'].values == 11152

    # generate pointing commands at exactly the same time as the decisions for testing purposes
    pointing_commands = flare_mission_sim.generate_pointing_commands(ar_targets, contacts, monte_carlo=False,
                                                        use_random_delay = True,
                                                        delay_params = [0.0, 0.0, 0.0])

    # check that an initial target is correct
    pos = flare_mission_sim.get_pointing('2011-02-04 00:00',pointing_commands,ar_data)
    assert pos[0] == 203.9568
    assert pos[1] == -268.7238
    assert pos[2] == 11150

    # check that the target changed correctly
    pos = flare_mission_sim.get_pointing('2011-02-04 01:00',pointing_commands,ar_data)
    assert_almost_equal(pos[0],-169.403,decimal=2) 
    assert_almost_equal(pos[1],-201.673,decimal=2)
    assert pos[2] == 11152

    # test that the pointing interpolates the position of the AR target over time
    pos = flare_mission_sim.get_pointing('2011-02-04 10:00',pointing_commands,ar_data)
    assert_almost_equal(pos[0],-96.7636,decimal=2) 
    assert_almost_equal(pos[1],-200.5707,decimal=2)
    assert pos[2] == 11152


def test_mcintosh_productivity():
    ar_data = flare_mission_sim.load_ar_data()
    num, noaa, hpc_x, hpc_y, mcintosh, hale, numspots = flare_mission_sim.get_ars_for_day('2011-02-15', ar_data.shift(1))
    ar_rank = [util.mcintosh_class_productivity(this_mcintosh) for this_mcintosh in mcintosh]
    assert np.argmax(ar_rank) == 3
    assert np.max(ar_rank) == 0.73
    assert noaa[np.argmax(ar_rank)] == 11158

    num, noaa, hpc_x, hpc_y, mcintosh, hale, numspots = flare_mission_sim.get_ars_for_day('2011-02-15', ar_data.shift(0))
    ar_rank = [util.mcintosh_class_productivity(this_mcintosh) for this_mcintosh in mcintosh]
    assert np.argmax(ar_rank) == 1
    assert np.max(ar_rank) == 1.0
    assert noaa[np.argmax(ar_rank)] == 11158


def test_hale_productivity():
    allowed_labels = util.find_unique_hale_classes()
    ar_data = flare_mission_sim.load_ar_data()
    num, noaa, hpc_x, hpc_y, mcintosh, hale, numspots = flare_mission_sim.get_ars_for_day('2011-02-15', ar_data.shift(1))
    ar_rank = [util.hale_class_productivity(this_hale, allowed_labels) for this_hale in hale]
    assert np.argmax(ar_rank) == 0
    assert np.argmin(ar_rank) == 1
    assert hale[0] == 'BETA-GAMMA'
    assert hale[1] == 'ALPHA'

    # now apply the tiebreaker (number of spots)
    ar_rank = list(np.array(ar_rank) + np.array(numspots, dtype='float'))
    assert np.argmax(ar_rank) == 3
    assert hale[np.argmax(ar_rank)] == 'BETA-GAMMA'
    assert noaa[np.argmax(ar_rank)] == 11158

    
def test_flare_index():

    flare_index = util.generate_flare_index()
    flare_index_entry = flare_index.loc['2011-03-14']
    assert np.argmax(flare_index_entry['flare_index']) == 1
    assert_almost_equal(max(flare_index_entry['flare_index']),5.44e-05)
    
    ar_data = flare_mission_sim.load_ar_data()
    num, noaa, hpc_x, hpc_y, mcintosh, hale, numspots = flare_mission_sim.get_ars_for_day('2011-03-14', ar_data.shift(0))
    flare_index_entry = flare_index.loc['2011-03-14']
    ar_rank = [util.noaanum_to_flare_index_rank(this_noaa, flare_index_entry) for this_noaa in noaa]

    assert np.argmax(ar_rank) == 2
    assert noaa[np.argmax(ar_rank)] == 11169
    
