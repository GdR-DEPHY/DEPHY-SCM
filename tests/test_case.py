#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Unit tests for :mod:`dephycf.Case`.

These tests cover the :class:`~dephycf.Case.Case` class: case
construction, global attribute management, adding initial-state and
forcing variables, computing derived variables, vertical extension,
completing a case with :meth:`~dephycf.Case.Case.add_missing_init_variables`
/ :meth:`~dephycf.Case.Case.add_missing_forcing_variables`, netCDF
read/write, and time/vertical interpolation.

As for :mod:`dephycf.Variable`, ``Case`` depends on modules that are
not part of the original deliverable (``attributes``, ``thermo``,
``constants``, in addition to ``Axis``, ``Variable``,
``variables_attributes`` and ``plotbasics``). Minimal, documented
stand-ins for all of them are shipped in the ``dephycf`` package next
to these tests; see ``README_tests.md`` for details and for an
important caveat about the ``thermo`` stub (it implements simplified,
NOT scientifically accurate, textbook formulas).

Run the tests (from the project root, next to this ``tests/``
directory and the ``dephycf/`` package):

    pytest -v
"""

import netCDF4 as nc
import numpy as np
import pytest

from dephycf.Case import Case
from dephycf import thermo


# ---------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------

@pytest.fixture
def basic_case():
    """A freshly constructed Case, with latitude/longitude set."""
    return Case('test/basic', lat=45.0, lon=5.0,
                startDate='2020-01-01 00:00:00',
                endDate='2020-01-01 06:00:00')


@pytest.fixture
def level_data_m():
    """A simple 5-level altitude axis (in meters)."""
    return [0., 500., 1000., 1500., 2000.]


@pytest.fixture
def case_with_init_state(basic_case, level_data_m):
    """A Case with a minimal, height-based, initial state.

    Includes ``ta`` (temperature), ``ua``/``va`` (wind, needed to
    infer the number of levels in
    :meth:`~dephycf.Case.Case.add_missing_init_variables`), ``ps``
    (surface pressure) and ``rv`` (water vapor mixing ratio).
    """
    c = basic_case
    nlev = len(level_data_m)

    ta_data = np.array([288., 280., 270., 260., 250.])
    c.add_variable('ta', ta_data.reshape(1, nlev),
                   name='air_temperature', units='K',
                   time=c.t0Axis, lev=level_data_m, levtype='altitude')
    c.add_variable('ua', np.zeros((1, nlev)),
                   name='eastward_wind', units='m s-1',
                   time=c.t0Axis, lev=level_data_m, levtype='altitude')
    c.add_variable('va', np.zeros((1, nlev)),
                   name='northward_wind', units='m s-1',
                   time=c.t0Axis, lev=level_data_m, levtype='altitude')
    c.add_init_ps(101325.)

    rv_data = np.array([0.01, 0.008, 0.006, 0.004, 0.002])
    c.add_variable('rv', rv_data.reshape(1, nlev),
                   name='water_vapor_mixing_ratio', units='kg kg-1',
                   time=c.t0Axis, lev=level_data_m, levtype='altitude')
    return c


# ---------------------------------------------------------------------
# Construction and global attributes
# ---------------------------------------------------------------------

class TestConstruction:

    def test_id_is_split_into_case_and_subcase(self):
        c = Case('DIURNAL/control', startDate='2020-01-01 00:00:00',
                 endDate='2020-01-01 06:00:00')
        assert c.id == 'DIURNAL/control'
        assert c._case == 'DIURNAL'
        assert c._subcase == 'control'

    def test_dates_accept_datetime_objects(self):
        from datetime import datetime
        start = datetime(2020, 6, 1, 0, 0, 0)
        end = datetime(2020, 6, 2, 0, 0, 0)
        c = Case('test/case', startDate=start, endDate=end)
        assert c.start_date == start
        assert c.end_date == end

    def test_dates_accept_14_digit_strings(self):
        c = Case('test/case', startDate='20200101000000',
                 endDate='20200101060000')
        assert c.start_date.year == 2020
        assert c.start_date.hour == 0
        assert c.end_date.hour == 6

    def test_dates_accept_iso_like_strings(self):
        c = Case('test/case', startDate='2020-01-01 00:00:00',
                 endDate='2020-01-01 06:00:00')
        assert c.tstart == 0.0
        assert c.tend == 6 * 3600.0

    def test_lat_lon_add_forcing_variables(self, basic_case):
        assert 'lat' in basic_case.var_forcing_list
        assert 'lon' in basic_case.var_forcing_list
        assert basic_case.lat == 45.0
        assert basic_case.lon == 5.0

    def test_orography_is_always_added(self):
        c = Case('test/case', startDate='2020-01-01 00:00:00',
                 endDate='2020-01-01 06:00:00', zorog=123.)
        assert 'orog' in c.var_forcing_list
        assert np.allclose(c.variables['orog'].data, 123.)

    def test_set_latlon_updates_case_and_variables(self, basic_case):
        basic_case.set_latlon(10.0, 20.0)
        assert basic_case.lat == 10.0
        assert basic_case.lon == 20.0
        assert np.allclose(basic_case.variables['lat'].data, 10.0)
        assert np.allclose(basic_case.variables['lon'].data, 20.0)


class TestAttributes:

    def test_set_attribute_known(self, basic_case):
        basic_case.set_attribute('case_type', 'custom')
        assert basic_case.attributes['case_type'] == 'custom'

    def test_set_attribute_unknown_still_sets_value(self, basic_case,
                                                      caplog):
        # An unknown attribute triggers a warning but is still set,
        # so that a caller relying on set_attribute is not silently
        # ignored.
        basic_case.set_attribute('my_custom_attribute', 42)
        assert basic_case.attributes['my_custom_attribute'] == 42
        assert 'my_custom_attribute' in basic_case.attlist

    def test_title_comment_reference_author(self, basic_case):
        basic_case.set_title('A title')
        basic_case.set_comment('A comment')
        basic_case.set_reference('A reference')
        basic_case.set_author('An author')
        basic_case.set_modifications('Some modifications')
        basic_case.set_script('build_case.py')
        basic_case.set_case_type('standard')

        assert basic_case.attributes['title'] == 'A title'
        assert basic_case.attributes['comment'] == 'A comment'
        assert basic_case.attributes['reference'] == 'A reference'
        assert basic_case.attributes['author'] == 'An author'
        assert basic_case.attributes['modifications'] == 'Some modifications'
        assert basic_case.attributes['script'] == 'build_case.py'
        assert basic_case.attributes['case_type'] == 'standard'


# ---------------------------------------------------------------------
# Adding / removing variables
# ---------------------------------------------------------------------

class TestVariableManagement:

    def test_add_variable_altitude(self, basic_case, level_data_m):
        nlev = len(level_data_m)
        data = np.zeros((1, nlev))
        basic_case.add_variable('ta', data, name='air_temperature',
                                units='K', time=basic_case.t0Axis,
                                lev=level_data_m, levtype='altitude')
        assert 'ta' in basic_case.var_init_list
        assert basic_case.variables['ta'].height is not None
        assert basic_case.variables['ta'].pressure is None

    def test_add_variable_pressure(self, basic_case):
        levdata_pa = [100000., 90000., 80000.]
        data = np.zeros((1, 3))
        basic_case.add_variable('ta', data, name='air_temperature',
                                units='K', time=basic_case.t0Axis,
                                lev=levdata_pa, levtype='pressure')
        assert basic_case.variables['ta'].pressure is not None
        assert basic_case.variables['ta'].height is None

    def test_add_variable_without_time_raises(self, basic_case):
        with pytest.raises(ValueError):
            basic_case.add_variable('ta', np.zeros((3,)), name='ta',
                                    units='K', time=None)

    def test_add_init_ps(self, basic_case):
        basic_case.add_init_ps(101325.)
        assert 'ps' in basic_case.var_init_list
        assert np.allclose(basic_case.variables['ps'].data, 101325.)

    def test_remove_variable(self, case_with_init_state):
        assert 'rv' in case_with_init_state.var_init_list
        case_with_init_state.remove_variable('rv')
        assert 'rv' not in case_with_init_state.var_init_list
        assert 'rv' not in case_with_init_state.variables


# ---------------------------------------------------------------------
# Case information
# ---------------------------------------------------------------------

class TestInfo:

    def test_info_prints_case_and_variable_information(
            self, case_with_init_state, capsys):
        case_with_init_state.info()
        captured = capsys.readouterr()
        assert 'case_type' in captured.out
        assert 'ta' in captured.out


# ---------------------------------------------------------------------
# Computing derived variables
# ---------------------------------------------------------------------

class TestComputeMethods:

    def test_compute_theta_matches_thermo(self, case_with_init_state,
                                           level_data_m):
        # ta has a height companion here, so theta is computed via
        # add_missing_init_variables's derived pressure; we exercise
        # compute_theta() directly using thermo.z2p for the pressure.
        ta = case_with_init_state.variables['ta']
        pressure_data = thermo.z2p(
            z=np.array(level_data_m), ps=101325.,
            ta=ta.data[0,:]).reshape(1, -1)
        from dephycf.Variable import Variable
        pressure_var = Variable('pa', data=pressure_data, name='air_pressure',
                                units='Pa', level=ta.level, time=ta.time)

        theta = case_with_init_state.compute_theta(pressure=pressure_var)

        expected = thermo.t2theta(pressure_data[0], ta.data[0])
        assert np.allclose(theta, expected)

    def test_compute_qv_from_rv(self, case_with_init_state):
        qv = case_with_init_state.compute_qv()
        rv = case_with_init_state.variables['rv'].data[0, :]
        expected = thermo.rt2qt(rv)
        assert np.allclose(qv, expected)

    def test_compute_qt_assumes_equal_to_qv_when_available(self,
                                                              case_with_init_state):
        qv = case_with_init_state.compute_qv()
        case_with_init_state.add_variable(
            'qv', qv.reshape(1, -1), name='specific_humidity',
            units='kg kg-1', time=case_with_init_state.t0Axis,
            lev=case_with_init_state.variables['rv'].level)
        qt = case_with_init_state.compute_qt()
        assert np.allclose(qt, qv)

    def test_compute_theta_raises_without_ta_or_thetal(self, basic_case):
        with pytest.raises(ValueError):
            basic_case.compute_theta(pressure=None)


# ---------------------------------------------------------------------
# Vertical extension
# ---------------------------------------------------------------------

class TestExtendVariable:

    def test_extend_init_temp_adds_one_level(self, case_with_init_state):
        before_shape = case_with_init_state.variables['ta'].data.shape
        case_with_init_state.extend_init_temp(temp=230., height=3000.)
        after_shape = case_with_init_state.variables['ta'].data.shape

        assert after_shape[1] == before_shape[1] + 1
        assert case_with_init_state.variables['ta'].data[0, -1] == 230.
        assert case_with_init_state.variables['ta'].height.data[0, -1] == 3000.


# ---------------------------------------------------------------------
# Completing a case (add_missing_*)
# ---------------------------------------------------------------------

class TestAddMissingVariables:

    def test_add_missing_init_variables_completes_the_case(
            self, case_with_init_state):
        case_with_init_state.add_missing_init_variables()

        for var in ('pa', 'zh', 'theta', 'thetal', 'qv', 'qt', 'rt', 'tke'):
            assert var in case_with_init_state.var_init_list

    def test_add_missing_forcing_variables_noop_without_forcing(
            self, case_with_init_state, caplog):
        # With no forcing variable defined, this should simply warn
        # and return without raising.
        case_with_init_state.var_forcing_list = []
        case_with_init_state.add_missing_forcing_variables()
        assert case_with_init_state.var_forcing_list == []


# ---------------------------------------------------------------------
# netCDF read/write, using real netCDF4 files
# ---------------------------------------------------------------------

class TestNetCDFReadWrite:

    def test_write_then_read_round_trip(self, case_with_init_state,
                                         tmp_path):
        case_with_init_state.add_missing_init_variables()
        # No forcing variables in this minimal case: this call is a
        # documented no-op (see add_missing_forcing_variables).
        case_with_init_state.var_forcing_list = []
        case_with_init_state.add_missing_forcing_variables()

        path = str(tmp_path / "case.nc")
        case_with_init_state.write(path)

        reread = Case('test/basic', startDate='2020-01-01 00:00:00',
                     endDate='2020-01-01 06:00:00')
        reread.read(path)

        assert np.allclose(reread.variables['ta'].data,
                            case_with_init_state.variables['ta'].data)
        assert np.allclose(reread.variables['qv'].data,
                            case_with_init_state.variables['qv'].data)
        assert reread.attributes['case_type'] == \
            case_with_init_state.attributes['case_type']


# ---------------------------------------------------------------------
# Interpolation and SCM conversion
# ---------------------------------------------------------------------

class TestInterpolateAndConvert:

    def test_interpolate_without_time_or_level_preserves_data(
            self, case_with_init_state):
        newcase = case_with_init_state.interpolate(time=None, lev=None)

        assert np.allclose(newcase.variables['ta'].data,
                            case_with_init_state.variables['ta'].data)
        assert 'ta' in newcase.var_init_list

    def test_convert2scm_runs_end_to_end(self, case_with_init_state):
        # This minimal case has no forcing variable, so
        # add_missing_forcing_variables (called internally by
        # convert2SCM) takes its documented no-op branch.
        case_with_init_state.var_forcing_list = []

        scm_case = case_with_init_state.convert2SCM()

        for var in ('ta', 'pa', 'zh', 'theta', 'qv'):
            assert var in scm_case.var_init_list
