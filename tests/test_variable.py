#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Unit tests for :mod:`dephycf.Variable`.

These tests cover the :class:`~dephycf.Variable.Variable` class as
well as the module-level functions :func:`~dephycf.Variable.read` and
:func:`~dephycf.Variable.interpol`.

The ``Variable`` module depends on ``Axis``, ``plotbasics`` and
``variables_attributes``. To avoid depending on any real project
data, this test file relies on:

- the ``dephycf`` package shipped alongside these tests, which
  contains a minimal ("stub") implementation of these three modules,
  sufficient to faithfully exercise ``Variable``'s logic;
- real, temporary netCDF files (created with ``netCDF4`` in a
  ``pytest`` ``tmp_path``) to exercise ``Variable.write()`` and the
  module-level ``read()`` function, instead of fake/mocked objects.

Run the tests (from the project root, next to this ``tests/``
directory and the ``dephycf/`` package):

    pytest

or, for a more detailed report:

    pytest -v
"""

import netCDF4 as nc
import numpy as np
import pytest

from dephycf.Variable import Variable, read, interpol
from dephycf.Axis import Axis


# ---------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------

@pytest.fixture
def time_axis():
    """Simple 3-timestep time axis."""
    return Axis('time', [0, 3600, 7200], units='seconds since 2020-01-01')


@pytest.fixture
def level_axis():
    """Simple 3-level (pressure) vertical axis."""
    return Axis('lev', [1000., 850., 700.], units='hPa')


@pytest.fixture
def sample_data():
    """3x3 (time, level) data array."""
    return np.arange(9, dtype=float).reshape(3, 3)


@pytest.fixture
def empty_netcdf_dataset(tmp_path):
    """A real, empty, writable netCDF4.Dataset in a temporary file.

    Used to test ``Variable.write()`` against an actual netCDF file
    instead of a fake/mocked one.
    """
    path = tmp_path / "output.nc"
    with nc.Dataset(path, mode="w") as ds:
        yield ds


@pytest.fixture
def sample_netcdf_dataset(tmp_path, sample_data):
    """A real, read-only netCDF4.Dataset with a sample variable.

    Contains a ``time`` dimension/variable, a ``lev`` dimension/
    variable, a ``zh_ta`` height companion variable, and a ``ta``
    variable with a ``coordinates`` attribute pointing to it. Used to
    test the module-level ``read()`` function against an actual
    netCDF file.
    """
    path = tmp_path / "sample.nc"
    with nc.Dataset(path, mode="w") as ds:
        ds.createDimension('time', 3)
        ds.createDimension('lev', 3)

        time_var = ds.createVariable('time', 'f8', ('time',))
        time_var[:] = [0., 1., 2.]
        time_var.standard_name = 'time'
        time_var.units = 'seconds since 2020-01-01'

        lev_var = ds.createVariable('lev', 'f8', ('lev',))
        lev_var[:] = [1000., 850., 700.]
        lev_var.standard_name = 'pressure'
        lev_var.units = 'hPa'

        height_var = ds.createVariable('zh_ta', 'f8', ('time', 'lev'))
        height_var[:] = np.tile([0., 100., 200.], (3, 1))
        height_var.units = 'm'

        ta_var = ds.createVariable('ta', 'f8', ('time', 'lev'))
        ta_var[:] = sample_data
        ta_var.standard_name = 'air_temperature'
        ta_var.units = 'K'
        ta_var.coordinates = 'time zh_ta lat lon'

    with nc.Dataset(path, mode="r") as ds:
        yield ds


@pytest.fixture
def sample_netcdf_dataset_no_standard_name(tmp_path, sample_data):
    """Same as `sample_netcdf_dataset`, but without a standard_name
    attribute on the ``ta`` variable, to exercise the fallback branch
    of ``read()``.
    """
    path = tmp_path / "sample_no_standard_name.nc"
    with nc.Dataset(path, mode="w") as ds:
        ds.createDimension('time', 3)
        ds.createDimension('lev', 3)

        time_var = ds.createVariable('time', 'f8', ('time',))
        time_var[:] = [0., 1., 2.]
        time_var.standard_name = 'time'
        time_var.units = 's'

        lev_var = ds.createVariable('lev', 'f8', ('lev',))
        lev_var[:] = [1000., 850., 700.]
        lev_var.standard_name = 'pressure'
        lev_var.units = 'hPa'

        height_var = ds.createVariable('zh_ta', 'f8', ('time', 'lev'))
        height_var[:] = np.tile([0., 100., 200.], (3, 1))
        height_var.units = 'm'

        ta_var = ds.createVariable('ta', 'f8', ('time', 'lev'))
        ta_var[:] = sample_data
        # No standard_name attribute set on purpose.
        ta_var.units = 'K'
        ta_var.coordinates = 'time zh_ta lat lon'

    with nc.Dataset(path, mode="r") as ds:
        yield ds


# ---------------------------------------------------------------------
# Construction
# ---------------------------------------------------------------------

class TestConstruction:

    def test_data_is_converted_to_float64(self, time_axis, level_axis):
        v = Variable('ta', data=np.zeros((3, 3), dtype=np.int32),
                     time=time_axis, level=level_axis)
        assert v.data.dtype == np.float64

    def test_axes_and_coord_from_time_level_kwargs(self, time_axis,
                                                     level_axis,
                                                     sample_data):
        v = Variable('ta', data=sample_data, time=time_axis,
                     level=level_axis)
        assert v.axlist == ['time', 'lev']
        assert v.axes == [time_axis, level_axis]
        assert v.coord == 'time lev lat lon'

    def test_axes_list_identifies_time_and_level(self, sample_data):
        t0 = Axis('t0', [0])
        lev = Axis('nlev', [1000.])
        v = Variable('ta', data=sample_data, axes=[t0, lev],
                     axlist=['t0', 'nlev'])
        assert v.time is t0
        assert v.level is lev

    def test_unexpected_axis_raises_value_error(self, sample_data):
        weird = Axis('weird_axis', [1, 2, 3])
        with pytest.raises(ValueError):
            Variable('ta', data=sample_data, axes=[weird],
                     axlist=['weird_axis'])

    def test_plotcoef_and_plotunits_from_var_attributes(self, time_axis,
                                                          level_axis,
                                                          sample_data):
        # 'ta' is defined in the variables_attributes stub with
        # plotcoef=1.0 and plotunits='K'.
        v = Variable('ta', data=sample_data, units='K',
                     time=time_axis, level=level_axis)
        assert v.plotcoef == 1.0
        assert v.plotunits == 'K'

    def test_plotcoef_and_plotunits_fallback(self, time_axis, level_axis,
                                              sample_data):
        # 'unknown_var' is not in variables_attributes: fall back to
        # the default / explicitly given values.
        v = Variable('unknown_var', data=sample_data, units='m/s',
                     time=time_axis, level=level_axis, plotcoef=2.5)
        assert v.plotcoef == 2.5
        assert v.plotunits == 'm/s'

    def test_lat_lon_are_stored_on_the_instance(self, time_axis,
                                                 level_axis, sample_data):
        v = Variable('ta', data=sample_data, time=time_axis,
                     level=level_axis, lat=48.8, lon=2.3)
        assert v.lat == 48.8
        assert v.lon == 2.3

    def test_lat_lon_default_to_none(self, time_axis, level_axis,
                                      sample_data):
        v = Variable('ta', data=sample_data, time=time_axis,
                     level=level_axis)
        assert v.lat is None
        assert v.lon is None

    def test_height_companion_created_from_raw_array(self, time_axis,
                                                       level_axis,
                                                       sample_data):
        height = np.tile([0., 100., 200.], (3, 1))
        v = Variable('ta', data=sample_data, time=time_axis,
                     level=level_axis, height=height, height_units='m')
        assert v.height is not None
        assert v.height.id == 'zh_ta'
        assert v.height.units == 'm'
        assert v.coord == 'time zh_ta lat lon'

    def test_height_companion_created_from_axis_instance(self, time_axis,
                                                           level_axis,
                                                           sample_data):
        height_axis = Axis('zh_custom', [0., 100., 200.], units='m',
                            name='height above surface')
        v = Variable('ta', data=sample_data, time=time_axis,
                     level=level_axis, height=height_axis)
        assert v.height.id == 'zh_custom'
        assert v.coord == 'time zh_custom lat lon'

    def test_pressure_companion_used_when_no_height(self, time_axis,
                                                     level_axis,
                                                     sample_data):
        pressure = np.tile([101325., 90000., 70000.], (3, 1))
        v = Variable('ta', data=sample_data, time=time_axis,
                     level=level_axis, pressure=pressure,
                     pressure_units='Pa')
        assert v.height is None
        assert v.pressure is not None
        assert v.pressure.id == 'pa_ta'
        assert v.pressure.units == 'Pa'

    def test_height_is_privileged_over_pressure(self, time_axis,
                                                 level_axis, sample_data):
        height = np.tile([0., 100., 200.], (3, 1))
        pressure = np.tile([101325., 90000., 70000.], (3, 1))
        v = Variable('ta', data=sample_data, time=time_axis,
                     level=level_axis, height=height, pressure=pressure)
        assert v.height is not None
        assert v.pressure is None


# ---------------------------------------------------------------------
# Utility methods: info, set_coordinates, set_level
# ---------------------------------------------------------------------

class TestUtilityMethods:

    def test_info_prints_summary(self, time_axis, level_axis, sample_data,
                                  capsys):
        v = Variable('ta', data=sample_data, name='air_temperature',
                     units='K', time=time_axis, level=level_axis)
        v.info()
        captured = capsys.readouterr()
        assert 'ta' in captured.out
        assert 'air_temperature' in captured.out

    def test_set_coordinates(self, time_axis, level_axis, sample_data):
        v = Variable('ta', data=sample_data, time=time_axis,
                     level=level_axis)
        v.set_coordinates('foo', 'bar', 'baz')
        assert v.coord == 'foo bar baz'

    def test_set_level_replaces_level_axis(self, time_axis, level_axis,
                                            sample_data):
        v = Variable('ta', data=sample_data, time=time_axis,
                     level=level_axis)
        new_level = Axis('lev2', [900., 800.], units='hPa')
        v.set_level(new_level)
        assert v.level is new_level
        assert v.axlist == ['time', 'lev2']
        assert v.axes == [time_axis, new_level]

    def test_set_level_with_none_logs_warning_and_does_nothing(
            self, time_axis, level_axis, sample_data, caplog):
        v = Variable('ta', data=sample_data, time=time_axis,
                     level=level_axis)
        original_axlist = list(v.axlist)
        v.set_level(None)
        assert v.axlist == original_axlist


# ---------------------------------------------------------------------
# netCDF read/write, using real netCDF4 files
# ---------------------------------------------------------------------

class TestNetCDFReadWrite:

    def test_write_creates_variable_with_expected_attributes(
            self, time_axis, level_axis, sample_data,
            empty_netcdf_dataset):
        v = Variable('ta', data=sample_data, name='air_temperature',
                     units='K', time=time_axis, level=level_axis)

        v.write(empty_netcdf_dataset)

        assert 'ta' in empty_netcdf_dataset.variables
        written = empty_netcdf_dataset.variables['ta']
        assert np.allclose(written[:], sample_data)
        assert written.standard_name == 'air_temperature'
        assert written.units == 'K'
        assert written.coordinates == v.coord
        # The time/level dimensions and coordinate variables should
        # also have been written by Variable.write().
        assert 'time' in empty_netcdf_dataset.dimensions
        assert 'lev' in empty_netcdf_dataset.dimensions

    def test_write_does_not_overwrite_existing_variable(
            self, time_axis, level_axis, sample_data,
            empty_netcdf_dataset):
        v = Variable('ta', data=sample_data, name='air_temperature',
                     units='K', time=time_axis, level=level_axis)
        v.write(empty_netcdf_dataset)

        other_data = sample_data + 100.
        v2 = Variable('ta', data=other_data, name='air_temperature',
                     units='K', time=time_axis, level=level_axis)
        v2.write(empty_netcdf_dataset)  # second write: should be a no-op

        assert np.allclose(
            empty_netcdf_dataset.variables['ta'][:], sample_data)

    def test_read_builds_variable_with_height_companion(
            self, sample_netcdf_dataset, sample_data):
        v = read('ta', sample_netcdf_dataset)

        assert v.id == 'ta'
        assert v.name == 'air_temperature'
        assert np.allclose(v.data, sample_data)
        assert v.height is not None
        assert v.height.id == 'zh_ta'

    def test_read_falls_back_when_standard_name_missing(
            self, sample_netcdf_dataset_no_standard_name):
        # Simulates a netCDF variable without a standard_name
        # attribute: read()'s "except AttributeError" branch should
        # be taken.
        v = read('ta', sample_netcdf_dataset_no_standard_name)

        assert v.id == 'ta'
        assert v.name is None


# ---------------------------------------------------------------------
# Time and vertical interpolation
# ---------------------------------------------------------------------

class TestInterpolation:

    def test_interpol_time_on_time_only_variable(self):
        t_in = Axis('time', [0, 1, 2], units='s')
        v = Variable('x', data=np.array([10., 20., 30.]), time=t_in)
        t_out = Axis('time', [0.5, 1.5], units='s')

        v_interp = v.interpol_time(time=t_out)

        assert np.allclose(v_interp.data, [15., 25.])

    def test_interpol_time_on_time_level_variable(self, time_axis,
                                                    level_axis):
        data = np.array([[0., 0., 0.],
                          [10., 20., 30.],
                          [20., 40., 60.]])
        v = Variable('ta', data=data, time=time_axis, level=level_axis)
        t_out = Axis('time', [1800, 5400], units='seconds since 2020-01-01')

        v_interp = v.interpol_time(time=t_out)

        assert v_interp.data.shape == (2, 3)
        assert np.allclose(v_interp.data[0], [5., 10., 15.])

    def test_interpol_time_returns_self_for_initial_state_variable(
            self, time_axis, level_axis):
        data = np.zeros((1, 3))
        v = Variable('ta', data=data, time=time_axis, level=level_axis)
        t_out = Axis('time', [0.5], units='s')

        assert v.interpol_time(time=t_out) is v

    def test_interpol_time_returns_self_when_no_target_time(
            self, time_axis, level_axis, sample_data):
        v = Variable('ta', data=sample_data, time=time_axis,
                     level=level_axis)
        assert v.interpol_time(time=None) is v

    def test_interpol_time_without_time_axis_raises(self, level_axis):
        v = Variable('ta', data=np.array([1., 2., 3.]), level=level_axis)
        with pytest.raises(ValueError):
            v.interpol_time(time=Axis('time', [0, 1], units='s'))

    def test_interpol_vert_with_height(self, time_axis, level_axis):
        data = np.array([[0., 10., 20.],
                          [0., 10., 20.],
                          [0., 10., 20.]])
        height = np.tile([0., 100., 200.], (3, 1))
        v = Variable('ta', data=data, time=time_axis, level=level_axis,
                     height=height, height_units='m')

        target_height = np.array([50., 150.])
        v_interp = v.interpol_vert(height=target_height)

        assert v_interp.data.shape == (3, 2)
        assert np.allclose(v_interp.data[0], [5., 15.])

    def test_interpol_vert_returns_self_without_level_axis(self):
        v = Variable('surface_var', data=np.array([1., 2., 3.]))
        result = v.interpol_vert(height=np.array([10.]))
        assert result is v

    def test_interpol_vert_returns_self_without_target(self, time_axis,
                                                          level_axis,
                                                          sample_data):
        v = Variable('ta', data=sample_data, time=time_axis,
                     level=level_axis)
        assert v.interpol_vert() is v


# ---------------------------------------------------------------------
# Module-level interpol() function
# ---------------------------------------------------------------------

class TestModuleInterpol:

    def test_interpol_raises_when_no_time_nor_level(self):
        v = Variable('scalar_var', data=np.array([1.]))
        with pytest.raises(ValueError):
            interpol(v)

    def test_interpol_time_and_level_branch(self, time_axis, level_axis):
        # Now that Variable stores self.lat/self.lon, interpol() can
        # read var.lat/var.lon without raising an AttributeError.
        data4d = np.zeros((3, 3, 1, 1))
        for it in range(3):
            data4d[it, :, 0, 0] = [10. * it, 20. * it, 30. * it]
        v = Variable('ta', data=data4d, time=time_axis, level=level_axis,
                     lat=45.0, lon=5.0)
        levout = Axis('lev', [900.], units='hPa')

        vout = interpol(v, levout=levout)

        assert vout.data.shape == (3, 1, 1, 1)
        assert vout.lat == 45.0
        assert vout.lon == 5.0
