import numpy as np
import pytest
from netCDF4 import Dataset
from dephycf.Axis import Axis


def test_axis_initialization():
    """Test initialization of Axis object with basic parameters."""
    axis = Axis("time", [0, 1, 2, 3], name="time", units="seconds")

    assert axis.id == "time"
    assert isinstance(axis.data, np.ndarray)
    assert axis.length == 4
    assert axis.name == "time"
    assert axis.units == "seconds"


def test_axis_with_extra_metadata():
    """Test that arbitrary metadata is stored in the Axis object."""
    axis = Axis("height", [0, 10, 20], name="altitude", units="m", source="radar")

    assert hasattr(axis, "source")
    assert axis.source == "radar"


def test_axis_info_prints_output(capsys):
    """Test that info() prints expected output."""
    axis = Axis("lat", [10, 20, 30], name="latitude", units="degrees_north")

    axis.info(data=True)
    captured = capsys.readouterr()

    assert "axis id: lat" in captured.out
    assert "length: 3" in captured.out
    assert "data:" in captured.out


def test_axis_write_to_netcdf(tmp_path):
    """Test writing axis to a temporary NetCDF file."""
    file_path = tmp_path / "test.nc"

    axis = Axis("lon", [0, 90, 180], name="longitude", units="degrees_east")

    with Dataset(file_path, "w", format="NETCDF4") as nc:
        axis.write(nc)

    # Reopen and check contents
    with Dataset(file_path, "r") as nc:
        assert "lon" in nc.dimensions
        assert "lon" in nc.variables

        var = nc.variables["lon"]
        np.testing.assert_array_equal(var[:], [0, 90, 180])
        assert var.standard_name == "longitude"
        assert var.units == "degrees_east"

