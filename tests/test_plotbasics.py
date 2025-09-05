# tests/test_plotting.py
import os
import numpy as np
import pytest
from dephycf import plotbasics


@pytest.fixture
def tmp_png(tmp_path):
    """Helper fixture that returns a temporary PNG filepath."""
    return tmp_path / "test.png"


def test_plot_creates_file(tmp_png):
    x = np.linspace(0, 10, 100)
    y = np.sin(x)

    plotbasics.plot(x, y, name=str(tmp_png))
    assert tmp_png.exists()
    assert tmp_png.stat().st_size > 0


def test_plot_with_second_curve(tmp_png):
    x = np.linspace(0, 10, 100)
    y = np.sin(x)
    y2 = np.cos(x)

    plotbasics.plot(x, y, x2=x, y2=y2, name=str(tmp_png), label="sin", label2="cos")
    assert tmp_png.exists()


def test_plot2d_creates_file(tmp_png):
    x = np.linspace(0, 2 * np.pi, 50)
    y = np.linspace(0, 1, 20)
    X, Y = np.meshgrid(x, y, indexing="ij")
    Z = np.sin(X) * np.cos(2 * np.pi * Y)

    plotbasics.plot2D(x, y, Z, name=str(tmp_png))
    assert tmp_png.exists()
    assert tmp_png.stat().st_size > 0


def test_plot2d_shape_mismatch(tmp_png):
    x = np.linspace(0, 10, 50)
    y = np.linspace(0, 5, 20)
    z = np.ones((10, 10))  # Wrong shape

    with pytest.raises(TypeError):
        plotbasics.plot2D(x, y, z, name=str(tmp_png))

