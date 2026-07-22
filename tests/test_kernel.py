import numpy as np

import pytest

from dephycf.kernel import (
    gaussian_weights,
    gaussian_kernel_mean,
)


def test_weights_mean():
    x = np.arange(10)

    w = gaussian_weights(x, 5, 2)

    assert np.isclose(w.mean(), 1.0)


def test_weights_positive():

    w = gaussian_weights(np.arange(20), 10, 3)

    assert np.all(w > 0)


def test_kernel_constant():

    x = np.arange(20)

    y = np.ones(20) * 7.

    out = gaussian_kernel_mean(y, x, sigma=2.)

    np.testing.assert_allclose(out, 7.0)


def test_kernel_linear():

    x = np.linspace(0, 10, 100)

    y = x

    out = gaussian_kernel_mean(y, x, sigma=0.2)

    np.testing.assert_allclose(
        out[20:-20],
        y[20:-20],
        atol=5e-2,
    )


def test_output_coordinates():

    x = np.arange(100)

    y = np.sin(x)

    xo = np.linspace(0, 99, 20)

    out = gaussian_kernel_mean(y, x, 5, xo)

    assert out.shape == xo.shape


def test_sigma_error():

    with pytest.raises(ValueError):
        gaussian_weights(np.arange(10), 5, 0)


def test_shape_error():

    with pytest.raises(ValueError):
        gaussian_kernel_mean(
            np.arange(10),
            np.arange(11),
            1,
        )
