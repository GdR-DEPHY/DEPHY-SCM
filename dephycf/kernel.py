"""
Gaussian kernel smoothing utilities.

This module provides simple functions for computing Gaussian weights and
Gaussian kernel moving averages.

The Gaussian kernel is defined as

.. math::

    w(x) = \\exp\\left(-\\frac{(x-\\mu)^2}{2\\sigma^2}\\right)

where :math:`\\mu` is the kernel center and :math:`\\sigma` its standard
deviation.

The weights are normalized so that their mean equals one, which preserves
the magnitude of weighted averages.
"""

import numpy as np
from scipy.stats import norm


def gaussian_weights(x, center, sigma):
    """
    Compute Gaussian weights.

    Parameters
    ----------
    x : array-like
        Coordinates.
    center : :class:`float`
        Center of the Gaussian kernel.
    sigma : :class:`float`
        Standard deviation of the Gaussian kernel.

    Returns
    -------
    ndarray
        Gaussian weights normalized by their mean.

    Notes
    -----
    The weights are computed as

    .. math::

        w_i =
        \\frac{\\mathcal N(x_i;\\mu,\\sigma)}
             {\\mathrm{mean}(\\mathcal N)}

    where :math:`\\mathcal N` denotes the Gaussian probability density.
    """
    if sigma <= 0:
        raise ValueError("sigma must be strictly positive")

    x = np.asarray(x, dtype=np.float32)

    weights = norm.pdf(x, loc=center, scale=sigma)

    return weights / weights.mean()


def gaussian_kernel_mean(x, y, sigma, yout=None):
    """
    Compute Gaussian kernel moving averages.

    Parameters
    ----------
    x : array-like
        Values to smooth.
    y : array-like
        Coordinates associated with the values.
    sigma : :class:`float`
        Gaussian kernel standard deviation.
    yout : array-like, optional
        Coordinates where the smoothed values are evaluated.
        If None, ``y`` is used.

    Returns
    -------
    ndarray
        Smoothed values.

    Notes
    -----
    The smoothed value at position :math:`y_0` is

    .. math::

        f(y_0)=
        \\frac{\\sum_i w_i(y_0)x_i}
             {\\sum_i w_i(y_0)}

    where :math:`w_i` are Gaussian weights.
    """
    if sigma <= 0:
        raise ValueError("sigma must be strictly positive")

    x = np.asarray(x, dtype=np.float32)
    y = np.asarray(y, dtype=np.float32)

    if x.shape != y.shape:
        raise ValueError(
            "values and coordinates must have identical shapes"
        )

    if yout is None:
        yout = y

    yout = np.asarray(yout, dtype=np.float32)

    return np.asarray(
        [
            np.average(
                x,
                weights=gaussian_weights(
                    y,
                    center,
                    sigma,
                ),
            )
            for center in yout
        ]
    )
