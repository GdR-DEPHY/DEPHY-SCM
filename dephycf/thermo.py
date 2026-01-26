#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Thermodynamical functions and other utilities.

This module provides conversions between water variables
(mixing ratio :math:`r_t`, specific humidity :math:`q_t`, relative humidity :math:`\\mathcal{H}`),
temperature variables (:math:`T`, :math:`\\theta`), and hydrostatic conversions (:math:`p` <-> :math:`z`).

Simple hypotheses are used and might be updated in the future.
"""

import logging
import math
import numpy as np

try:
    from . import constants as cc
except ImportError:
    import constants as cc

logger = logging.getLogger(__name__)

# =============================================================================
# Optional dependency: MetPy
# =============================================================================
try:
    from metpy.calc import (
        mixing_ratio_from_relative_humidity,
        relative_humidity_from_mixing_ratio,
        specific_humidity_from_dewpoint,
    )
    from metpy.units import units as Munits

    HAS_METPY = True
except ImportError:
    HAS_METPY = False
    logger.debug("MetPy not available: RH-based functions disabled")


# =============================================================================
# Helpers
# =============================================================================
def _check_units(units: str, allowed: tuple[str, ...]) -> None:
    """Validate unit string."""
    if units not in allowed:
        raise ValueError(f"Unknown units '{units}', expected one of {allowed}")


def _asarray(x, dtype=np.float64):
    """Convert input to numpy array (float64)."""
    return np.asarray(x, dtype=dtype)


def _init_temp(z, T0=300, lapse_rate=-6.5e-3):
    """
    Initialize a temperature profile starting à T0 (300 K by default) at the surface
    and using a given lapse rate (-6.5 K km-1 by default)
    """

    return T0 + lapse_rate * z

def _init_qv(z, qv0=0.010, Hq=2000.0):
    """
    Initialize a specific humidity profile starting à qv0 (10 g kg-1 by default) at the surface
    and exponentially decreasing with a length scale Hq (2 km by default)
    """

    return qv0 * np.exp(-z / Hq)


# =============================================================================
# Conversion specific mass of total water / total water mixing ratio
# =============================================================================
def rt2qt(rt, units="kg kg-1"):
    """
    Convert total water mixing ratio to specific mass of total water.

    Parameters
    ----------
    rt : :class:`float` or array-like
        Total water mixing ratio.
    units : :class:`str`, optional
        Units for total water mixing ratio.
        Should be 'kg kg-1' (default) or 'g kg-1'.

    Returns
    -------
    :class:`float` or :class:`ndarray`
        Specific mass of total water (in **units**).
    """
    _check_units(units, ("kg kg-1", "g kg-1"))
    rt = _asarray(rt)

    if units == "kg kg-1":
        return rt / (1.0 + rt)

    return rt / (1.0 + rt/1000.)


def qt2rt(qt, units="kg kg-1"):
    """
    Convert specific mass of total water to total water mixing ratio.

    Parameters
    ----------
    qt : :class:`float` or array-like
        Specific mass of total water.
    units : :class:`str`, optional
        Units for total water mixing ratio.
        Should be 'kg kg-1' (default) or 'g kg-1'.

    Returns
    -------
    :class:`float` or :class:`ndarray`
        Total water mixing ratio (in **units**).
    """
    _check_units(units, ("kg kg-1", "g kg-1"))
    qt = _asarray(qt)

    if units == "kg kg-1":
        return qt / (1.0 - qt)

    return qt / (1.0 - qt/1000.)


# ====================================================================================
# Conversion between relative humidity and specific mass / mixing ratio of total water
# (MetPy required)
# ====================================================================================
def hur2rt(hur, pressure, temp, units="kg kg-1"):
    """
    Convert relative humidity to total water mixing ratio.
    Requires the MetPy module.

    Parameters
    ----------
    hur : :class:`float` or array-like
        Relative humidity (0–1).
    pressure : :class:`float` or array-like
        Pressure (Pa).
    temp : :class:`float` or array-like
        Temperature (K).
    units : :class:`str`, optional
        Units for total water mixing ratio.
        Should be 'kg kg-1' (default) or 'g kg-1'.

    Returns
    -------
    :class:`float` or :class:`ndarray`
        Total water mixing ratio.
    """
    if not HAS_METPY:
        raise RuntimeError("MetPy is required for RH-based conversions")

    _check_units(units, ("kg kg-1", "g kg-1"))

    hur = _asarray(hur)
    pressure = _asarray(pressure) * Munits.Pa
    temp = _asarray(temp) * Munits.kelvin

    rv = mixing_ratio_from_relative_humidity(
        pressure, temp, hur
    ).magnitude

    if units == "kg kg-1":
        return rv

    return rv * 1000.0


def hur2qt(hur, pressure, temp, units="kg kg-1"):
    """
    Convert relative humidity to specific mass of total water.
    Requires the MetPy module.

    Parameters
    ----------
    hur : :class:`float` or array-like
        Relative humidity (0–1).
    pressure : :class:`float` or array-like
        Pressure (Pa).
    temp : :class:`float` or array-like
        Temperature (K).
    units : :class:`str`, optional
        Units for the specific mass of total water.
        Should be 'kg kg-1' (default) or 'g kg-1'.

    Returns
    -------
    :class:`float` or :class:`ndarray`
        Specific mass of total water.
    """
    if not HAS_METPY:
        raise RuntimeError("MetPy is required for RH-based conversions")

    _check_units(units, ("kg kg-1", "g kg-1"))

    rt = hur2rt(hur, pressure, temp, units=units)

    if units == "kg kg-1":
        return rt / (1.0 + rt)

    return (rt / (1000.0 + rt))


def rt2hur(rt, pressure, temp, units="kg kg-1"):
    """
    Convert total water mixing ratio to relative humidity.
    Require the MetPy module.

    Parameters
    ----------
    rt : :class:`float` or array-like
        Total water mixing ratio.
    pressure : :class:`float` or array-like
        Pressure (Pa).
    temp : :class:`float` or array-like
        Temperature (K).
    units : :class:`str`, optional
        Units for the total water mixing ratio.
        Should be 'kg kg-1' (default) or 'g kg-1'.

    Returns
    -------
    :class:`float` or :class:`ndarray`
        Relative humidity (0-1).
    """
    if not HAS_METPY:
        raise RuntimeError("MetPy is required for RH-based conversions")

    _check_units(units, ("kg kg-1", "g kg-1"))

    rt = _asarray(rt)
    pressure = _asarray(pressure) * Munits.Pa
    temp = _asarray(temp) * Munits.kelvin

    rv = rt if units == "kg kg-1" else rt / 1000.0

    return relative_humidity_from_mixing_ratio(
        pressure, temp, rv
    ).magnitude


def qt2hur(qt, pressure, temp, units="kg kg-1"):
    """
    Convert specific mass of total water to relative humidity.

    Parameters
    ----------
    qt : :class:`float` or array-like
        Specific mass of total water.
    pressure : :class:`float` or array-like
        Pressure (Pa).
    temp : :class:`float` or array-like
        Temperature (K).
    units : :class:`str`, optional
        Units for the specific mass of total water.
        Should be 'kg kg-1' (default) or 'g kg-1'.

    Returns
    -------
    :class:`float` or :class:`ndarray`
        Relative humidity (0-1).
    """
    if not HAS_METPY:
        raise RuntimeError("MetPy is required for RH-based conversions")

    _check_units(units, ("kg kg-1", "g kg-1"))

    qt = _asarray(qt)
    qt = qt if units == "kg kg-1" else qt / 1000.0
    rt = qt / (1.0 - qt)

    return rt2hur(rt, pressure, temp, units="kg kg-1")
    
# =============================================================================
# Advection conversions
# =============================================================================
def advrt2advqt(rt, advrt, rt_units="kg kg-1"):
    """
    Convert advection of the total water mixing ratio
    to the advection of the specific mass of total water.

    .. math:: \\left.\\frac{\\partial q_t}{\\partial t}\\right|_\\textrm{adv}
              = \\frac{1}{(1+r_t)^2} \\times 
              \\left.\\frac{\\partial r_t}{\\partial t}\\right|_\\textrm{adv}

    Parameters
    ----------
    rt : :class:`float` or array-like
        Total water mixing ratio (in **rt_units**).
    advrt : :class:`float` or array-like
        Advection of total water mixing ratio (kg kg-1 s-1).
    rt_units : :class:`str`, optional
        Units for the total water mixing ratio.
        Should be 'kg kg-1' (default) or 'g kg-1'.

    Returns
    -------
    :class:`float` or :class:`ndarray`
        Advection of the specific mass of total water (kg kg-1 s-1).
    """

    _check_units(rt_units, ("kg kg-1", "g kg-1"))

    rt = _asarray(rt)
    advrt = _asarray(advrt)

    if rt_units == "kg kg-1":
        return advrt / (1.0 + rt) ** 2

    return advrt / (1.0 + rt / 1000.0) ** 2


def advqt2advrt(qt, advqt, qt_units="kg kg-1"):
    """
    Convert advection of the specific mass of total water
    to the advection of the total water mixing ratio, following:

    .. math:: \\left.\\frac{\\partial r_t}{\\partial t}\\right|_\\textrm{adv}
              = \\frac{1}{(1-q_t)^2} \\times 
              \\left.\\frac{\\partial q_t}{\\partial t}\\right|_\\textrm{adv}

    Parameters
    ----------
    qt : :class:`float` or array-like
        Specific mass of total water (in **qt_units**).
    advqt : :class:`float` or array-like
        Advection of the specific mass of total water (kg kg-1 s-1).
    qt_units : :class:`str`, optional
        Units for the specific mass of total water.
        Should be 'kg kg-1' (default) or 'g kg-1'.

    Returns
    -------
    :class:`float` or :class:`ndarray`
        Advection of the specific mass of total water (kg kg-1 s-1).
    """
    _check_units(qt_units, ("kg kg-1", "g kg-1"))

    qt = _asarray(qt)
    advqt = _asarray(advqt)

    if qt_units == "kg kg-1":
        return advqt / (1.0 - qt) ** 2

    return advqt / (1.0 - qt / 1000.0) ** 2


# =============================================================================
# Temperature / potential temperature
# =============================================================================
def theta2t(p, theta, p0=cc.p0, kappa=cc.kappa):
    """
    Convert potential temperature to temperature
    using air pressure according to:

    .. math:: T = \\theta \\left( \\frac{p}{p_0}  \\right)^\\kappa

    Default values for :math:`p_0` and :math:`\\kappa` comes from :mod:`.constants`.

    Parameters
    ----------
    p : :class:`float` or array-like
        Air pressure :math:`p` (Pa).
    theta : :class:`float` or array-like
        Potential temperature :math:`\\theta` (K).
    p0 : :class:`float`, optional
        Reference pressure :math:`p_0` (Pa). Value from :data:`.constants.p0`
    kappa : :class:`float`, optional
        Ratio :math:`\\kappa` of the gas constant to the gas heat capacity (no units).
        Default is for dry air taken from :data:`.constants.kappa`.

    Returns
    -------
    :class:`float` or :class:`ndarray`
        Temperature :math:`T` (K).
    """
    p = _asarray(p)
    theta = _asarray(theta)
    return theta * (p / p0) ** kappa


def t2theta(p, temp, p0: float = cc.p0, kappa: float = cc.kappa):
    """
    Convert temperature to potential temperature
    using air pressure according to:

    .. math:: \\theta = T \\left( \\frac{p_0}{p}  \\right)^\\kappa

    Default values for :math:`p_0` and :math:`\\kappa` comes from :mod:`.constants`.

    Parameters
    ----------
    p : :class:`float` or array-like
        Air pressure :math:`p` (Pa).
    temp : :class:`float` or array-like
        Temperature :math:`T` (K).
    p0 : :class:`float`, optional
        Reference pressure :math:`p_0` (Pa). Value from :data:`.constants.p0`
    kappa : :class:`float`, optional
        Ratio :math:`\\kappa` of the gas constant to the gas heat capacity (no units).
        Default is for dry air taken from :data:`.constants.kappa`.

    Returns
    -------
    :class:`float` or :class:`ndarray`
        Potential temperature :math:`\\theta` (K).
    """

    p = _asarray(p)
    temp = _asarray(temp)
    return temp * (p0 / p) ** kappa


# =============================================================================
# Hydrostatic conversions
# =============================================================================
def z2p(
    z,
    ps,
    theta=None, ta=None, thetal=None,
    qv=None,
    g=cc.g, Rd=cc.Rd, Rv=cc.Rv, p0=cc.p0, kappa=cc.kappa,
):
    """
    Compute the pressure of a set of given altitudes using the hydrostatic balance.
    At least one temperature in ta (:math:`T`), theta (:math:`\\theta`),
    or thetal (:math:`\\theta_l`) should be given.
    If :math:`\\theta_l` is given, one simply assumes that :math:`\\theta = \\theta_l`. 
    If :math:`\\theta` is given, it is first converted to temperature :math:`T` 
    using :func:`.thermo.theta2t`.
    
    Basically, in case temperature :math:`T` is provided, at level :math:`k`, one uses:

    .. math:: \\ln\\frac{p_k}{p_{k-1}} = -\\int_{z_k-1}^{z_{k}} \\frac{g}{R T} dz
    .. math:: \\ln\\frac{p_k}{p_{k-1}} \\approx -0.5 g\\times
                                               \\left(\\frac{1}{R_k T_k} 
                                               + \\frac{1}{R_{k-1}T_{k-1}}\\right)
                                                 \\times (z_k - z_{k-1})
    
    In case potential temperature :math:`\\theta` is provided, at level :math:`k`, one uses:

    .. math:: To be completed

    In case both temperature and potential temperature are provided, potential temperature is
    use as a priority.

    If specific humidity :math:`q_v` is given, :math:`R = R_d + q_v(R_v-R_d)`,
    otherwise :math:`R = R_d`.
    
    Default values for :math:`g`, :math:`R_d`, :math:`R_v`, :math:`p_0` and :math:`\\kappa`
    are set in :mod:`.constants`.

    Parameters
    ----------
    z : array-like of size N
        Altitude :math:`z` (m).
    ps : :class:`float`
        Surface pressure (Pa).
    ta : array-like of size N, optional
        Temperature :math:`T` (K).
    theta : array-like of size N, optional
        Potential temperature :math:`\\theta` (K).
    thetal : array-like of size N, optional
        Liquid water potential temperature :math:`\\theta_l` (K).
    qv : array-like of size N, optional
        Specific humidity (kg kg-1). If given, the weight of water vapor is accounted for.
    g : :class:`float`, optional
        Gravity acceleration :math:`g` (m s-2). Value from :data:`.constants.g`
    Rd : :class:`float`, optional
        Dry air gas constant :math:`R_d` (J kg-1 K-1). Value from :data:`.constants.Rd`
    Rv : :class:`float`, optional
        Water vapor gas constant :math:`R_v` (J kg-1 K-1). Value from :data:`.constants.Rv`
    p0 : :class:`float`, optional
        Reference pressure :math:`p_0` (Pa). Value from :data:`.constants.p0`
    kappa : :class:`float`, optional
        Ratio :math:`\\kappa` of the gas constant to the gas heat capacity (no units).
        Default is for dry air taken from :data:`.constants.kappa`.

    Returns
    -------
    :class:`ndarray` of size N
        Air pressure :math:`p` (Pa).
    """
    if theta is None and ta is None and thetal is None:
        raise ValueError("theta, thetal or ta must be provided")

    if theta is None and thetal is not None:
        logger.warning("Using thetal as theta (approximation)")
        theta = thetal

    z = _asarray(z)
    nlev, = z.shape
    ps = float(ps)

    if qv is None:
        R = np.zeros_like(z) + Rd
    else:
        R = Rd + _asarray(qv) * (Rv - Rd)

    p = np.zeros_like(z)
    p[0] = ps

    if theta is not None: # use theta as a priority
        theta = _asarray(theta)
        integ = 0.
        for k in range(1,nlev):
            dz = z[k] - z[k - 1]
            integ = integ + (g/(R[k-1]*theta[k-1])+g/(R[k]*theta[k]))/2*dz
            tmp = ps**kappa-p0**kappa*kappa*integ
            try: p[k] = math.exp(math.log(tmp)/kappa)
            except:
                print(k, z[k-1], z[k], theta[k-1], theta[k], R[k-1], R[k], integ, tmp)
                print(math.log(tmp), kappa)
                raise ValueError('Error in z2p vertical integration under theta case')

    else: # use ta instead
        ta = _asarray(ta)
        for k in range(1, nlev):
            dz = z[k] - z[k - 1]
            p[k] = p[k - 1] * math.exp(
                -(g/(R[k]*ta[k]) + g/(R[k-1]*ta[k - 1]))*0.5 * dz
            )

    return p


def p2z(
    p,
    zs=0.0,
    theta=None, ta=None, thetal=None,
    qv=None,
    g=cc.g, Rd=cc.Rd, Rv=cc.Rv, p0=cc.p0, kappa=cc.kappa
):
    """
    Compute the altitude of a set of given pressure levels using the hydrostatic balance.
    At least one temperature in ta (:math:`T`), theta (:math:`\\theta`),
    or thetal (:math:`\\theta_l`) should be given.
    If :math:`\\theta_l` is given, one simply assumes that :math:`\\theta = \\theta_l`. 
    If :math:`\\theta` is given, it is first converted to temperature :math:`T` 
    using :func:`.thermo.theta2t`.
    
    Basically, at level :math:`k`, one uses:

    .. math:: z_k - z_{k-1} = -\\int_{p_k-1}^{p_{k}} \\frac{R T}{g} d\\ln p
    .. math:: z_k - z_{k-1}  \\approx -\\frac{0.5\\times(R_k T_k+R_{k-1}T_{k-1})}{g}
                                                 \\times\\ln\\frac{p_k}{p_{k-1}}
    
    If specific humidity :math:`q_v` is given, :math:`R = R_d + q_v(R_v-R_d)`,
    otherwise :math:`R = R_d`.
    Default values for :math:`g`, :math:`R_d`, :math:`R_v`, :math:`p_0` and :math:`\\kappa`
    are set in :mod:`.constants`.

    Parameters
    ----------
    p : array-like of size N
        Air pressure :math:`p` (m).
    zs : :class:`float`, optional
        Surface altitude (m). Default is 0.
    ta : array-like of size N, optional
        Temperature :math:`T` (K).
    theta : array-like of size N, optional
        Potential temperature :math:`\\theta` (K).
    thetal : array-like of size N, optional
        Liquid water potential temperature :math:`\\theta_l` (K).
    qv : array-like of size N, optional
        Specific humidity (kg kg-1). If given, the weight of water vapor is accounted for.
    g : :class:`float`, optional
        Gravity acceleration :math:`g` (m s-2). Value from :data:`.constants.g`
    Rd : :class:`float`, optional
        Dry air gas constant :math:`R_d` (J kg-1 K-1). Value from :data:`.constants.Rd`
    Rv : :class:`float`, optional
        Water vapor gas constant :math:`R_v` (J kg-1 K-1). Value from :data:`.constants.Rv`
    p0 : :class:`float`, optional
        Reference pressure :math:`p_0` (Pa). Value from :data:`.constants.p0`
    kappa : :class:`float`, optional
        Ratio :math:`\\kappa` of the gas constant to the gas heat capacity (no units).
        Default is for dry air taken from :data:`.constants.kappa`.

    Returns
    -------
    :class:`ndarray` of size N
        Altitude :math:`z` (m).
    """
    if theta is None and ta is None and thetal is None:
        raise ValueError("theta, thetal or ta must be provided")

    p = _asarray(p)

    if ta is None:
        if theta is None:
            logger.warning("Using thetal as theta (approximation)")
            theta = thetal
        ta = theta2t(p, theta, p0=p0, kappa=kappa)

    ta = _asarray(ta)

    if qv is None:
        R = np.zeros_like(p) + Rd
    else:
        R = Rd + _asarray(qv) * (Rv - Rd)

    z = np.zeros_like(p)
    z[0] = zs

    for k in range(1, len(p)):
        z[k] = z[k - 1] + (
            0.5 * (R[k]*ta[k] + R[k-1]*ta[k - 1]) / g
        ) * math.log(p[k - 1]/p[k])

    return z


# =============================================================================
# Interpolation utilities
# =============================================================================
def zlev2plev(zlev, z, p):
    """
    Linear interpolation of pressure at a given altitude
    The function :func:`numpy.interp` is used.

    Parameters
    ----------
    zlev : :class:`float`
        Altitude at which one wants to interpolate pressure (m).
    z : array-like of size N
        Altitudes at which pressure is known (m).
    p : array-like of size N
        Air pressure at levels **z** (Pa).

    Returns
    -------
    :class:`float`
        Air pressure at level **zlev** (Pa).
    """
    return np.interp(zlev, z, p)


def plev2zlev(plev, z, p):
    """
    Linear interpolation of altitude at a given pressure level.
    The function :func:`numpy.interp` is used.

    Parameters
    ----------
    plev : :class:`float`
        Pressure level at which one wants to interpolate altitude (Pa).
    z : array-like of size N
        Altitudes at which pressure is known (m).
    p : array-like of size N
        Air pressure at levels **z** (Pa).

    Returns
    -------
    :class:`float`
        Altitude of the pressure **plev** (m).
    """
    return np.interp(plev, p[::-1], z[::-1])


# =============================================================================
# Humidity from dew-point / Goff–Gratch
# =============================================================================
def rh2qv_GG(hur, pressure, temp):
    """
    Convert relative humidity to specific humidity
    using the Goff-Gratch formulation.
    Adapted from http://climatologie.u-bourgogne.fr/data/matlab/goff_gratch.m

    Parameters
    ----------
    hur : :class:`float` or array-like
        Relative humidity (%).
    pressure : :class:`float` or array-like
        Pressure (Pa).
    temp : :class:`float` or array-like
        Temperature (K).

    Returns
    -------
    :class:`float` or :class:`ndarray`
        Specific humidity.
    """
    hur = _asarray(hur)
    temp = _asarray(temp)
    pressure = _asarray(pressure)

    es_ice = 6.1071 * np.exp(
        (-9.09718 * (273.16 / temp - 1)
         - 3.5654 * np.log10(273.16 / temp)
         + 0.876793 * (1 - temp / 273.16)) * math.log(10)
    )

    es_water = 1013.246 * np.exp(
        (-7.90298 * (373.16 / temp - 1)
         + 5.02808 * np.log10(373.16 / temp)) * math.log(10)
    )

    es = np.where(temp < 273.15, es_ice, es_water)
    ws = 0.62197 * es / (pressure / 100.0 - 0.378 * es)

    return (hur / 100.0) * ws


def td2qv(td, pressure, units: str = "kg kg-1"):
    """
    Convert dew-point temperature to specific humidity.
    Requires the MetPy module.

    Parameters
    ----------
    td : :class:`float` or array-like
        Dew-point temperature (K)
    pressure : :class:`float` or array-like
        Pressure (Pa).
    units : :class:`str`, optional
        Units for the specific humidity.
        Should be 'kg kg-1' (default) or 'g kg-1'.

    Returns
    -------
    :class:`float` or :class:`ndarray`
        Specific humidity.
    """
    if not HAS_METPY:
        raise RuntimeError("MetPy is required for td2qv")

    _check_units(units, ("kg kg-1", "g kg-1"))

    td = _asarray(td)
    pressure = _asarray(pressure) * Munits.Pa
    td = td * Munits.kelvin

    qv = specific_humidity_from_dewpoint(pressure, td).magnitude

    if units == "kg kg-1":
        return qv

    return qv * 1000.0


if __name__ == "__main__":
    print("Running basic self-tests for thermo.py")

    # ------------------------------------------------------------------
    # Realistic atmospheric profile
    # ------------------------------------------------------------------
    z = np.arange(0.0, 12001.0, 500.0)  # m

    # Temperature profile
    T0 = 288.15       # K
    lapse_rate = -6.5e-3  # K/m
    temp = _init_temp(z)#, T0, lapse_rate)

    print("Temperature")
    print(temp)

    # Specific humidity profile
    qv0 = 0.010       # kg/kg
    Hq = 2000.0       # m
    qv = _init_qv(z)#, qv0, Hq=Hq)
    print("Specific humidity")
    print(qv)

    print("Water vapor mixing ratio")
    print(qt2rt(qv))

    # Surface pressure
    ps = 100000.0 # Pa


    pressure = z2p(z, ps, ta=temp, qv=qv)
    print("Pressure")
    print(pressure)

    z2 = p2z(pressure, ta=temp, qv=qv)
    print("Recomputed altitude")
    print(z2)

    theta = t2theta(pressure, temp)
    print("Potential temperature")
    print(theta)
    
    z2 = p2z(pressure, theta=theta, qv=qv)
    print("Recomputed altitude using theta")
    print(z2)

    p2 = z2p(z, ps, theta=theta, qv=qv)
    print("Recomputed pressure using theta")
    print(p2)

    print("Relative humidity")
    print(qt2hur(qv, pressure, temp))
