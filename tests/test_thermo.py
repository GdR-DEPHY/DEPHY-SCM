import numpy as np
import pytest

from dephycf import thermo
from dephycf import constants as cc


# =============================================================================
# Helpers
# =============================================================================
RTOL = 1e-6


# =============================================================================
# rt / qt conversions
# =============================================================================
def test_rt_qt_inverse_scalar():
    rt = 0.01
    qt = thermo.rt2qt(rt)
    rt2 = thermo.qt2rt(qt)

    assert np.isclose(rt, rt2, rtol=RTOL)


def test_rt_qt_inverse_array():
    rt = np.array([0.0, 0.01, 0.02])
    qt = thermo.rt2qt(rt)
    rt2 = thermo.qt2rt(qt)

    assert np.allclose(rt, rt2, rtol=RTOL)


def test_rt_qt_units_gkg():
    rt = np.array([10.0, 20.0])  # g kg-1
    qt = thermo.rt2qt(rt, units="g kg-1")
    rt2 = thermo.qt2rt(qt, units="g kg-1")

    assert np.allclose(rt, rt2, rtol=RTOL)


def test_rt_qt_invalid_units():
    with pytest.raises(ValueError):
        thermo.rt2qt(0.01, units="foo")


# =============================================================================
# Advection conversions
# =============================================================================
def test_adv_rt_qt_inverse():
    rt = np.array([0.01, 0.02])
    advrt = np.array([1e-7, 2e-7])

    advqt = thermo.advrt2advqt(rt, advrt)
    advrt2 = thermo.advqt2advrt(thermo.rt2qt(rt), advqt)

    assert np.allclose(advrt, advrt2, rtol=RTOL)


# =============================================================================
# Temperature / potential temperature
# =============================================================================
def test_theta_t_inverse():
    p = np.array([100000.0, 90000.0, 80000.0])
    theta = np.array([300.0, 305.0, 310.0])

    t = thermo.theta2t(p, theta)
    theta2 = thermo.t2theta(p, t)

    assert np.allclose(theta, theta2, rtol=RTOL)


def test_theta2t_scalar():
    t = thermo.theta2t(100000.0, 300.0)
    assert np.isclose(t, 300.0, rtol=RTOL)


# =============================================================================
# Hydrostatic balance
# =============================================================================
def test_z2p_p2z_consistency_isothermal():
    z = np.linspace(0.0, 1000.0, 11)
    ta = np.full_like(z, 300.0)
    ps = 100000.0

    p = thermo.z2p(z, ps, ta=ta)
    z2 = thermo.p2z(p, zs=0.0, ta=ta)

    assert np.allclose(z, z2, rtol=RTOL)

def test_z2p_p2z_consistency_isotheta():
    z = np.linspace(0.0, 1000.0, 11)
    theta = np.full_like(z, 300.0)
    ps = 100000.0

    p = thermo.z2p(z, ps, theta=theta)
    z2 = thermo.p2z(p, zs=0.0, theta=theta)

    assert np.allclose(z, z2, rtol=RTOL)

def test_z2p_requires_temperature():
    z = np.array([0.0, 100.0])
    with pytest.raises(ValueError):
        thermo.z2p(z, 100000.0)


def test_p2z_requires_temperature():
    p = np.array([100000.0, 90000.0])
    with pytest.raises(ValueError):
        thermo.p2z(p)


# =============================================================================
# Interpolation utilities
# =============================================================================
def test_zlev2plev_linear():
    z = np.array([0.0, 100.0, 200.0])
    p = np.array([100000.0, 90000.0, 80000.0])

    p50 = thermo.zlev2plev(50.0, z, p)
    assert np.isclose(p50, 95000.0, rtol=RTOL)


def test_plev2zlev_linear():
    z = np.array([0.0, 100.0, 200.0])
    p = np.array([100000.0, 90000.0, 80000.0])

    z90 = thermo.plev2zlev(90000.0, z, p)
    assert np.isclose(z90, 100.0, rtol=RTOL)


# =============================================================================
# MetPy-dependent tests
# =============================================================================
@pytest.mark.skipif(not thermo.HAS_METPY, reason="MetPy not available")
def test_hur_rt_qt_roundtrip():
    hur = 0.5
    p = 100000.0
    t = 300.0

    rt = thermo.hur2rt(hur, p, t)
    hur2 = thermo.rt2hur(rt, p, t)

    assert np.isclose(hur, hur2, rtol=1e-3)


@pytest.mark.skipif(not thermo.HAS_METPY, reason="MetPy not available")
def test_td2qv_units():
    td = 290.0
    p = 100000.0

    qv = thermo.td2qv(td, p)
    qv_g = thermo.td2qv(td, p, units="g kg-1")

    assert np.isclose(qv * 1000.0, qv_g, rtol=RTOL)

