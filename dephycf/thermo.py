#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Provide a set of useful thermodynamical functions

Functions
---------
rt2qt(rt, units='kg kg-1')
    Compute total water knowing total water mixing ratio
qt2rt(qt, units='kg kg-1')
    Compute total water mixing ratio knowing total water
hur2qt(hur, pressure, temp, units='kg kg-1')
    Compute total water knowing relative humidity
hur2rt(hur, pressure, temp, units='kg kg-1')
    Compute total water mixing ratio knowing relative humidity
qt2hur(qt, pressure, temp, units='kg kg-1')
    Compute relative humidity knowin total water
rt2hur(rt, pressure, temp, units='kg kg-1')
    Compute relative humidity knowing total water mixing ratio
advrt2advqt(rt=None, advrt=None, rt_units='kg kg-1')
    Compute total water advection knowing total water mixing ratio advection
advqt2advrt(qt=None, advqt=None, qt_units='kg kg-1')
    Compute total water mixing ratio advection knowing total water advection
theta2t(p=None, theta=None, p0=cc.p0, kappa=cc.kappa)
    Compute temperature knowing potential temperature
t2theta(p=None, temp=None, p0=cc.p0, kappa=cc.kappa)
    Compute potential temperature knowing temperature
z2p(thetal=None, theta=None, ta=None, z=None, ps=None, qv=None, g=cc.g, Rd=cc.Rd, Rv=cc.Rv, p0=cc.p0, kappa=cc.kappa)
    Compute the pressure of a set of given altitudes
p2z(thetal=None, theta=None, ta=None, p=None, zs=0., qv=None, g=cc.g, Rd=cc.Rd, Rv=cc.Rv, p0=cc.p0, kappa=cc.kappa)
    Compute the altitude of a given set of pressure levels
zlev2plev(zlev, z, p)
    Compute the pressure of a given altitude (interpolation)
plev2zlev(plev, z, p)
    Compute the altitude of a given pressure level (interpolation)
rh2qv(rh,temp,pres)
    Compute the specific humidity knowing the relative humidity

Created on 27 November 2019
@author: Romain Roehrig
"""

import logging
logger = logging.getLogger(__name__)

import math
import numpy as np

try:
    from metpy.calc import mixing_ratio_from_relative_humidity, relative_humidity_from_mixing_ratio
    from metpy.units import units as Munits
except ImportError:
    logger.debug('Cannot load metpy library')
    logger.debug('Using relative humidity is not yet possible without this library')
except:
    raise

from . import constants as cc

#############################
def rt2qt(rt, units='kg kg-1'):
    """Compute total water knowing total water mixing ratio

    Parameters
    ----------
    rt : float, array
        Total water mixing ratio
    units : str, optional
        Total water mixing ratio units (default is kg kg-1)

    Returns
    -------
    float, array
        Total water (in units)
    """

    if units == 'kg kg-1':
        return rt/(1+rt)
    elif units == 'g kg-1':
        return rt/(1.+rt/1000.)
    else:
        logger.error('units unknown: {0}'.forma(units))
        raise ValueError('units unknown: {0}'.forma(units))

#############################
def qt2rt(qt, units='kg kg-1'):
    """Compute total water mixing ratio knowing total water

    Parameters
    ----------
    qt : float, array
        Total water
    units : str, optional
        Total water units (default is kg kg-1)

    Returns
    -------
    float, array
        Total water mixing ratio (in units)
    """

    if units == 'kg kg-1':
        return qt/(1-qt)
    elif units == 'g kg-1':
        return qt/(1.-qt/1000.)
    else:
        logger.error('units unknown: {0}'.forma(units))
        raise ValueError('units unknown: {0}'.forma(units))

#############################
def hur2qt(hur, pressure, temp, units='kg kg-1'):
    """Compute total water knowing relative humidity

    Parameters
    ----------
    hur : float, array
        Relative humidity (no units)
    pressure : float, array
        Air pressure (Pa)
    temp : float, array
        Air temperature (K)
    units : str, optional
        Total water mixing ratio units (default is kg kg-1)

    Returns
    -------
    float, array
        Total water (in units)
    """

    p_loc = pressure * Munits.Pa
    temp_loc = temp * Munits.kelvin
    
    if isinstance(hur,float):
        rv = float(mixing_ratio_from_relative_humidity(p_loc, temp_loc, hur))
    else:
        nlev, = hur.shape
        rv = np.array([mixing_ratio_from_relative_humidity(p_loc[i], temp_loc[i], hur[i]) for i in range(nlev)])

    if units == 'kg kg-1':
        return rv/(1+rv)
    elif units == 'g kg-1':
        return rv/(1.+rv)*1000.
    else:
        logger.error('units unknown: {0}'.forma(units))
        raise ValueError('units unknown: {0}'.forma(units))

#############################
def hur2rt(hur, pressure, temp, units='kg kg-1'):
    """Compute total water mixing ratio knowing relative humidity

    Parameters
    ----------
    hur : float, array
        Relative humidity (no units)
    pressure : float, array
        Air pressure (Pa)
    temp : float, array
        Air temperature (K)
    units : str, optional
        Total water mixing ratio units (default is kg kg-1)

    Returns
    -------
    float, array
        Total water mixing ratio (in units)
    """

    p_loc = pressure * Munits.Pa
    temp_loc = temp * Munits.kelvin
    
    if isinstance(hur,float):
        rv = float(mixing_ratio_from_relative_humidity(p_loc, temp_loc, hur))
    else:
        nlev, = hur.shape
        rv = np.array([mixing_ratio_from_relative_humidity(p_loc[i], temp_loc[i], hur[i]) for i in range(nlev)])

    if units == 'kg kg-1':
        return rv
    elif units == 'g kg-1':
        return rv*1000.
    else:
        logger.error('units unknown: {0}'.forma(units))
        raise ValueError('units unknown: {0}'.forma(units))

#############################
def qt2hur(qt, pressure, temp, units='kg kg-1'):
    """Compute relative humidity knowing total water

    Parameters
    ----------
    qt : float, array
        Total water (in units)
    pressure : float, array
        Air pressure (Pa)
    temp : float, array
        Air temperature (K)
    units : str, optional
        Total water units (default is kg kg-1)

    Returns
    -------
    float, array
        Relative humidityr (no units)
    """

    p_loc = pressure * Munits.Pa
    temp_loc = temp * Munits.kelvin

    if units == 'kg kg-1':
        rv_loc = qt/(1.-qt)
    elif units == 'g kg-1':
        rv_loc = qt/1000./(1.+rv/1000.)
    else:
        logger.error('units unknown: {0}'.forma(units))
        raise ValueError('units unknown: {0}'.forma(units))

    if isinstance(qt,float):
        hur = float(relative_humidity_from_mixing_ratio(p_loc, temp_loc, rv_loc))
    else:
        nlev, = qt.shape
        hur = np.array([mixing_ratio_from_relative_humidity(p_loc[i], temp_loc[i], rv_loc[i]) for i in range(nlev)])

    return hur

#############################
def rt2hur(rt, pressure, temp, units='kg kg-1'):
    """Compute relative humidity knowing total water mixing ratio

    Parameters
    ----------
    rt : float, array
        Total water mixing ratio (in units)
    pressure : float, array
        Air pressure (Pa)
    temp : float, array
        Air temperature (K)
    units : str, optional
        Total water mixing ratio units (default is kg kg-1)

    Returns
    -------
    float, array
        Relative humidityr (no units)
    """

    p_loc = pressure * Munits.Pa
    temp_loc = temp * Munits.kelvin

    if units == 'kg kg-1':
        rv_loc = rt*1.
    elif units == 'g kg-1':
        rv_loc = rt/1000.
    else:
        logger.error('units unknown: {0}'.forma(units))
        raise ValueError('units unknown: {0}'.forma(units))

    if isinstance(qt,float):
        hur = float(relative_humidity_from_mixing_ratio(p_loc, temp_loc, rv_loc))
    else:
        nlev, =rt.shape
        hur = np.array([mixing_ratio_from_relative_humidity(p_loc[i], temp_loc[i], rv_loc[i]) for i in range(nlev)])

    return hur

#############################
def advrt2advqt(rt=None, advrt=None, rt_units='kg kg-1'):
    """Compute total water advection knowing total water mixing ratio advection

    Parameters
    ----------
    rt : float, array
        Total water mixing ratio
    advrt : float, array
        Total water advection
    rt_units : str, optional
        Units of mixing ratio (default is kg kg-1)

    Returns
    -------
    float, array
        Total water mixing ratio advection (units is rt_units per unit of time)
    """

    if rt is None:
        logger.error("rt is missing")
        raise ValueError("rt is missing")
    if advrt is None:
        logger.error("advrt is missing")
        raise ValueError("advrt is missing")

    if rt_units == 'kg kg-1':
        return advrt/((1+rt)*(1+rt))
    elif rt_units == 'g kg-1':
        return advrt/((1.+rt/1000.)*(1.+rt/1000.))
    else:
        logger.error('units unknown for rt: {0}'.forma(units))
        raise ValueError('units unknown for rt: {0}'.forma(units))

#############################
def advqt2advrt(qt=None, advqt=None, qt_units='kg kg-1'):
    """Compute total water mixing ratio advection knowing total water advection

    Parameters
    ----------
    qt : float, array
        Total water
    advqt : float, array
        Total water advection
    qt_units : str, optional
        Units of total water (default is kg kg-1)

    Returns
    -------
    float, array
        Total water mixing ratio advection (units is qt_units per unit of time)
    """

    if qt is None:
        logger.error("qt is missing")
        raise ValueError("qt is missing")
    if advqt is None:
        logger.error("advqt is missing")
        raise ValueError("advqt is missing")

    if qt_units == 'kg kg-1':
        return advqt/((1-qt)*(1-qt))
    elif qt_units == 'g kg-1':
        return advqt/((1.-qt/1000.)*(1.-qt/1000.))
    else:
        logger.error('units unknown for qt: {0}'.forma(units))
        raise ValueError('units unknown for qt: {0}'.forma(units))

#############################
def theta2t(p=None, theta=None, p0=cc.p0, kappa=cc.kappa):
    """Compute temperature knowing potential temperature

    Parameters
    ----------
    theta : float, array
        Potential temperature
    p : float, array
        Pressure
    p0 : float, optional
        Reference pressure
    kappa : float, optional
        Rd/Cp, used to compute theta

    Returns
    -------
    float, array
        Temperature
    """

    if theta is None:
        logger.error("theta is missing")
        raise ValueError("theta is missing")
    if p is None:
        logger.error("p is missing")
        raise ValueError("p is missing")

    temp = theta*(p/p0)**kappa
    return temp

#############################
def t2theta(p=None, temp=None, p0=cc.p0, kappa=cc.kappa):
    """Compute potential temperature knowing temperature

    Parameters
    ----------
    temp : float, array
        Temperature
    p : float, array
        Pressure
    p0 : float, optional
        Reference pressure
    kappa : float, optional
        Rd/Cp, used to compute theta

    Returns
    -------
    float, array
        Potential temperature
    """

    if temp is None:
        logger.error("temp is missing")
        raise ValueError("temp is missing")
    if p is None:
        logger.error("p is missing")
        raise ValueError("p is missing")

    theta = temp*(p0/p)**kappa
    return theta

#############################
def z2p(thetal=None, theta=None, ta=None,
        z=None, ps=None, qv=None,
        g=cc.g, Rd=cc.Rd, Rv=cc.Rv, p0=cc.p0, kappa=cc.kappa):
    """Compute the pressure of a set of given altitudes

    Parameters
    ----------
    thetal/theta/ta : array, optional
        Liquid water potential temperature/potential temperature/temperature
        At last one of these must be given
    z : array
        Altitude
    ps : float
        Surface pressure
    qv : array
        Specific humidity
    g : float, optional
        Gravity acceleration
    Rd : float, optional
        Gaz constant for dry air
    Rv : float, optional
        Gaz constant for water vapor
    p0 : float, optional
        Reference pressure
    kappa : float, optional
        Rd/Cp, used to compute theta

    Returns
    -------
    array
        Pressure of the given altitude levels
    """

    if (thetal is None) and (theta is None) and (ta is None):
        logger.error("thetal, theta and ta are missing. At least one of them should be given")
        raise ValueError("thetal, theta and ta are missing. At least one of them should be given")
    if z is None:
        logger.error("z is missing")
        raise ValueError("z is missing")
    if ps is None:
        logger.error("ps is missing")
        raise ValueError("ps is missing")

    if thetal is not None and theta is None and ta is None:
        theta = thetal

    if theta is not None:
        nlev, = theta.shape

        if qv is None:
            R = theta*0.+Rd
        else:
            R = Rd+qv*(Rv-Rd)

        integ = 0.

        p = np.zeros((nlev,),dtype=np.float64)
        p[0] = ps

        for ilev in range(1,nlev):
            dz = z[ilev]- z[ilev-1]    
            integ = integ + (g/(R[ilev-1]*theta[ilev-1])+g/(R[ilev]*theta[ilev]))/2*dz
            tmp = ps**kappa-p0**kappa*kappa*integ
            p[ilev] = math.exp(math.log(tmp)/kappa)
    else: # Use ta instead
        nlev, = ta.shape

        if qv is None:
            R = ta*0.+Rd
        else:
            R = Rd+qv*(Rv-Rd)

        integ = 0.

        p = np.zeros((nlev,),dtype=np.float64)
        p[0] = ps

        for ilev in range(1,nlev):
            dz = z[ilev]- z[ilev-1]    
            integ = integ + (g/(R[ilev-1]*ta[ilev-1])+g/(R[ilev]*ta[ilev]))/2*dz
            tmp = math.log(ps)-integ
            p[ilev] = math.exp(tmp)


    return p

#############################
def p2z(thetal=None, theta=None, ta=None,
        p=None, zs=0., qv=None,
        g=cc.g, Rd=cc.Rd, Rv=cc.Rv, p0=cc.p0, kappa=cc.kappa):
    """Compute the altitude of a given set of pressure levels

    Start from the surface

    Parameters
    ----------
    thetal/theta/ta : array, optional
        Liquid water potential temperature/potential temperature/temperature
        At last one of these must be given
    p : array
        Pressure
    zs : float, optional
        Surface altitude (default is 0)
    qv : array
        Specific humidity
    g : float, optional
        Gravity acceleration
    Rd : float, optional
        Gaz constant for dry air
    Rv : float, optional
        Gaz constant for water vapor
    p0 : float, optional
        Reference pressure
    kappa : float, optional
        Rd/Cp, used to compute theta

    Returns
    -------
    array
        Altitude of the given pressure levels
    """

    if (thetal is None) and (theta is None) and (ta is None):
        logger.error("thetal, theta and ta are missing. At least one of them should be given")
        raise ValueError("thetal, theta and ta are missing. At least one of them should be given")
    if p is None:
        logger.error("p is missing")
        raise ValueError("p is missing")

    if theta is not None and ta is None:
        ta = theta2t(theta=theta,p=p)
    if thetal is not None and ta is None:
        # we assume thetal=theta for ta calculation
        ta = theta2t(theta=thetal,p=p)

    nlev, = ta.shape

    if qv is None:
        R = ta*0.+Rd
    else:
        R = Rd+qv*(Rv-Rd)

    z = np.zeros((nlev,),dtype=np.float64)
    z[0] = zs

    for ilev in range(1,nlev):
        dz = (R[ilev-1]*ta[ilev-1]+R[ilev]*ta[ilev])/(2.0*g)*(math.log(p[ilev-1])-math.log(p[ilev]))
        z[ilev] = z[ilev-1] + dz

    return z

#############################
def zlev2plev(zlev, z, p):
    """Compute the pressure of a given altitude (interpolation)

    Parameters
    ----------
    zlev : float
        Altitude for which the pressure is sought
    z : array
        Altitude
    p : array
        Pressure

    Returns
    -------
    float
        Pressure of the given altitude
    """
    
    nlev, = z.shape
    plev = 110000
    if z[0] < z[1]:
        for ilev in range(0,nlev-1):
            if z[ilev] <= zlev and z[ilev+1] > zlev:
                plev = p[ilev] + (p[ilev+1]-p[ilev])/(z[ilev+1]-z[ilev])*(zlev-z[ilev])
    else:
        for ilev in range(0,nlev-1):
            if z[ilev+1] <= zlev and z[ilev] > zlev:
                plev = p[ilev] + (p[ilev+1]-p[ilev])/(z[ilev+1]-z[ilev])*(zlev-z[ilev])

    return plev

#############################
def plev2zlev(plev, z, p):
    """Compute the altitude of a given pressure level (interpolation)

    Parameters
    ----------
    plev : float
        Pressure level for which the altitude is sought
    z : array
        Altitude
    p : array
        Pressure

    Returns
    -------
    float
        Altitude of the given pressure level
    """

    nlev, = z.shape
    zlev = 0
    if p[0] > p[1]:
        for ilev in range(0,nlev-1):
            if p[ilev] >= plev and p[ilev+1] < plev:
                zlev = z[ilev] + (z[ilev+1]-z[ilev])/(p[ilev+1]-p[ilev])*(plev-p[ilev])
    else:
        for ilev in range(0,nlev-1):
            if p[ilev+1] >= plev and p[ilev] < plev:
                zlev = z[ilev] + (z[ilev+1]-z[ilev])/(p[ilev+1]-p[ilev])*(plev-p[ilev])

    return zlev

#############################
def rh2qv_GG(rh,temp,pres):
    """Compute the specific humidity knowing the relative humidity

    Based on Goff-Gratch equation
    From http://climatologie.u-bourgogne.fr/data/matlab/goff_gratch.m

    Parameters
    ----------
    rh : array
        Relative humidity in %, wrt ice if temp < 273.15, wrt liquid water if temp >= 273.15
    temp : array
        Temperature in K
    pres : array
        Pressure in Pa

    Returns
    -------
    float
        Specific humidity in kg kg-1
    """

    # Saturation water vapor pressure against ice
    eilog = -9.09718 * ((273.16/temp) -1.)
    eilog2 = -3.5654 * np.log10(273.16/temp)
    eilog3 = 0.876793 * (1. - (temp/273.16))
    es1=6.1071*np.exp((eilog+eilog2+eilog3)*math.log(10.))

    # against liquid water
    eilog=-7.90298*((373.16/temp) - 1.)
    eilog2=5.02808*np.log10((373.16/temp))
    eilog3=-1.3816e-7*(np.exp((11.344*(1.-(temp/373.16)))*math.log(10.)) -1.)
    eilog4=8.1328e-3*(np.exp((-3.49149*((373.16/temp) - 1.) )*math.log(10)) -1.)
    es2=1013.246*np.exp((eilog+eilog2+eilog3+eilog4)*math.log(10.))

    es = np.where(temp < 273.15, es1,es2)

    ws = 0.62197* es/(pres/100. - 0.378*es) # qsat

    # specific humidity in kg kg-1
    return (rh/100.0)*ws
