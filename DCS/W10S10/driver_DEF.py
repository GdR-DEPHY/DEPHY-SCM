#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 30 May 2024

@author: Gaston Bidoux, Romain Roehrig, Catherine Rio

Modification
"""

## LBA/REF original case definition

import netCDF4 as nc
import math
import numpy as np

from dephycf.Case import Case
import dephycf.constants as CC

################################################
# 0. General configuration of the present script
################################################

lplot = True # plot all the variables
lverbose = False # print information about variables and case

################################################
# 1. General information about the case
################################################


case = Case('DCS/W10S10',
        lat=45.,
        lon=0.,
        startDate="2007-01-31 10:00:00",
        endDate="2007-01-31 16:00:00",
        surfaceType='land',
        zorog=0.)


case.set_title("Forcing and initial conditions for Deep Convection and Shear case - Initial idealized profile from Weisman and Klemp (1982) ; Wind maximum = 10 m s-1 and shear = 10 m s-1 km-1 (from surface)")
case.set_reference("Weisman and Klemp (1982, MWR)")
case.set_author("G. Bidoux, R. Roehrig, C. Rio")
case.set_script("DEPHY-SCM/DCS/W10S10/driver_DEF.py")

################################################
# 2. Initial state
################################################

# Surface pressure
ps = 100000.
case.add_init_ps(ps)

# Wind initial profile

hshear = 1000
r1 = math.sqrt(10)
r2 = 0.0003

def wind(z, u0=0, v0=0):
    nz, = z.shape
    u, v = np.zeros(nz,dtype=np.float32), np.zeros(nz,dtype=np.float32)
    u = np.where(z <= hshear, u0 - r1*np.cos(np.pi*z/(2*hshear)) + r1, u0 + r2*(z-hshear) + r1)
    v = np.where(z <= hshear, v0 + r1*np.sin(np.pi*z/(2*hshear)), v0 + r1)
    """
    if z <= hshear:
        u = u0 - r1*np.cos(np.pi*z/(2*hshear)) + r1
        v = v0 + r1*np.sin(np.pi*z/(2*hshear))
    else:
        u = u0 + r2*(z-zshear) + r1
        v = v0 + r1
    """
    return u, v

zwind = np.arange(0,20001,5)
u, v = wind(zwind)

case.add_init_wind(u=u, v=v, lev=zwind, levtype='altitude')

ztr = 12000.
theta0 = 300.
thetatr = 343
Ttr = 213
def theta(z):
    if z <= ztr:
        tmp = theta0 + (thetatr-theta0)*(z/ztr)**(5./4.)
    else:
        tmp = thetatr * math.exp(CC.g/(CC.Cpd*Ttr)*(z-ztr))

    return tmp

def hur(z):
    if z <= ztr:
       tmp = 1 - 3./4.*(z/ztr)**(5./4.)
    else:
       tmp = 0.25

    return tmp

zthermo = np.arange(0,20001,5)
theta_init = [theta(zloc) for zloc in zthermo]
hur_init = [hur(zloc) for zloc in zthermo]

# Temperature and relative humidity
case.add_init_theta(theta_init, lev=zthermo, levtype='altitude')
case.add_init_hur(hur_init, lev=zthermo, levtype='altitude')

################################################
# 3. Forcing
################################################

# Forcing time axis
timeForc = [0.,36000.]
shf = [200., 200.]
lhf = [350., 350.]
ts = [theta0, theta0]

case.add_surface_fluxes(shf,lhf,time=timeForc,forc_wind='z0',z0=0.035)
case.add_surface_temp(ts, time=timeForc)

################################################
# 4. Writing file
################################################

case.write('DCS_W10S10_DEF_driver.nc')

if lverbose:
    case.info()

################################################
# 6. Ploting, if asked
################################################

if lplot:
    case.plot(rep_images='./images/driver_DEF/',timeunits='hours')
