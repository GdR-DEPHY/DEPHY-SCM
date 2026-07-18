#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 30 October 2020

@author: Fleur Couvreux

Modification
  2021/01/06, R. Roehrig: update for improved case definition interface.
  2026/07/18, N. Villefranque: clean for publication.
"""

## IHOP original case definition
## a realistic convective boundary layer growth over Oklahoma
## From 10.1256/qj.04.167

import numpy as np

from datetime import datetime, timedelta
from dephycf.Case import Case

################################################
# 0. General configuration of the present script
################################################

lplot    = True  # plot all the variables
lverbose = False # print information about variables and case

################################################
# 1. General information about the case
################################################

duration = 7
tmin = datetime(2002,6,14,12,00)
tmax = tmin + timedelta(hours=duration)

case = Case('IHOP/REF',
        lat=36.56,
        lon=-100.61,
        startDate=tmin,
        endDate=tmax,
        surfaceType='land',
        zorog=0.)

case.set_title("Forcing and initial conditions for IHOP case - Original definition")
case.set_reference("Couvreux et al. (2005, QJRMS)")
case.set_author("F. Couvreux")
case.set_script("DEPHY-SCM/IHOP/REF/driver_DEF.py")

################################################
# 2. Initial state
################################################

# Surface pressure
ps = 91800.
case.add_init_ps(ps)

# Initial profiles : read from aux file
z, theta, rv, u, v = np.genfromtxt('../aux/init.txt', dtype=float).transpose()

# Zonal and meridional wind
case.add_init_wind(u=u, v=v, lev=z, levtype='altitude')

# Potential temperature
case.add_init_theta(theta, lev=z, levtype='altitude')

# Water vapor mixing ratio
case.add_init_rv(rv/1000., lev=z, levtype='altitude')

################################################
# 3. Forcing
################################################

# Atmospheric forcing (12h, 15h, 18h UTC)
timeF= [0.,10800.,21600.]
ntf = len(timeF)

zforc, ug, vg, w, theta_adv, rv_adv = np.genfromtxt('../aux/atm_forcing.txt',
                                      dtype=float).transpose().reshape((6,ntf,-1))
zforc = zforc[0]

# Geostrophic wind
case.add_geostrophic_wind(ug=ug,vg=vg,time=timeF,lev=zforc,levtype='altitude')

# Vertical velocity
case.add_vertical_velocity(w=w,time=timeF,lev=zforc,levtype='altitude')

# Advection of potential temperature
case.add_theta_advection(theta_adv*1e-3, time=timeF, lev=zforc,
                         levtype='altitude', include_rad=True)

# Advection of water vapor mixing ratio
case.add_rv_advection(rv_adv*1e-6, time=timeF, lev=zforc, levtype='altitude')

# Surface Forcing
timeSfc, shf, lhf = np.genfromtxt('../aux/surface_flux_forcings.txt',
                                  dtype=float).transpose()

case.add_surface_fluxes(sens=shf, lat=lhf, time=timeSfc, 
                        forc_wind='z0', z0=0.1)

################################################
# 4. Writing file
################################################

case.write('IHOP_REF_DEF_driver.nc')

if lverbose:
    case.info()

################################################
# 5. Ploting, if asked
################################################

if lplot:
    case.plot(rep_images='./images/driver_DEF/',timeunits='hours')
