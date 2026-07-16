#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 24 May 2023

@author: Fleur Couvreux

Modification
  2023/06/05, R. Roehrig, some fixes and cleaning
  2026/07/16, N. Villefranque: clean for publication.
"""

## LBA original case definition
# from Grabowski et al., 2006
# https://rmets.onlinelibrary.wiley.com/doi/full/10.1256/qj.04.147

import numpy as np

from datetime import datetime, timedelta
from dephycf.Case import Case

################################################
# 0. General configuration of the present script
################################################

lplot = True # plot all the variables
lverbose = False # print information about variables and case

################################################
# 1. General information about the case
################################################

duration=7
tmin = datetime(1999, 2, 23, 7, 30)
tmax = tmin + timedelta(hours=duration)

case = Case('LBA/OTHER_FLUX',
        lat=-8.,
        lon=-63.,
        startDate=tmin,
        endDate=tmax,
        surfaceType='land',
        zorog=0.)

case.set_title("Forcing and initial conditions for LBA case - Slightly modified fluxes")
case.set_reference("Grabowski et al. (QJRMS, 2006)")
case.set_author("F. Couvreux")
case.set_script("DEPHY-SCM/LBA/REF/driver_DEF.py")
case.set_modifications("Surface fluxes are read from a file instead of computed as in the REF paper.")

################################################
# 2. Initial state
################################################

# Surface pressure
ps = 99130.
case.add_init_ps(ps)

# Read all initial profiles from init.txt aux file
z, theta, rv, u, v = np.genfromtxt('../aux/init.txt', dtype=float).transpose()

# Zonal and meridional wind
case.add_init_wind(u=u, v=v, lev=z, levtype='altitude')

# Potential temperature
case.add_init_theta(theta, lev=z, levtype='altitude')

# Water vapor mixing ratio
case.add_init_rv(rv, lev=z, levtype='altitude')

################################################
# 3. Forcing
################################################

# Forcing time axis
timeForc = [0.,3600.,7200.,10800.,14400.,18000.,21600.]
nt = len(timeForc)

height_forc = np.genfromtxt('../aux/LBA_formatcommun_ZFR_1',dtype=float,usecols=0)
nlev, = height_forc.shape

data = {}
for var in ['u','v','dthdt']:
    data[var] = np.zeros((nt, nlev), dtype=float)

for it in range(0,nt):
    fin = f'../aux/LBA_formatcommun_ZFR_{it+1}'
    data['u'][it,:] = np.genfromtxt(fin,dtype=float,usecols=1)
    data['v'][it,:] = np.genfromtxt(fin,dtype=float,usecols=2)
    data['dthdt'][it,:] = np.genfromtxt(fin,dtype=float,usecols=6)

# Potential temperature advection, which includes radiation
case.add_theta_advection(data['dthdt'], lev=height_forc, levtype='altitude',
                         time=timeForc, include_rad=True)

# Wind nudging
case.add_wind_nudging(unudg=data['u'], vnudg=data['v'], timescale=3600.,
                      time=timeForc, lev=height_forc, levtype='altitude')

# Surface fluxes
timeflux, sens, flat = np.genfromtxt('../aux/flux_LBA',dtype=float).transpose()
case.add_surface_fluxes(sens, flat, time=timeflux, forc_wind='z0', z0=0.035)

################################################
# 4. Writing file
################################################

case.write('LBA_OTHER_FLUX_DEF_driver.nc')

if lverbose:
    case.info()

################################################
# 6. Ploting, if asked
################################################

if lplot:
    case.plot(rep_images='./images/driver_DEF/',timeunits='hours')
