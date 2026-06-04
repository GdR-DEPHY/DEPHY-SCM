#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 16 November 2020

@author: Catherine Rio

Modification
  2024/01/26, F Couvreux: modification du cas car il a été défini en th/rv et comme tous les champs sont disponibles dans le fichier AMMA j'ai choisi de l'écrire en th & rv
  2026/05/07, N. Villefranque: clean for publication
"""

## AMMA/REF original case definition
## From https://rmets.onlinelibrary.wiley.com/doi/full/10.1002/qj.903

import netCDF4 as nc
import numpy as np

from datetime import datetime, timedelta
from dephycf.Case import Case

################################################
# 0. General configuration of the present script
################################################

lplot = True     # plot all the variables
lverbose = False # print information about variables and case

################################################
# 1. General information about the case
################################################


duration=18
tmin = datetime(2006, 7, 10, 6, 0)
tmax = tmin + timedelta(hours=duration)

case = Case('AMMA/REF',
        lat=13.47,
        lon=2.18,
        startDate=tmin,
        endDate=tmax,
        surfaceType='land',
        zorog=0.)


case.set_title("Forcing and initial conditions for AMMA case - Original definition")
case.set_reference("Couvreux et al. (2012, QJRMS)")
case.set_author("C. Rio")
case.set_script("DEPHY-SCM/AMMA/REF/driver_DEF.py")
case.set_comment("Use of file amma.nc")

################################################
# 2. Input netCDF file
################################################

fin = nc.Dataset('amma.nc','r')

################################################
# 2. Initial state
################################################

# Surface pressure
ps = 98800.
case.add_init_ps(ps)

# Height
height = fin['zz'][:]
case.add_init_height(height,lev=height,levtype='altitude')

# Pressure
pressure  = fin['pp'][:]
case.add_init_pressure(pressure,lev=height,levtype='altitude')

# Zonal and meridional wind
u  = fin['u'][:]
v  = fin['v'][:]
case.add_init_wind(u=u,v=v,lev=height,levtype='altitude')

# Potential temperature
theta = fin['theta'][:]
case.add_init_theta(theta,lev=height,levtype='altitude')

# Vapor water content
rv = fin['rv'][:]
case.add_init_rv(rv,lev=height,levtype='altitude')

################################################
# 3. Forcing
################################################

# Forcing time axis
# convert to hours from the beginning of simulation
timeForc = fin['time'][:]*86400-6.*3600.

# large-scale velocity
w = fin['dw'][:,:]
case.add_vertical_velocity(w=w,lev=height,levtype='altitude',time=timeForc)

#Moisture advection
hrv = fin['drv'][:,:]
case.add_rv_advection(hrv,time=timeForc,lev=height,levtype='altitude')

# Potential temperature advection and radiation
hth = fin['dth'][:,:]
case.add_theta_advection(hth,include_rad=True,lev=height,levtype='altitude',time=timeForc)

# Surface fluxes
sens=fin['sens'][:]
flat=fin['flat'][:]
case.add_surface_fluxes(sens,flat,time=timeForc,forc_wind='z0',z0=0.01)

################################################
# 4. Writing file
################################################

case.write('AMMA_REF_DEF_driver.nc')

if lverbose:
    case.info()

################################################
# 6. Ploting, if asked
################################################

if lplot:
    case.plot(rep_images='./images/driver_DEF/',timeunits='hours')
