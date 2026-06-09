#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 09 June 2026

@author: DEPHY team

Modification
"""

## Cirrus case (idealized) - Borella 2025

import netCDF4 as nc
import numpy as np

from datetime import datetime, timedelta
from dephycf.Case import Case
from dephycf import constants

################################################
# 0. General configuration of the present script
################################################

lplot    = True  # plot all the variables
lverbose = False # print information about variables and case

################################################
# 1. General information about the case
################################################

duration=12
tmin = datetime(2022, 12, 27, 00, 00)
tmax = tmin + timedelta(hours=duration)

case = Case('CIRRUS/REF',
        lat=48.4, # BREST 
        lon=-4.5,
        startDate=tmin,
        endDate=tmax,
        surfaceType='land',
        zorog=0)

case.set_title("Forcing and initial conditions for CIRRUS")
case.set_reference("Borella et al. 2025")
case.set_author("DEPHY Team")
case.set_script("DEPHY-SCM/CIRRUS/REF/driver_DEF.py")

################################################
# 2. Initial state
################################################

# Surface pressure
ps = 101600.
case.add_init_ps(ps)

htop = 12000

# Zonal and meridional wind
zu = zv = [0, 500, 2500, 3000, 5000, 8000, htop]
u  = [0, 10, 10, 14, 14, 30, 30]
v  = [0]*len(u)

case.add_init_wind(u=u,ulev=zu,v=v,vlev=zv,levtype='altitude')

# Potential temperature
ztheta = [0,   250,   1100, 2500, 4300, 6500, 10500, htop]
theta  = [280, 283.5, 283.5, 289,  301,  311,  320,   335 ]
case.add_init_theta(theta,lev=ztheta,levtype='altitude')

# Total humidity
zqv = [0, 2000, 3000, 3800, 4200, 8000, htop]
qv =  [5.2,  2, 0.15, 0.15, 2, 0.1 , 0.1]
case.add_init_qv(np.array(qv)/1000.,lev=zqv,levtype='altitude') 

################################################
# 3. Forcing
################################################

# vertical velocities
wasc = 0.5 # m/s
zw = [0., 7700, 7800., 8500, 8600, htop]
w  = [0., 0,    wasc,  wasc, 0,       0]
case.add_vertical_velocity(w=w,lev=zw,levtype='altitude')

# add zero tendencies 
case.add_theta_advection([0]*len(zw), lev=zw, levtype="altitude", include_rad=False)

# Surface Forcing
case.add_surface_fluxes(sens=0,lat=0,forc_wind='z0',z0=0.01)

tskin = 280
case.add_surface_skin_temp(tskin)

################################################
# 4. Writing file
################################################

case.write('CIRRUS_REF_DEF_driver.nc')

if lverbose:
    case.info()

################################################
# 5. Ploting, if asked
################################################

if lplot:
    case.plot(rep_images='./images/driver_DEF/',timeunits='hours')
