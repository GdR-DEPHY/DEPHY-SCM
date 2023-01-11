#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 11 December 2019

@author: Romain Roehrig

Modification
  2020/11/12, R. Roehrig: update for improved case definition interface.
  2020/11/16, C. Rio: AMMA case
"""

## AMMA/REF original case definition

import netCDF4 as nc
import numpy as np

from dephycf.Case import Case

################################################
# 0. General configuration of the present script
################################################

lplot = True # plot all the variables
lverbose = False # print information about variables and case

################################################
# 1. General information about the case
################################################


case = Case('AMMA/REF',
        lat=13.47,
        lon=2.18,
        startDate="20060710060000",
        endDate="20060711000000",
        surfaceType='land',
        zorog=0.)


case.set_title("Forcing and initial conditions for AMMA case - Original definition")
case.set_reference(" https://rmets.onlinelibrary.wiley.com/doi/full/10.1002/qj.903; Couvreux et al. (2012, QJRMS)")
case.set_author("C. Rio")
case.set_script("DEPHY-SCM/AMMA/REF/driver_DEF.py")
case.set_comment("Use of file amma.nc")

# time units are expected to be seconds since startDate
#t0 = 0 # 18:00 UTC, 15 July 2006
#t1 = 72*3600 # 72-hour long simulation

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
# Modify first surface level value - Seems more consistent
height[0] = 0

case.add_init_height(height,lev=height,levtype='altitude')

# Pressure
pressure  = fin['pp'][:]
case.add_init_pressure(pressure,lev=height,levtype='altitude')

# Zonal and meridional wind
u  = fin['u'][:]
v  = fin['v'][:]
case.add_init_wind(u=u,v=v,lev=height,levtype='altitude')

# potential temperature
temp = fin['temp'][:]
case.add_init_temp(temp,lev=height,levtype='altitude')

# Vapor water content
qv = fin['qv'][:]
case.add_init_qv(qv,lev=height,levtype='altitude')

################################################
# 3. Forcing
################################################


# Forcing time axis
timeForc = fin['time'][:]*86400-6.*3600.

# large-scale velocity
w = fin['dw'][:,:]
case.add_vertical_velocity(w=w,lev=height,levtype='altitude',time=timeForc)

#Moisture advection
hq = fin['dq'][:,:]
case.add_qv_advection(hq,time=timeForc,lev=height,levtype='altitude')

#temperature advection and radiation
hT = fin['dt'][:,:]
case.add_theta_advection(hT,include_rad=True,lev=height,levtype='altitude',time=timeForc)

# Surface fluxes
sens=fin['sens'][:]
flat=fin['flat'][:]
case.add_surface_fluxes(sens,flat,time=timeForc,forc_wind='z0',z0=0.001)

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
