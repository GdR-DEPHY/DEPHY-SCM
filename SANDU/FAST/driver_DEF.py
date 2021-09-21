#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on 09 January 2021

@author: Romain Roehrig

"""

## SANDU/FAST original case definition

import sys
sys.path = ['../../utils/',] + sys.path

import netCDF4 as nc
import numpy as np

from Case import Case

################################################
# 0. General configuration of the present script
################################################

lplot = True # plot all the variables
lverbose = False # print information about variables and case

################################################
# 1. General information about the case
################################################

case = Case('SANDU/FAST',
        lat=25,
        lon=-125,
        startDate="20060715180000",
        endDate="20060718180000",
        surfaceType='ocean',
        zorog=0.)

case.set_title("Forcing and initial conditions for SANDU composite fast transition case - Original definition")
case.set_reference("Sandu and Stevens (2011, JAS); Input file composite_fast_init.nc downloaded at http://gop.meteo.uni-koeln.de/~neggers/transitions/")
case.set_author("R. Roehrig")
case.set_script("DEPHY-SCM/SANDU/FAST/driver_DEF.py")
case.set_comment("Use of file composite_fast_init.nc downloaded at http://gop.meteo.uni-koeln.de/~neggers/transitions/")

################################################
# 2. Input netCDF file
################################################

fin = nc.Dataset('composite_fast_init.nc','r')

################################################
# 2. Initial state
################################################

# Surface pressure
ps = fin['Ps'][0]
case.add_init_ps(ps)

# Height
height = fin['height'][:]
# Modify first surface level value - Seems more consistent
height[0] = 0

case.add_init_height(height,lev=height,levtype='altitude')

# Pressure
pressure  = fin['lev'][:]
case.add_init_pressure(pressure,lev=height,levtype='altitude')

# Zonal and meridional wind
u  = fin['u'][:]
v  = fin['v'][:]
case.add_init_wind(u=u,v=v,lev=height,levtype='altitude')

# Liquid potential temperature
thetal = fin['thetal'][:]
case.add_init_thetal(thetal,lev=height,levtype='altitude')

# Temperature
temp = fin['T'][:]
case.add_init_temp(temp,lev=height,levtype='altitude')

# Total water content
qt = fin['qt'][:]
case.add_init_qt(qt,lev=height,levtype='altitude')

################################################
# 3. Forcing
################################################

# Constant geostrophic wind
case.add_geostrophic_wind(ug=u,vg=v,lev=height,levtype='altitude')

# Nudging of temperature and humidity at high levels
case.add_temp_nudging(temp,timescale=3600.,z_nudging=3000,lev=height,levtype='altitude')
case.add_thetal_nudging(thetal,timescale=3600.,z_nudging=3000,lev=height,levtype='altitude')
case.add_qt_nudging(qt,timescale=3600.,z_nudging=3000,lev=height,levtype='altitude')

# Constant large-scale velocity
w = fin['w'][:]
case.add_vertical_velocity(w=w,lev=height,levtype='altitude')

# Ozone
ro3 = fin['o3mmr'][:]
case.add_ozone(ro3=ro3,lev=height,levtype='altitude')

# Forcing time axis
timeForc = fin['time'][:]

# Sea surface temperature
ts = fin['Tg'][:]
case.add_forcing_ts(ts,time=timeForc)

################################################
# 4. Writing file
################################################

case.write('SANDU_FAST_DEF_driver.nc')

if lverbose:
    case.info()

################################################
# 6. Ploting, if asked
################################################

if lplot:
    case.plot(rep_images='./images/driver_DEF/',timeunits='hours')
