#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on 11 December 2019

@author: Romain Roehrig
"""

## ARM-Cumulus original case definition
## From http://projects.knmi.nl/eurocs/ARM/case_ARM_html/

import os
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

case = Case('SANDU/REF',
        lat=25,
        lon=-125,
        startDate="20060715180000",
        endDate="20060718180000",
        zorog=0.)

case.set_title("Forcing and initial conditions for SANDU composite reference transition case - Original definition")
case.set_reference("Sandu and Stevens (2011, JAS); Input file composite_ref_init.nc downloaded at http://gop.meteo.uni-koeln.de/~neggers/transitions/")
case.set_author("R. Roehrig")
case.set_script("DEPHY-SCM/SANDU/REF/driver_DEF.py")
case.set_comment("Use of file composite_ref_init.nc downloaded at http://gop.meteo.uni-koeln.de/~neggers/transitions/")

# time units are expected to be seconds since startDate
t0 = 0 # 00:00 UTC, 16 December 2004
t1 = 72*3600 # 72-hour long simulation

################################################
# 2. Input netCDF file
################################################

fin = nc.Dataset('composite_ref_init.nc','r')

################################################
# 2. Initial state
################################################

# Surface pressure
ps = fin['Ps'][0]
case.add_variable('ps',[ps,])

# Pressure profile
pressure = fin['lev'][:]
case.add_variable('pressure',pressure,lev=pressure,levtype='pressure')

# Height
height = fin['height'][:]
case.add_variable('height',height,lev=pressure,levtype='pressure')

# Zonal wind
u  = fin['u'][:]
case.add_variable('u',u,lev=pressure,levtype='pressure')

# Meridional wind
v  = fin['v'][:]
case.add_variable('v',v,lev=pressure,levtype='pressure')

# Temperature
temp = fin['T'][:]
case.add_variable('temp',temp,lev=pressure,levtype='pressure')

# Liquid potential temperature
thetal = fin['thetal'][:]
case.add_variable('thetal',thetal,lev=pressure,levtype='pressure')

# Total water content
qt = fin['qt'][:]
case.add_variable('qt',qt,lev=pressure,levtype='pressure')

################################################
# 3. Forcing
################################################

# Constant wind nudging
case.add_variable('u_nudging',[u,u],time=[t0,t1],lev=pressure,levtype='pressure')
case.add_variable('v_nudging',[v,v],time=[t0,t1],lev=pressure,levtype='pressure')

# Nudging of temperature and humidity at high levels
case.add_variable('temp_nudging',   [temp,temp],    time=[t0,t1],lev=pressure,levtype='pressure')
case.add_variable('thetal_nudging', [thetal,thetal],time=[t0,t1],lev=pressure,levtype='pressure')
case.add_variable('qt_nudging',     [qt,qt],        time=[t0,t1],lev=pressure,levtype='pressure')

# Forcing time axis
timeForc = fin['time'][:]

# Sea surface temperature
ts = fin['Tg'][:]
case.add_variable('ts',ts,time=timeForc)

# Constant large-scale velocity
w = fin['w'][:]
case.add_variable('w',[w,w],time=[t0,t1],lev=pressure,levtype='pressure')

# Ozone
o3 = fin['o3mmr'][:]
case.add_variable('o3',[o3,o3],time=[t0,t1],lev=pressure,levtype='pressure')

################################################
# 4. Attributes
################################################

# Radiation should be activated
case.set_attribute("trad",0)
# Wind nudging from the surface, nudging time = 1 hour
case.set_attribute("nudging_u",3600.)
case.set_attribute("nudging_v",3600.)
case.set_attribute("z_nudging_u",0)
case.set_attribute("z_nudging_v",0)
# Temperature and humidity nudging above 3 km
case.set_attribute("nudging_t",3600.)
case.set_attribute("nudging_thl",3600.)
case.set_attribute("nudging_qt",3600.)
case.set_attribute("z_nudging_t",3000)
case.set_attribute("z_nudging_thl",3000)
case.set_attribute("z_nudging_qt",3000)
# Vertical velocity
case.set_attribute("forc_w",1)
# Surface flux forcing, wind stress is computed using z0
case.set_attribute("surfaceType","ocean")
case.set_attribute("surfaceForcing","ts")

################################################
# 5. Writing file
################################################

case.write('SANDU_REF_DEF_driver.nc',verbose=False)

if lverbose:
    case.info()

################################################
# 6. Ploting, if asked
################################################

if lplot:
    case.plot(rep_images='./images/driver_DEF/',timeunits='hours')
