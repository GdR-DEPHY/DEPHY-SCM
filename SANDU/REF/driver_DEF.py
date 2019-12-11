#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on 11 December 2019

@author: Romain Roehrig
"""

## SANDU/REF original case definition

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
t0 = 0 # 18:00 UTC, 15 July 2006
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

# Height
height = fin['height'][:]

# add surface level value
nlev, = height.shape
heightnew = np.zeros(nlev+1)
heightnew[1:] = height[:]

case.add_variable('height',heightnew,lev=heightnew,levtype='altitude')

# Zonal wind
u  = fin['u'][:]
unew = np.zeros(nlev+1)
unew[0] = -2.25
unew[1:] = u[:]

case.add_variable('u',unew,lev=heightnew,levtype='altitude')

# Meridional wind
v  = fin['v'][:]
vnew = np.zeros(nlev+1)
vnew[0] = -5.69
vnew[1:] = v[:]

case.add_variable('v',vnew,lev=heightnew,levtype='altitude')

# Liquid potential temperature
thetal = fin['thetal'][:]
thetalnew = np.zeros(nlev+1)
thetalnew[0] = 290.97
thetalnew[1:] = thetal[:]

case.add_variable('thetal',thetalnew,lev=heightnew,levtype='altitude')

# Total water content
qt = fin['qt'][:]
qtnew = np.zeros(nlev+1)
qtnew[0] = 10.46e-3
qtnew[1:] = qt[:]

case.add_variable('qt',qtnew,lev=heightnew,levtype='altitude')

################################################
# 3. Forcing
################################################

# Constant geostrophic wind
case.add_variable('ug',[unew,unew],time=[t0,t1],lev=heightnew,levtype='altitude')
case.add_variable('vg',[vnew,vnew],time=[t0,t1],lev=heightnew,levtype='altitude')

# Nudging of temperature and humidity at high levels
case.add_variable('thetal_nudging', [thetalnew,thetalnew],time=[t0,t1],lev=heightnew,levtype='altitude')
case.add_variable('qt_nudging',     [qtnew,qtnew],        time=[t0,t1],lev=heightnew,levtype='altitude')

# Forcing time axis
timeForc = fin['time'][:]

# Sea surface temperature
ts = fin['Tg'][:]
case.add_variable('ts',ts,time=timeForc)

# Constant large-scale velocity
w = fin['w'][:]
wnew = np.zeros(nlev+1)
wnew[0] = 0.
wnew[1:] = v[:]

case.add_variable('w',[wnew,wnew],time=[t0,t1],lev=heightnew,levtype='altitude')

# Ozone
o3 = fin['o3mmr'][:]
o3new = np.zeros(nlev+1)
o3new[0] = o3[0]
o3new[1:] = v[:]

case.add_variable('o3',[o3new,o3new],time=[t0,t1],lev=heightnew,levtype='altitude')

################################################
# 4. Attributes
################################################

# Radiation should be activated
case.set_attribute("trad",0)
# Geostrophic wind 
case.set_attribute("forc_geo",1)
# Temperature and humidity nudging above 3 km
case.set_attribute("nudging_thl",3600.)
case.set_attribute("nudging_qt",3600.)
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
