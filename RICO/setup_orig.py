#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on 27 November 2019

@author: Romain Roehrig
"""

## RICO composite original case definition
## From http://projects.knmi.nl/rico/setup1d_composite.html

import os
import sys
sys.path = ['../utils/',] + sys.path

import time

import netCDF4 as nc
import numpy as np

import thermo
import SCM_utils as utils

data = {}

###############################
# 1. General information about the case
###############################


lat = 18
lon = -61.5

startDate = "20041216000000"
endDate =   "20041219000000"

tunits = 'seconds since 2004-12-16 0:0:0.0'
t0 = 0 # 00:00 UTC, 16 December 2004
t1 = t0 + 72*3600 # 72 hours

ps = 101540. # Surface pressure (Pa)
zorog = 0 # Altitude above sea level (m)

ts = 299.8 # Sea surface temperature (K)

# Shared axes
t0Axis = utils.Axis('t0',[t0,],name='Initial time',units=tunits,calendar='gregorian')
latAxis = utils.Axis('lat',[lat,],name='Latitude',units='degrees_north')
lonAxis = utils.Axis('lon',[lon,],name='Longitude',units='degrees_east')

levAxes = {}
data = {}

###############################
# 2. Initial state
###############################

data['ps'] = utils.Variable('ps',name='Surface Pressure',units='Pa',data=[ps,],time=t0Axis,lat=latAxis,lon=lonAxis)

# Zonal Wind
zu = range(0,30000+1,100)
nu = len(zu)
u = np.zeros(nu,dtype=np.float64)

for iz in range(0,nu):
    if zu[iz] < 4000:
        u[iz] = -9.9 + (-1.9+9.9)*zu[iz]/4000.
    elif zu[iz] < 12000:
        u[iz] = -1.9+(30+1.9)/(12000.-4000.)*(zu[iz]-4000.)
    elif zu[iz] < 13000:
        u[iz] = 30.
    elif zu[iz] < 20000:
        u[iz] = 30. - 30./(20000.-13000.)*(zu[iz]-13000.)
    else:
        u[iz] = 0.

levAxes['u'] = utils.Axis('lev_u',zu,name='Altitude for variable u',units='m')
data['u'] = utils.Variable('u',name='Zonal wind',units='m s-1',data=u,level=levAxes['u'],time=t0Axis,lat=latAxis,lon=lonAxis)

# Meridional Wind
zv = range(0,30000+1,30000)
nv = len(zv)
v = np.zeros(nv,dtype=np.float64) - 3.8

levAxes['v'] = utils.Axis('lev_v',zv,name='Altitude for variable v',units='m')
data['v'] = utils.Variable('v',name='Meridional wind',units='m s-1',data=v,level=levAxes['v'],time=t0Axis,lat=latAxis,lon=lonAxis)

# Temperature
zt = range(0,60000+1,10)
nt = len(zt)
temp = np.zeros(nt,dtype=np.float64)

for iz in range(0,nt):
    if zt[iz] < 740:
        temp[iz] = 299.2 + (292.0-299.2)/740.*zt[iz]
    elif zt[iz] < 4000:
        temp[iz] = 292.0 + (278.0-292.0)/(4000.-740.)*(zt[iz]-740.)
    elif zt[iz] < 15000:
        temp[iz] = 278.0 + (203.0-278.0)/(15000.-4000.)*(zt[iz]-4000.)
    elif zt[iz] < 17500:
        temp[iz] = 203.0 + (194.0 - 203.0)/(17500.-15000.)*(zt[iz]-15000.)
    elif zt[iz] < 20000:
        temp[iz] = 194.0 + (206.0 - 194.0)/(20000.-17500.)*(zt[iz]-17500.)
    else:
        temp[iz] = 206.0 + (270.0 - 206.0)/(60000.-20000.)*(zt[iz]-20000.)

levAxes['temp'] = utils.Axis('lev_temp',zt,name='Altitude for variable temp',units='m')
data['temp'] = utils.Variable('temp',name='Temperature',units='K',data=temp,level=levAxes['temp'],time=t0Axis,lat=latAxis,lon=lonAxis)

# Specific humidity
zqv = range(0,10000+1,10)
nqv = len(zqv)
qv = np.zeros(nqv,dtype=np.float64)

for iz in range(0,nqv):
    if zqv[iz] < 740:
        qv[iz] = 16.0 + (13.8-16.0)/740.*zqv[iz]
    elif zqv[iz] < 3260:
        qv[iz] = 13.8 + (2.4-13.8)/(3260.-740.)*(zqv[iz]-740.)
    elif zqv[iz] < 4000:
        qv[iz] = 2.4 + (1.8-2.4)/(4000.-3260.)*(zqv[iz]-3260.)
    elif zqv[iz] < 9000:
        qv[iz] = 1.8 + (0.-1.8)/(10000.-4000.)*(zqv[iz]-4000.) # 9000 instead of 10000 ??
    else:
        qv[iz] = 0.0

levAxes['qv'] = utils.Axis('lev_qv',zqv,name='Altitude for variable qv',units='m')
data['qv'] = utils.Variable('qv',name='Specific humidity',units='K',data=qv/1000.,level=levAxes['qv'],time=t0Axis,lat=latAxis,lon=lonAxis)


###############################
# 3. Forcing
###############################

timeAxes = {}

# Constant Geostrophic wind across the simulation
ug = np.zeros((2,nu),dtype=np.float64)
ug[0,:] = u[:]
ug[1,:] = u[:]

vg = np.zeros((2,nv),dtype=np.float64)
vg[0,:] = v[:]
vg[1,:] = v[:]

timeAxes = {}
for var in ['ug','vg']:
    timeAxes[var] = utils.Axis('time_{0}'.format(var),[t0,t1],name='Forcing time for {0}'.format(var),units=tunits)
    
levAxes['ug'] = utils.Axis('lev_ug',zu,name='Altitude for variable ug',units='m')
levAxes['vg'] = utils.Axis('lev_vg',zv,name='Altitude for variable vg',units='m')

data['ug']     = utils.Variable('ug',    name='Geostrophic Zonal Wind',     units='m s-1', data=ug,level=levAxes['ug'],time=timeAxes['ug'],lat=latAxis,lon=lonAxis)
data['vg']     = utils.Variable('vg',    name='Geostrophic Meridional Wind',units='m s-1', data=vg,level=levAxes['vg'],time=timeAxes['vg'],lat=latAxis,lon=lonAxis)

# Surface Forcing
timeAxes['ts'] = utils.Axis('time_ts',[t0,t1],name='Forcing time for ts',units=tunits)
data['ts'] = utils.Variable('ts',name='Sea surface temperature',units='K',data=[ts,ts],time=timeAxes['ts'],lat=latAxis,lon=lonAxis)

# Large-scale vertical velocity
zw = range(0,5000+1,10)
nw = len(zw)
w = np.zeros((2,nw),dtype=np.float64)

for iz in range(0,nw):
    if zw[iz] < 2260:
        w[:,iz] = -0.005/2260.*zw[iz]
    elif zw[iz] < 4000.:
        w[:,iz] = -0.005
    elif zw[iz] < 5000:
        w[:,iz] = -0.005 + 0.005/(5000.-4000.)*(zw[iz]-4000.)
    else:
        w[:,iz] = 0.0

timeAxes['w'] = utils.Axis('time_w',[t0,t1],name='Forcing time for w',units=tunits)
levAxes['w'] = utils.Axis('lev_w',zw,name='Altitude for variable w',units='m')
data['w'] = utils.Variable('w',name='Vertical velocity',units='m s-1',data=w,level=levAxes['w'],time=timeAxes['w'],lat=latAxis,lon=lonAxis)

# Large-scale advection of temperature + radiative tendency
ztadv = range(0,5000+1,100)
ntadv = len(ztadv)
tadv = np.zeros((2,ntadv),dtype=np.float64)

for iz in range(0,ntadv):
    if ztadv[iz] < 4000:
        tadv[:,iz] = -2.51/86400. + (-2.18+2.51)/(86400.*4000.)*ztadv[iz]
    elif ztadv[iz] < 5000.:
        tadv[:,iz] = -2.18/86400. + 2.18/(86400.*(5000.-4000.))*(ztadv[iz]-4000.)
    else:
        tadv[:,iz] = 0.0

timeAxes['tadv'] = utils.Axis('time_tadv',[t0,t1],name='Forcing time for tadv',units=tunits)
levAxes['tadv'] = utils.Axis('lev_tadv',ztadv,name='Altitude for variable tadv',units='m')
data['tadv'] = utils.Variable('tadv',name='Large-scale advection of temperature + radiative tendency',units='K s-1',data=tadv,level=levAxes['tadv'],time=timeAxes['tadv'],lat=latAxis,lon=lonAxis)

# Large-scale advection of specific humidity
zqvadv = range(0,5000+1,100)
nqvadv = len(zqvadv)
qvadv = np.zeros((2,nqvadv),dtype=np.float64)

for iz in range(0,nqvadv):
    if zqvadv[iz] < 3000:
        qvadv[:,iz] = -1.0/86400. + (0.345+1.0)/(86400.*3000.)*zqvadv[iz]
    elif zqvadv[iz] < 4000.:
        qvadv[:,iz] = 0.345/86400.
    elif zqvadv[iz] < 5000.:
        qvadv[:,iz] = 0.345/86400. - 0.345/(86400.*(5000.-4000.))*(zqvadv[iz]-4000.)
    else:
        qvadv[:,iz] = 0.0

timeAxes['qvadv'] = utils.Axis('time_qvadv',[t0,t1],name='Forcing time for qvadv',units=tunits)
levAxes['qvadv'] = utils.Axis('lev_qvadv',zqvadv,name='Altitude for variable qvadv',units='m')
data['qvadv'] = utils.Variable('qvadv',name='Large-scale advection of specific humidity',units='kg kg-1 s-1',data=qvadv/1000.,level=levAxes['qvadv'],time=timeAxes['qvadv'],lat=latAxis,lon=lonAxis)

###############################
# 4. Writing file
###############################


g = nc.Dataset('RICO_REF_orig.nc','w',format='NETCDF3_CLASSIC')

for var in ['ps','u','v','temp','qv','ug','vg','w','tadv','qvadv','ts']:
    print var
    #data[var].info()
    #data[var].plot(rep_images=rep_images)
    data[var].write(g)

g.Conventions = "CF-1.0" 
g.comment = "Forcing and initial conditions for RICO composite case - Original definition" 
g.reference = "http://projects.knmi.nl/rico/setup1d_composite.html" 
g.author = "R. Roehrig" 
g.modifications = ""
g.script = 'DEPHY-SCM/RICO/setup_orig.py'
g.history = "Created " + time.ctime(time.time())
g.case = "RICO/REF" 
g.startDate = startDate
g.endDate = endDate
g.tadv = 1 
g.qvadv = 1 
g.trad = "adv"
g.forc_omega = 0
g.forc_w = 1
g.forc_geo = 1
g.nudging_u = 0
g.nudging_v = 0
g.nudging_t = 0
g.nudging_qv = 0
g.zorog = zorog
g.surfaceType = "ocean"
g.surfaceForcing = "ts"

g.close()


