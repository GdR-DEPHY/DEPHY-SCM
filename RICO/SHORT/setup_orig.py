#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on 30 November 2019

@author: Romain Roehrig
"""

## ARM-Cumulus original case definition
## From http://projects.knmi.nl/eurocs/ARM/case_ARM_html/

import os
import sys
sys.path = ['../../utils/',] + sys.path

import netCDF4 as nc
import numpy as np

import SCM_utils as utils

################################################
# 0. General configuration of the present script
################################################

lplot = True

################################################
# 1. General information about the case
################################################

case = utils.Case('RICO/SHORT',
        lat=18,
        lon=-61.5,
        startDate="20041216000000",
        endDate="20041219000000",
        zorog=0.)

case.set_comment("Forcing and initial conditions for RICO composite short case - Original definition")
case.set_reference("http://projects.knmi.nl/rico/setup1d_composite.html")
case.set_author("R. Roehrig")
case.set_script("DEPHY-SCM/RICO/SHORT/setup_orig.py")


# time units are expected to be seconds since startDate
t0 = 0 # 00:00 UTC, 16 December 2004
t1 = 72*3600 # 72-hour long simulation


################################################
# 2. Initial state
################################################

# Surface pressure
ps = 101540.
case.add_variable('ps',[ps,])

# Zonal wind
zu = [   0, 4000., 12000.,13000., 20000.]
u  = [-9.9,   -1.9,   30.,   30.,     0.]

case.add_variable('u',u,lev=zu,levtype='altitude')

# Meridional wind
zv = [ 0.,   20000. ]
v  = [-3.8,     -3.8]

case.add_variable('v',v,lev=zv,levtype='altitude')

# Temperature
ztemp = [  0.,   740., 4000., 15000., 17500., 20000., 60000.]
temp  = [299.2,  292.,  278.,   203.,   194.,  206.,   270.]

case.add_variable('temp',temp,lev=ztemp,levtype='altitude')

# Specific humidity
# put 9000 m instead of 10000 m which seems more relevant 
# according to the formula at http://projects.knmi.nl/rico/setup1d_composite.html
zqv =[ 0., 740., 3260., 4000., 9000.] 
qv = [16.,  13.8,   2.4,   1.8,   0.] # in g kg-1

case.add_variable('qv',np.array(qv)/1000.,lev=zqv,levtype='altitude') # converted in kg kg-1 (array type required)

################################################
# 3. Forcing
################################################

# Constant Geostrophic wind across the simulation
ug = np.zeros((2,len(u)),dtype=np.float64)
for i in range(0,2):
    ug[i,:] = u[:]

case.add_variable('ug',ug,time=[t0,t1],lev=zu,levtype='altitude')

vg = np.zeros((2,len(v)),dtype=np.float64)
for i in range(0,2):
    vg[i,:] = v[:]

case.add_variable('vg',vg,time=[t0,t1],lev=zv,levtype='altitude')

# Surface Forcing. Constant sea surface temperature
ts = 299.8

case.add_variable('ts',[ts,ts],time=[t0,t1])

# Large-scale velocity - constant
zw = [0., 2260,      4000.,   5000.]
w  = [0.,   -0.005,    -0.005,   0.]

case.add_variable('w',[w,w],time=[t0,t1],lev=zw,levtype='altitude')

# Large-scale advection of temperature + radiative tendency - constant 
ztadv = [ 0.,  4000.,  5000.]
tadv  = [-2.51,  -2.18,   0.] # in K day-1

case.add_variable('tadv',np.array([tadv,tadv])/86400.,time=[t0,t1],lev=ztadv,levtype='altitude') # converted in K s-1 (array type required)

# Large-scale advection of specific humidity - constant
zqvadv = [ 0.,  3000.,   4000.,   5000.]
qvadv  = [-1.0,    0.345,   0.345,   0.] # in g kg-1 day-1

case.add_variable('qvadv',np.array([qvadv,qvadv])/86400./1000.,time=[t0,t1],lev=zqvadv,levtype='altitude') # converted in kg kg-1 s-1 (array type required)


################################################
# 4. Attributes
################################################

# advection of theta and rt
case.set_attribute("tadv",1)
case.set_attribute("qvadv",1)
# potential temperature radiative tendency is included in advection
case.set_attribute("trad","adv")
# Geostrophic wind forcing
case.set_attribute("forc_geo",1)
# Vertical velocity
case.set_attribute("forc_w",1)
# Surface flux forcing, wind stress is computed using z0
case.set_attribute("surfaceType","ocean")
case.set_attribute("surfaceForcing","ts")

################################################
# 5. Writing file
################################################

case.write('RICO_SHORT_orig.nc',verbose=False)

case.info()

################################################
# 6. Ploting, if asked
################################################

if lplot:
    case.plot(rep_images='./images/setup_orig/',timeunits='hours')
