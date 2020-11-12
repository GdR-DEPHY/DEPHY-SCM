#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on 30 November 2019

@author: Romain Roehrig

Modification
  2020/11/12, R. Roehrig: update for improved case definition interface.
"""

## RICO Composite short case original case definition

import sys
sys.path = ['../../utils/',] + sys.path

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

case = Case('RICO/SHORT',
        lat=18,
        lon=-61.5,
        startDate="20041216000000",
        endDate="20041219000000",
        surfaceType="ocean",
        zorog=0.)

case.set_title("Forcing and initial conditions for RICO composite short case - Original definition")
case.set_reference("http://projects.knmi.nl/rico/setup1d_composite.html")
case.set_author("R. Roehrig")
case.set_script("DEPHY-SCM/RICO/SHORT/driver_DEF.py")

################################################
# 2. Initial state
################################################

# Surface pressure
ps = 101540.
case.add_init_ps(ps)

# Zonal and meridional wind
zu = [   0, 4000., 12000.,13000., 20000.]
u  = [-9.9,   -1.9,   30.,   30.,     0.]

zv = [ 0.,   20000. ]
v  = [-3.8,     -3.8]

case.add_init_wind(u=u,ulev=zu,v=v,vlev=zv,levtype='altitude')

# Temperature
ztemp = [  0.,   740., 4000., 15000., 17500., 20000., 60000.]
temp  = [299.2,  292.,  278.,   203.,   194.,  206.,   270.]

case.add_init_temp(temp,lev=ztemp,levtype='altitude')

# Specific humidity
# put 9000 m instead of 10000 m which seems more relevant 
# according to the formula at http://projects.knmi.nl/rico/setup1d_composite.html
zqv =[ 0., 740., 3260., 4000., 9000.] 
qv = [16.,  13.8,   2.4,   1.8,   0.] # in g kg-1

case.add_init_qv(np.array(qv)/1000.,lev=zqv,levtype='altitude') # converted in kg kg-1 (array type required)

################################################
# 3. Forcing
################################################

# Constant Geostrophic wind across the simulation
case.add_geostrophic_wind(ug=u,uglev=zu,vg=v,vglev=zv,levtype='altitude')

# Large-scale velocity - constant
zw = [0., 2260,      4000.,   5000.]
w  = [0.,   -0.005,    -0.005,   0.]

case.add_vertical_velocity(w=w,lev=zw,levtype='altitude')

# Large-scale advection of temperature + radiative tendency - constant 
ztadv = [ 0.,  4000.,  5000.]
tadv  = [-2.51,  -2.18,   0.] # in K day-1

case.add_temp_advection(np.array(tadv)/86400.,lev=ztadv,levtype='altitude',include_rad=True) # converted in K s-1 (array type required)

# Large-scale advection of specific humidity - constant
zqvadv = [ 0.,  3000.,   4000.,   5000.]
qvadv  = [-1.0,    0.345,   0.345,   0.] # in g kg-1 day-1

case.add_qv_advection(np.array(qvadv)/86400./1000.,lev=zqvadv,levtype='altitude') # converted in kg kg-1 s-1 (array type required)

# Surface Forcing. Constant sea surface temperature
ts = 299.8
case. add_forcing_ts(ts)

################################################
# 4. Writing file
################################################

case.write('RICO_SHORT_DEF_driver.nc',verbose=False)

if lverbose:
    case.info()

################################################
# 6. Ploting, if asked
################################################

if lplot:
    case.plot(rep_images='./images/driver_DEF/',timeunits='hours')
