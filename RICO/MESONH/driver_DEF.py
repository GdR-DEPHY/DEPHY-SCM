#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on 30 November 2019

@author: Romain Roehrig
"""

## RICO composite short case, Meso-NH definition

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

case = Case('RICO/MESONH',
        lat=18,
        lon=-61.5,
        startDate="20041227000000",
        endDate="20041228000000",
        zorog=0.)

case.set_title("Forcing and initial conditions for RICO composite short case - Meso-NH definition")
case.set_reference("http://projects.knmi.nl/rico/setup1d_composite.html")
case.set_author("R. Roehrig, F. Couvreux")
case.set_script("DEPHY-SCM/RICO/MESONH/driver_DEF.py")

# time units are expected to be seconds since startDate
t0 = 0 # 00:00 UTC, 16 December 2004
t1 = 24*3600 # 72-hour long simulation


################################################
# 2. Initial state
################################################

# Surface pressure
ps = 101540.
case.add_variable('ps',[ps,])

# Zonal wind
zu = [   0, 4000. ]
u  = [-9.9,   -1.9]

case.add_variable('u',u,lev=zu,levtype='altitude')

# Meridional wind
zv = [ 0.,   4000. ]
v  = [-3.8,    -3.8]

case.add_variable('v',v,lev=zv,levtype='altitude')

# Potential temperature
ztheta = [  0.,   740., 3260.,   4000.]
theta  = [297.9,  297.9, 312.664, 317.]

case.add_variable('theta',theta,lev=ztheta,levtype='altitude')

# Specific humidity
zrv =[0.,     740.,     3260.,     4000.     ] 
rv = [0.01626,  0.01399,   0.00241,   0.00180]

case.add_variable('rv',np.array(rv),lev=zrv,levtype='altitude')

################################################
# 3. Forcing
################################################

# Constant Geostrophic wind across the simulation
zforc = [0., 2260., 2980., 4000.]
ug = [-9.9, -5.4, -3.9, -1.9]

case.add_variable('ug',np.array([ug,ug]),time=[t0,t1],lev=zforc,levtype='altitude')

vg = [-3.8, -3.8, -3.8, -3.8]

case.add_variable('vg',np.array([vg,vg]),time=[t0,t1],lev=zforc,levtype='altitude')

# Surface Forcing. Constant sea surface temperature
ts = 299.8

case.add_variable('ts',[ts,ts],time=[t0,t1])

# Large-scale velocity - constant
w  = [0.,   -0.005,    -0.005,   -0.005]

case.add_variable('w',[w,w],time=[t0,t1],lev=zforc,levtype='altitude')

# Large-scale advection of temperature + radiative tendency - constant 
thadv  = [-2.89e-5, -2.89e-5, -2.89e-5, -2.89e-5]

case.add_variable('theta_adv',np.array([thadv,thadv]),time=[t0,t1],lev=zforc,levtype='altitude')

# Large-scale advection of specific humidity - constant
rvadv  = [-1.16e-8, 0.02e-8, 0.40e-8, 0.40e-8]
case.add_variable('rv_adv',np.array([rvadv,rvadv]),time=[t0,t1],lev=zforc,levtype='altitude')


################################################
# 4. Attributes
################################################

# advection of theta and rt
case.set_attribute("adv_theta",1)
case.set_attribute("adv_rv",1)
# potential temperature radiative tendency is included in advection
case.set_attribute("rad_temp","adv")
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

case.write('RICO_MESONH_DEF_driver.nc',verbose=False)

if lverbose:
    case.info()

################################################
# 6. Ploting, if asked
################################################

if lplot:
    case.plot(rep_images='./images/driver_DEF/',timeunits='hours')
