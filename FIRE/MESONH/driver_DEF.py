#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 21 Septembre 2021

@author: Romain Roehrig

"""

## EUROCS FIRE straotumulus case original case definition

import os

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

case = Case('FIRE/MESONH',
        lat=33.3,
        lon=-119.5,
        startDate="19870714080000",
        endDate="19870715210000",
        surfaceType="ocean",
        zorog=0.)

case.set_title("Forcing and initial conditions for FIRE Reference case - Original definition")
case.set_reference("Dyunkerke et al. (2004, QJRMS) Brient et al (2019,GRL) modifié en haut")
case.set_author("F. Couvreux")
case.set_script("DEPHY-SCM/FIRE/MESONH/driver_DEF.py")

################################################
# 2. Initial state
################################################

# Surface pressure
ps = 101250.
case.add_init_ps(ps)

# Zonal and meridional wind
zwind = [   0, 1200. ]
u  =    [ 3.4415 for z in zwind]
v  =    [-4.9149 for z in zwind]

case.add_init_wind(u=u,ulev=zwind,v=v,vlev=zwind,levtype='altitude')

# Liquid-water potential temperature
zthetal = [  0., 595., 605., 650., 800., 1000., 1200.]
thetal  = [ 287.5, 287.5, 299.57, 299.91, 301.03, 302.54, 304.04]

case.add_init_thetal(thetal,lev=zthetal,levtype='altitude')

# Total water content
zqt = [  0., 595., 605., 650., 800., 1000., 1200.]
qt  = [  9.6e-3, 9.6e-3, 6.57e-3, 6.43e-3, 5.99e-3, 5.39e-3, 4.79e-3]

case.add_init_qt(qt,lev=zqt,levtype='altitude')

################################################
# 3. Forcing
################################################

# Constant Geostrophic wind across the simulation
case.add_geostrophic_wind(ug=u,uglev=zwind,vg=v,vglev=zwind,levtype='altitude')

# Large-scale velocity - constant
zw = [0., 100., 300., 500., 595., 605., 650., 800., 900., 1050., 1100., 1200.]
w  = [0., -0.001, -0.003, -0.005, -0.00595, -0.00605, -0.0065, -0.008, -0.009, -0.0105, -0.011,0.]

case.add_vertical_velocity(w=w,lev=zw,levtype='altitude')

# Large-scale advection of potential temperature  
zthetal_adv = [0., 500., 595., 605., 650., 800., 900., 1050., 1100., 1200.]
thetal_adv  = [-3.75e-5, -3.75e-5, -4.462e-5, -4.537e-5, -4.875e-5, -6.e-5, -6.75e-5, -7.875e-5,-8.25e-5,0.]

case.add_thetal_advection(thetal_adv,lev=zthetal_adv,levtype='altitude')

# Large-scale advection of specific humidity - constant
zqt_adv = [0., 500., 595., 605., 650., 800., 900., 1050., 1100., 1200.]
qt_adv  = [1.5e-8, 1.5e-8, 1.785e-8, 1.815e-8, 1.95e-8, 2.4e-8, 2.7e-8, 3.15e-8,3.3e-8,0.]

case.add_qt_advection(qt_adv,lev=zqt_adv,levtype='altitude') # converted in kg kg-1 s-1 (array type required)

# Surface Forcing. Constant sea surface temperature
ts = 289.
case.add_forcing_ts(ts)

################################################
# 4. Writing file
################################################

case.write('FIRE_MESONH_DEF_driver.nc')

if lverbose:
    case.info()

################################################
# 5. Ploting, if asked
################################################

if lplot:
    case.plot(rep_images='./images/driver_DEF/',timeunits='hours')
