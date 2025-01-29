#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 21 Septembre 2021

@author: Romain Roehrig

"""

## EUROCS FIRE straotumulus case original case definition

from datetime import datetime, timedelta

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

tmin = datetime(1987, 7, 14, 8)
tmax = tmin + timedelta(hours=120)

case = Case('FIRE/REF',
        lat=33.3,
        lon=-119.5,
        startDate=tmin,
        endDate=tmax,
        surfaceType="ocean",
        zorog=0.)

case.set_title("Forcing and initial conditions for FIRE Reference case - Original definition")
case.set_reference("Dyunkerke et al. (2004, QJRMS), Chlond et al. (2004, QJRMS)")
case.set_author("R. Roehrig")
case.set_script("DEPHY-SCM/FIRE/REF/driver_DEF.py")

################################################
# 2. Initial state
################################################

# Surface pressure
ps = 101250.
case.add_init_ps(ps)

# Zonal and meridional wind
zwind = [   0, 1200. ]
u  =    [ 3.4 for z in zwind]
v  =    [-4.9 for z in zwind]

case.add_init_wind(u=u,ulev=zwind,v=v,vlev=zwind,levtype='altitude')

# Liquid-water potential temperature
zthetal = [  0.,   595., 605., 1200.]

def f_thetal(z):
    if z <= 595:
        return 287.5
    elif z <= 605:
        return 287.5+12.*(z-595.)/(605.-595.)
    else:
        return 287.5+12 + 0.0075*(z-605)

thetal  = [f_thetal(z) for z in zthetal]

case.add_init_thetal(thetal,lev=zthetal,levtype='altitude')

# Total water content
zqt = [  0.,   595., 605., 1200.]

def f_qt(z):
    if z <= 595:
        return 9.6
    elif z <= 605:
        return 9.6-3.0*(z-595.)/(605.-595.)
    else:
        return 9.6-3.0 - 0.003*(z-605)

qt = [f_qt(z)/1000. for z in zqt] # in kg kg-1

case.add_init_qt(qt,lev=zqt,levtype='altitude')

################################################
# 3. Forcing
################################################

# Constant Geostrophic wind across the simulation
case.add_geostrophic_wind(ug=u,uglev=zwind,vg=v,vglev=zwind,levtype='altitude')

# Large-scale velocity - constant
zw = [0., 1200.]
w  = [-1.e-5*z for z in zw]

case.add_vertical_velocity(w=w,lev=zw,levtype='altitude')

# Large-scale advection of temperature + radiative tendency - constant 
zthetal_adv = [ 0., 500., 1200.]
thetal_adv  = [-7.5e-8*max(z, 500) for z in zthetal_adv] # in K s-1

case.add_thetal_advection(thetal_adv,lev=zthetal_adv,levtype='altitude')

# Large-scale advection of specific humidity - constant
zqt_adv = [ 0., 500., 1200.]
qt_adv  = [3.0e-11*max(z,500) for z in zqt_adv] # in kg kg-1 s-1

case.add_qt_advection(qt_adv,lev=zqt_adv,levtype='altitude') # converted in kg kg-1 s-1 (array type required)

# Surface Forcing. Constant sea surface temperature
ts = 289.
case. add_forcing_ts(ts)

################################################
# 4. Writing file
################################################

case.write('FIRE_REF_DEF_driver.nc')

if lverbose:
    case.info()

################################################
# 5. Ploting, if asked
################################################

if lplot:
    case.plot(rep_images='./images/driver_DEF/',timeunits='hours')
