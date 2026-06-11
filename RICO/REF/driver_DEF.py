#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 30 November 2019

@author: Romain Roehrig

Modification
  2021/01/03, R. Roehrig: update for improved case definition interface.
  2026/01/26, N. Villefranque: merge with MESONH and clean for publication.
"""

## RICO original case definition 

import numpy as np

from datetime import datetime, timedelta
from dephycf.Case import Case

################################################
# 0. General configuration of the present script
################################################

lplot = True # plot all the variables
lverbose = False # print information about variables and case

################################################
# 1. General information about the case
################################################

# RICO is a composite case from observations acquired during the RICO campaign
# between 2004/12/16 and 2005/01/08

duration=24
tmin = datetime(2004, 12, 23, 00)
tmax = tmin + timedelta(hours=duration)

case = Case('RICO/REF',
        lat=18,
        lon=-61.5,
        startDate=tmin,
        endDate=tmax,
        surfaceType="ocean",
        zorog=0.)

case.set_title("Forcing and initial conditions for RICO composite case - Original case")
case.set_reference("VanZanten et al. (2011, JAMES)")
case.set_author("R. Roehrig, F. Couvreux")
case.set_script("DEPHY-SCM/RICO/REF/driver_DEF.py")
case.set_modifications("No surface scheme prescribed (only SST)")

################################################
# 2. Initial state
################################################

# Surface pressure
ps = 101540.
case.add_init_ps(ps)

# Zonal and meridional wind
zu = zv = [   0, 4000.]
u  = [-9.9,   -1.9]
v  = [-3.8,    -3.8]

case.add_init_wind(u=u,ulev=zu,v=v,vlev=zv,levtype='altitude')

# Potential temperature
ztheta = [  0.,   740., 3260.,   4000.]
theta  = [297.9,  297.9, 312.664, 317.]

case.add_init_theta(theta,lev=ztheta,levtype='altitude')

# Specific humidity
## in the paper, qv values are provided, here converted to rv
zrv =[0.,     740.,     3260.,     4000.     ] 
rv = [0.01626,  0.01399,   0.00241,   0.00180] # in kg kg-1 

case.add_init_rv(rv,lev=zrv,levtype='altitude') 

################################################
# 3. Forcing
################################################

# Constant Geostrophic wind across the simulation
zforc = [0., 2260., 2980., 4000.]

ug = [-9.9, -5.4, -3.9, -1.9]
vg = [-3.8, -3.8, -3.8, -3.8]
case.add_geostrophic_wind(ug=ug,vg=vg,lev=zforc,levtype='altitude')

# Surface Forcing. Constant sea surface temperature
ts = 299.8
case.add_forcing_ts(ts)

# Large-scale velocity - constant in time
w  = [0.,   -0.005,    -0.005,   -0.005]
case.add_vertical_velocity(w=w,lev=zforc,levtype='altitude')

# Large-scale advection of temperature + radiative tendency - constant 
thadv  = [-2.89e-5, -2.89e-5, -2.89e-5, -2.89e-5] # in K s-1 
case.add_theta_advection(thadv,lev=zforc,levtype='altitude',include_rad=True) 

# Large-scale advection of specific humidity - constant in time
# in the paper, dqv/dt values are provided, here converted to drv/dt
rvadv  = [-1.16e-8, 0.02e-8, 0.40e-8, 0.40e-8] # in kg kg-1 s-1
case.add_rv_advection(rvadv,lev=zforc,levtype='altitude')

################################################
# 4. Writing file
################################################

case.write('RICO_REF_DEF_driver.nc')

if lverbose:
    case.info()

################################################
# 5. Ploting, if asked
################################################

if lplot:
    case.plot(rep_images='./images/driver_DEF/',timeunits='hours')
