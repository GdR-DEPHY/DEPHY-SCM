#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 22 February 2024

@author: Romain Roehrig

Modifications
"""

## Revisiting the BLLAST original case defined in 
## Darbieu et al, ACP, 2015 Turbulence vertical structure of the boundary layer during the 
## afternoon transition doi:10.5194/acp-15-10071-2015

## No advection
## Revisit of idealized initial profiles

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

case = Case('BLLAST/B2024',
        lat=43.1,
        lon=0.36,
        startDate="20110620050000",
        endDate="20110620180000",
        surfaceType='land',
        zorog=588.)

case.set_title("Forcing and initial conditions for BLLAST case - B2024 definition")
case.set_reference("Darbieu et al. (2015, ACP); Bernard et al. (2024)")
case.set_author("R. Roehrig")
case.set_script("DEPHY-SCM/BLLAST/B2024/driver_DEF.py")

################################################
# 2. Initial state
################################################

# Surface pressure
ps = 95000. # Approximate value from ERA5
case.add_init_ps(ps)


# Wind initial profiles
# Altitude above the ground
zu = [0., 3500., 5000., 14000., 17500.]
u =  [0.,   15.,   18.,    18.,     5.]
zv = [0.,  400., 4000., 5000., 17500.]
v =  [0.5,   0.,   -4.,    0.,     0.]

case.add_init_wind(u=u, v=v, ulev=zu, vlev=zv, levtype='altitude')

# Potential temperature initial profile
ztheta = [ 0., 65.,  85., 150., 500.,  1500., 2500., 4000., 11000., 17500.]
theta =  [20., 21., 23.5,  24.,  25.5,   32.,   35.,   42.,    66.,   155.] # in degrees Celcius
theta = [x+273.15 for x in theta]
case.add_init_theta(theta, lev=ztheta, levtype='altitude')

# Relative humidity
zhur = [ 0., 45., 85., 150., 550., 6000., 12000., 17500.]
hur =  [78., 72., 62.,  58.,  50.,   50.,     0.,     0.] 
hur = [x/100. for x in hur]

case.add_init_hur(hur, lev=zhur, levtype='altitude')

################################################
# 3. Forcing
################################################
# No advection

# No wind forcing

# No radiation
case.deactivate_radiation()

# Surface Forcing
timeSfc = np.genfromtxt('Surface_flux.txt',dtype=None,skip_header=6,usecols=0)
timeref=timeSfc[0]*3600.
for it in range(0,timeSfc.shape[0]):
    timeSfc[it]=timeSfc[it]*3600.-timeref

shf = np.genfromtxt('Surface_flux.txt',dtype=None,skip_header=6,usecols=1)
lhf = np.genfromtxt('Surface_flux.txt',dtype=None,skip_header=6,usecols=2)

# No information on the roughness length. Take some default (similar to IHOP)
case.add_surface_fluxes(sens=shf,lat=lhf,time=timeSfc,forc_wind='z0',z0=0.1)

################################################
# 4. Writing file
################################################

case.write('BLLAST_B2024_DEF_driver.nc')

if lverbose:
    case.info()

################################################
# 5. Ploting, if asked
################################################

if lplot:
    case.plot(rep_images='./images/driver_DEF/',timeunits='hours')
