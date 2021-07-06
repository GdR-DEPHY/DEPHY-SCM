#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on 30 October 2020

@author: Fleur Couvreux

Modification
  2021/01/06, R. Roehrig: update for improved case definition interface.
"""

## IHOP original case definition: a realistic convective boundary layer growth over Oklahoma
## Couvreux et al, QJRMS, 2005 Water-vapour variability within a convective boundary-layer assessed by
## large-eddy simulations and IHOP 2002 observations
## doi: 10.1256/qj.04.167

import os
import sys
sys.path = ['../../utils/',] + sys.path

import netCDF4 as nc
import numpy as np

import constants
from Case import Case

################################################
# 0. General configuration of the present script
################################################

lplot = True # plot all the variables
lverbose = False # print information about variables and case

################################################
# 1. General information about the case
################################################

case = Case('IHOP/REF',
        lat=36.56,
        lon=-100.61,
        startDate="20020614120000",
        endDate="20020614190000",
        surfaceType='land',
        zorog=0.)

case.set_title("Forcing and initial conditions for IHOP case - Meso-NH definition")
case.set_reference("Couvreux et al. (2005, QJRMS)")
case.set_author("F. Couvreux")
case.set_script("DEPHY-SCM/IHOP/REF/driver_DEF.py")

################################################
# 2. Initial state
################################################

# Surface pressure
ps = 91800.
case.add_init_ps(ps)

# Initial profiles
# Altitude above the ground
z = np.genfromtxt('init.txt',dtype=None,skip_header=1,usecols=0)

# Zonal and meridional wind
u = np.genfromtxt('init.txt',dtype=None,skip_header=1,usecols=1)
v = np.genfromtxt('init.txt',dtype=None,skip_header=1,usecols=2)

case.add_init_wind(u=u, v=v, lev=z, levtype='altitude')

# Potential temperature
theta = np.genfromtxt('init.txt',dtype=None,skip_header=1,usecols=3)

case.add_init_theta(theta, lev=z, levtype='altitude')

# Water vapor mixing ratio
rv = np.genfromtxt('init.txt',dtype=None,skip_header=1,usecols=4)

case.add_init_rv(rv, lev=z, levtype='altitude')

################################################
# 3. Forcing
################################################

# Atmospheric forcing
timeF= [0.,10800.,21600.]
ntf = len(timeF)

zforc = np.genfromtxt('atm_forcing.txt',dtype=None,skip_header=1,usecols=0)
nzf, = zforc.shape
nzf = nzf/ntf
zforc = zforc[0:nzf]

# Geostrophic wind
ug = np.genfromtxt('atm_forcing.txt',dtype=None,skip_header=1,usecols=1)
ug = np.reshape(ug,(ntf,nzf))

vg = np.genfromtxt('atm_forcing.txt',dtype=None,skip_header=1,usecols=2)
vg = np.reshape(vg,(ntf,nzf))

case.add_geostrophic_wind(ug=ug,vg=vg,time=timeF,lev=zforc,levtype='altitude')

# Vertical velocity
w = np.genfromtxt('atm_forcing.txt',dtype=None,skip_header=1,usecols=3)
w = np.reshape(w,(ntf,nzf))

case.add_vertical_velocity(w=w,time=timeF,lev=zforc,levtype='altitude')

# Advection of potential temperature
theta_adv = np.genfromtxt('atm_forcing.txt',dtype=None,skip_header=1,usecols=4)
theta_adv = np.reshape(theta_adv,(ntf,nzf))

case.add_theta_advection(theta_adv,time=timeF,lev=zforc,levtype='altitude',include_rad=True)

# Advection of water vapor mixing ratio
rv_adv = np.genfromtxt('atm_forcing.txt',dtype=None,skip_header=1,usecols=5)
rv_adv = np.reshape(rv_adv,(ntf,nzf))

case.add_rv_advection(rv_adv,time=timeF,lev=zforc,levtype='altitude')

# Surface Forcing
timeSfc = np.genfromtxt('surface_forcing.txt',dtype=None,skip_header=1,usecols=0)
shf = np.genfromtxt('surface_forcing.txt',dtype=None,skip_header=1,usecols=1)
lhf = np.genfromtxt('surface_forcing.txt',dtype=None,skip_header=1,usecols=2)

case.add_surface_fluxes(sens=shf,lat=lhf,time=timeSfc,forc_wind='z0',z0=0.1)

################################################
# 4. Writing file
################################################

case.write('IHOP_REF_DEF_driver.nc')

if lverbose:
    case.info()

################################################
# 5. Ploting, if asked
################################################

if lplot:
    case.plot(rep_images='./images/driver_DEF/',timeunits='hours')
