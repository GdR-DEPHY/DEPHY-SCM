#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on 5 November 2020

@author: Fleur Couvreux

Modification
  2021/01/06, R. Roehrig: update for improved case definition interface.
"""

## SCMS original case definition
## From Neggers et al 2003

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

case = Case('SCMS/REF',
        lat=28.7,
        lon=-81.,
        startDate="19950805120000",
        endDate="19950806000000",
        surfaceType='land',
        zorog=0.)

case.set_title("Forcing and initial conditions for SCMS-Cumulus case - Meso-NH definition")
case.set_reference(" Neggers et al. (2003, QJRMS)")
case.set_author("F. Couvreux, N. Villefranque")
case.set_script("DEPHY-SCM/SCMS/MESONH/driver_MESONH.py")

################################################
# 2. Initial state
################################################

# Surface pressure
ps = 102100.
case.add_init_ps(ps)

#case.add_variable('ps',[ps,])

# Wind initial profiles
#         z (m) u (m s-1) v (m s-1)
init = [  0.0,      -4.0,     4.0,\
       4980.0,      -4.0,     4.0]

case.add_init_wind(u=init[1::3], v=init[2::3], lev=init[0::3], levtype='altitude')

# Thermodynamical initial profiles
# Altitude above the ground
z = np.genfromtxt('init.txt',dtype=None,skip_header=1,usecols=0)

# Potential temperature
theta = np.genfromtxt('init.txt',dtype=None,skip_header=1,usecols=1)

case.add_init_theta(theta, lev=z, levtype='altitude')

# Water vapor mixing ratio
rv = np.genfromtxt('init.txt',dtype=None,skip_header=1,usecols=2)

case.add_init_rv(rv, lev=z, levtype='altitude')

################################################
# 3. Forcing
################################################

# Constant Geostrophic wind across the simulation
zforc = [0., 1500., 5000.]
ug = [-4., -4., -4.]
vg = [ 4.,  4.,  4.]
case.add_geostrophic_wind(ug=ug,vg=vg,lev=zforc,levtype='altitude')

# Potential temperature and water vapor mixing ratio advection
zadv     =  [ 0.,       1500., 3000., 5000.]
theta_adv = [-0.0000347,   0.,     0.,   0.]
rv_adv    = [0.,           0.,     0.,   0.]

case.add_theta_advection(theta_adv,lev=zadv,levtype='altitude',include_rad=True)
case.add_rv_advection(rv_adv,lev=zadv,levtype='altitude')


# Surface Forcing
timeSfc = np.genfromtxt('surface_forcing.txt',dtype=None,skip_header=1,usecols=0)
shf = np.genfromtxt('surface_forcing.txt',dtype=None,skip_header=1,usecols=1)
lhf = np.genfromtxt('surface_forcing.txt',dtype=None,skip_header=1,usecols=2)

case.add_surface_fluxes(sens=shf,lat=lhf,time=timeSfc,forc_wind='z0',z0=0.035)

################################################
# 4. Writing file
################################################

case.write('SCMS_REF_DEF_driver.nc',verbose=False)

if lverbose:
    case.info()

################################################
# 5. Ploting, if asked
################################################

if lplot:
    case.plot(rep_images='./images/driver_DEF/',timeunits='hours')
