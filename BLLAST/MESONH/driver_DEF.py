#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on 26 August 2021

@author: Fleur Couvreux

Modifications
  2022/07/04, Romain Roehrig, Make start/end dates consistent with Darbieu et al.
                              Add P2OA altitude, and a consistent surface pressure (from ERA5)

"""

## BLLAST original case definition: a reaslitic convective boundary layer over Southern of France
## Darbieu et al, ACP, 2015 Turbulence vertical structure of the boundary layer during the 
## afternoon transition doi:10.5194/acp-15-10071-2015

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

case = Case('BLLAST/MESONH',
        lat=43.1,
        lon=0.21,
        startDate="20110620060000",
        endDate="20110620180000",
        surfaceType='land',
        zorog=0.)

case.set_title("Forcing and initial conditions for BLLAST case - Meso-NH definition")
case.set_reference("Darbieu et al. (2015, ACP)")
case.set_author("F. Couvreux")
case.set_script("DEPHY-SCM/BLLAST/MESONH/driver_REF.py")

################################################
# 2. Initial state
################################################

# Surface pressure
ps = 96000. # Approximate value from ERA5
case.add_init_ps(ps)


# Thermodynamical initial profiles
# Altitude above the ground
z = np.genfromtxt('Initial_profile.txt',dtype=None,skip_header=6,usecols=0)
u = np.genfromtxt('Initial_profile.txt',dtype=None,skip_header=6,usecols=3)
v = np.genfromtxt('Initial_profile.txt',dtype=None,skip_header=6,usecols=4)
theta = np.genfromtxt('Initial_profile.txt',dtype=None,skip_header=6,usecols=1)
rv = np.genfromtxt('Initial_profile.txt',dtype=None,skip_header=6,usecols=2)

# Wind initial profiles
case.add_init_wind(u=u, v=v, lev=z, levtype='altitude')


# Potential temperature
case.add_init_theta(theta, lev=z, levtype='altitude')

# Water vapor mixing ratio
rv = rv*0.001 # conversion g/kg en kg/kg

case.add_init_rv(rv, lev=z, levtype='altitude')

################################################
# 3. Forcing
################################################
zadv=[100.,225.,350.,475.,600.,725.,850.,975,1100.,1225.,1350.,1475.,1600.,1725.,1850.,1975.,2100.,2225.,2350.,2475.]
#print('z adv=',zadv)
timeadv=[0.,3600.,7200.,10800.,14400.,18000.,21600.,25200.,28800.,32400.,36000.,39600.,43200.,46800.]

advth=np.genfromtxt('advection_profiles_Theta.txt',dtype=None,skip_header=8)
#print('size advth',advth.shape)
advth=advth/86400. #conversion K/day en K/s
advrv=np.genfromtxt('advection_profiles_qv.txt',dtype=None,skip_header=8)
#print('size advrv',advrv.shape)
advrv=advrv/86400./1000. #conversion g/kg/day en kg/kg/s


# No advection
#case.add_theta_advection(advth,time=timeadv,lev=zadv,levtype='altitude')
#case.add_rv_advection(advrv,time=timeadv,lev=zadv,levtype='altitude')

# No radiation
case.deactivate_radiation()

# Surface Forcing
timeSfc = np.genfromtxt('Surface_flux.txt',dtype=None,skip_header=6,usecols=0)
timeref=timeSfc[0]*3600.
for it in range(0,timeSfc.shape[0]):
    timeSfc[it]=timeSfc[it]*3600.-timeref
#print('timeSfc',timeSfc)

shf = np.genfromtxt('Surface_flux.txt',dtype=None,skip_header=6,usecols=1)
lhf = np.genfromtxt('Surface_flux.txt',dtype=None,skip_header=6,usecols=2)

case.add_surface_fluxes(sens=shf,lat=lhf,time=timeSfc,forc_wind='z0',z0=0.1)
# pas d'info sur la rugosité on a laissé celui de IHOP

################################################
# 4. Writing file
################################################

case.write('BLLAST_MESONH_DEF_driver.nc')

if lverbose:
    case.info()

################################################
# 5. Ploting, if asked
################################################

if lplot:
    case.plot(rep_images='./images/driver_DEF/',timeunits='hours')
