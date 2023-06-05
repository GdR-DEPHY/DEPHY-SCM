#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 24 May 2023

@author: Fleur Couvreux

Modification
  2023/06/05, R. Roehrig, some fixes and cleaning
"""

## LBA/REF original case definition

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


case = Case('LBA/REF',
        lat=-8.,
        lon=-63.,
        startDate="19990223073000",
        endDate="19990223150000",
        surfaceType='land',
        zorog=0.)


case.set_title("Forcing and initial conditions for LBA case - Original definition")
case.set_reference(" https://rmets.onlinelibrary.wiley.com/doi/full/10.1256/qj.04.147; Grabowski et al. (2006, QJRMS)")
case.set_author("F. Couvreux")
case.set_script("DEPHY-SCM/LBA/REF/driver_DEF.py")

################################################
# 2. Initial state
################################################

# Surface pressure
ps = 99130.
case.add_init_ps(ps)

# Initial profile in init.txt
# Altitude above the ground
z = np.genfromtxt('init.txt',dtype=np.float32,usecols=0)

# Zonal and meridional wind
u = np.genfromtxt('init.txt',dtype=np.float32,usecols=3)
v = np.genfromtxt('init.txt',dtype=np.float32,usecols=4)

case.add_init_wind(u=u, v=v, lev=z, levtype='altitude')

# Potential temperature
theta = np.genfromtxt('init.txt',dtype=np.float32,usecols=1)

case.add_init_theta(theta, lev=z, levtype='altitude')

# Water vapor mixing ratio
rv = np.genfromtxt('init.txt',dtype=np.float32,usecols=2)

case.add_init_rv(rv, lev=z, levtype='altitude')

################################################
# 3. Forcing
################################################

# Forcing time axis
timeForc = [0.,3600.,7200.,10800.,14400.,18000.,21600.]
nt = len(timeForc)

height_forc = np.genfromtxt('LBA_formatcommun_ZFR_1',dtype=np.float32,usecols=0)
nlev, = height_forc.shape

data = {}
for var in ['u','v','dthdt']:
    data[var] = np.zeros((nt, nlev), dtype=np.float32)

for it in range(0,nt):
    fin = f'LBA_formatcommun_ZFR_{it+1}'
    data['u'][it,:] = np.genfromtxt(fin,dtype=np.float32,usecols=1)
    data['v'][it,:] = np.genfromtxt(fin,dtype=np.float32,usecols=2)
    data['dthdt'][it,:] = np.genfromtxt(fin,dtype=np.float32,usecols=6)

# Potential temperature advection, which includes radiation
case.add_theta_advection(data['dthdt'],include_rad=True,lev=height_forc,levtype='altitude',time=timeForc)

# Wind nudging
case.add_wind_nudging(unudg=data['u'],vnudg=data['v'],timescale=3600.,time=timeForc,lev=height_forc,levtype='altitude')

# Surface fluxes
timeflux = np.genfromtxt('flux_LBA',dtype=np.float32,usecols=0)
sens = np.genfromtxt('flux_LBA',dtype=np.float32,usecols=1)
flat = np.genfromtxt('flux_LBA',dtype=np.float32,usecols=2)
case.add_surface_fluxes(sens,flat,time=timeflux,forc_wind='z0',z0=0.035)

################################################
# 4. Writing file
################################################

case.write('LBA_REF_DEF_driver.nc')

if lverbose:
    case.info()

################################################
# 6. Ploting, if asked
################################################

if lplot:
    case.plot(rep_images='./images/driver_DEF/',timeunits='hours')
