#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 11 July 2025

@author: Fleur Couvreux

## CASS original case definition: a shallow cumulus cases over the SGP obtained from a composite from a lot of fair-weather day observations
## Zhang et al, JAS, 2017 Large-Eddy Simulation of shallow cumulus over land: a composite case based on AROM long-term observations at its southern great plains site
## doi: 10.1175/JAS-D-16-0317.1
"""

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

case = Case('CASS/REF',
        lat=36.56,
        lon=-100.61,
        startDate="20010514120000",
        endDate="20010515030000",
        surfaceType='land',
        zorog=0.)

case.set_title("Forcing and initial conditions for CASE case - Meso-NH definition")
case.set_reference("Zhang et al. (2017, JAS)")
case.set_author("F. Couvreux")
case.set_script("DEPHY-SCM/CASS/REF/driver_DEF.py")

################################################
# 2. Initial state
################################################

# Surface pressure
ps = 97000.
case.add_init_ps(ps)

# Initial profiles
# Altitude above the ground
z = np.genfromtxt('init.txt',dtype=None,skip_header=1,usecols=0)

# Zonal and meridional wind
u = np.genfromtxt('init.txt',dtype=None,skip_header=1,usecols=4)
v = np.genfromtxt('init.txt',dtype=None,skip_header=1,usecols=5)

case.add_init_wind(u=u, v=v, lev=z, levtype='altitude')

# Potential temperature
theta = np.genfromtxt('init.txt',dtype=None,skip_header=1,usecols=2)

case.add_init_theta(theta, lev=z, levtype='altitude')

# Water vapor mixing ratio
rv = np.genfromtxt('init.txt',dtype=None,skip_header=1,usecols=3)

case.add_init_rv(rv, lev=z, levtype='altitude')

################################################
# 3. Forcing
################################################

# Atmospheric forcing
timeF=np.linspace(0.,14.,15)*3600.
ntf = len(timeF)
print('ntf',ntf)

zforc = np.genfromtxt('forcing.txt',dtype=None,skip_header=1,usecols=0)
nzf, = zforc.shape
print('nzf',nzf)
nzf = nzf//ntf
print('nzf',nzf)
zforc = zforc[0:nzf]

# Geostrophic wind
ug = np.genfromtxt('forcing.txt',dtype=None,skip_header=1,usecols=4)
ug = np.reshape(ug,(ntf,nzf))

vg = np.genfromtxt('forcing.txt',dtype=None,skip_header=1,usecols=5)
vg = np.reshape(vg,(ntf,nzf))

case.add_geostrophic_wind(ug=ug,vg=vg,time=timeF,lev=zforc,levtype='altitude')

# Vertical velocity
w = np.genfromtxt('forcing.txt',dtype=None,skip_header=1,usecols=6)
w = np.reshape(w,(ntf,nzf))

case.add_vertical_velocity(w=w,time=timeF,lev=zforc,levtype='altitude')

# Advection of potential temperature
theta_adv = np.genfromtxt('forcing.txt',dtype=None,skip_header=1,usecols=2)
theta_adv = np.reshape(theta_adv,(ntf,nzf))

case.add_theta_advection(theta_adv,time=timeF,lev=zforc,levtype='altitude',include_rad=True)

# Advection of water vapor mixing ratio
rv_adv = np.genfromtxt('forcing.txt',dtype=None,skip_header=1,usecols=3)
rv_adv = np.reshape(rv_adv,(ntf,nzf))

case.add_rv_advection(rv_adv,time=timeF,lev=zforc,levtype='altitude')

# Surface Forcing
timeSfc = np.genfromtxt('surface_forcing.txt',dtype=None,skip_header=1,usecols=0)
timeSfc=(timeSfc-timeSfc[0])*86400.
print('timeSfc',timeSfc)
shf = np.genfromtxt('surface_forcing.txt',dtype=None,skip_header=1,usecols=2)
lhf = np.genfromtxt('surface_forcing.txt',dtype=None,skip_header=1,usecols=3)
tau = np.genfromtxt('surface_forcing.txt',dtype=None,skip_header=1,usecols=4)
ts = np.genfromtxt('surface_forcing.txt',dtype=None,skip_header=1,usecols=1)
ustar = np.sqrt(tau)
print('shf',shf,'lhf',lhf,'ustar',ustar)

case.add_surface_fluxes(sens=shf,lat=lhf,time=timeSfc,forc_wind='ustar',ustar=ustar[0])
case.add_forcing_ts(ts,time=timeSfc)

################################################
# 4. Writing file
################################################

case.write('CASS_REF_DEF_driver.nc')

if lverbose:
    case.info()

################################################
# 5. Ploting, if asked
################################################

if lplot:
    case.plot(rep_images='./images/driver_DEF/',timeunits='hours')
