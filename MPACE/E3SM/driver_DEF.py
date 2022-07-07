#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on 01 October 2020

@author: Etienne Vignon

Modifications:
  2021/01/03, R. Roehrig: update for improved case definition interface.
  2021/01/25, E. Vignon: change subcase name since this case is not the "def" one but forced by an intermedirary file from the E3SM library

MPACE case, version forced by an intermediary file from the E3SM library
From Klein et al. 2009; QJRMS, DOI: 10.1002/qj.416
In  the  baseline  simulation, longitude is 209.0, latitude is 71.75. The lower initial boundary condition is specified as an ocean surface with temperature 274.01 K Models were asked to simulate the 12h starting from 1700 UTC 9 October 2004. Initial profiles of ice-liquid water temperature and total water are prescribed and correspond to a cloudy (purely liquid) convective boundary layer topped by an inversion. Surfaces heat and latent fluxes are assumed constant throughout the simulation and vertical velocity and horizontal advections of heat and water vapor are prescribed.
"""

import os
import sys
sys.path = ['../../utils/',] + sys.path

import netCDF4 as nc
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
# This is an idealized case so date is arbritrary
# lat/lon are fixed to a lat representative of Arctic conditions
# 8h duration with a constant surface cooling

case = Case('MPACE/E3SM',
        lat=71.75,
        lon=209.00,
        startDate="20041009170000",
        endDate="20041010050000",
        surfaceType="ocean",
        zorog=0.)

case.set_title("Forcing and initial conditions for MPACE case - Original definition")
case.set_reference("Klein et al. (2009, QJRMS)")
case.set_author("E. Vignon")
case.set_script("driver_DEF.py")
case.set_comment("Use of forcing file from E3SM, https://github.com/E3SM-Project/scmlib/wiki/E3SM-Single-Column-Model-Case-Library")

################################################
# 2. Input netCDF file
################################################

fin = nc.Dataset('MPACE_forcing_dephy.nc','r')

################################################
# 3. Initial state
################################################

# Surface pressure
ps  =  fin['Ps'][0,0,0]
case.add_init_ps(ps)

# Pressure
pressure  = fin['lev'][-1:0:-1]
pressure[0] = ps

case.add_init_pressure(pressure,lev=pressure,levtype='pressure',levid='lev')

# Temperature in K and moisture (vapor and total mixing ratio)
# Note that in the ref paper, the ice-liquid potential temperature is initially set 
# but the initial state corresponds to a purely-liquid cloud

temp  = fin['T'][0,-1:0:-1,0,0]

case.add_init_temp(temp,lev=pressure,levtype='pressure',levid='lev')

#  water vapor mixing ratio in kg/kg
rv =  fin['q'][0,-1:0:-1,0,0]
qv=rv/(1.+rv)

case.add_init_qv(qv,lev=pressure,levtype='pressure',levid='lev')


#  total water mixing ratio in kg/kg
rl =  fin['CLDLIQ'][0,-1:0:-1,0,0]
ri =  fin['CLDICE'][0,-1:0:-1,0,0]
rt=rl+ri+rv
qt=rt/(1.+rt)

case.add_init_qt(qt,lev=pressure,levtype='pressure',levid='lev')

# Zonal and meridional wind
u  = fin['u'][0,-1:0:-1,0,0]
v  = fin['v'][0,-1:0:-1,0,0]

case.add_init_wind(u=u,v=v,lev=pressure,levtype='pressure',levid='lev')

################################################
# 3. Forcing
################################################

# pressure levels

case.add_surface_pressure_forcing(ps,timeid='time')

case.add_pressure_forcing(pressure,timeid='time',lev=pressure,levtype='pressure',levid='lev')

# nudging of wind towards constant values
case.add_wind_nudging(unudg=u,vnudg=v,timescale=3600.,p_nudging=110000.,timeid='time',lev=pressure,levtype='pressure',levid='lev')

# Large scale vertical pressure velocity
omega  = fin['omega'][0,-1:0:-1,0,0]

case.add_vertical_velocity(omega=omega,timeid='time',lev=pressure,levtype='pressure',levid='lev')

# Large-scale advection of temperature in K s-1
tadv  = fin['divT'][0,-1:0:-1,0,0]

case.add_temp_advection(tadv,timeid='time',lev=pressure,levtype='pressure',levid='lev')

# Large-scale advection of specific humidity in kg kg-1 s-1
qvadv=fin['divq'][0,-1:0:-1,0,0]

case.add_qv_advection(qvadv,timeid='time',lev=pressure,levtype='pressure',levid='lev')

# Surface Forcing
# constant surface latent and sensible fluxes [W m-2] (be careful sign convention)
hs=fin['shflx'][0,0,0]
hl=fin['lhflx'][0,0,0]

case.add_surface_fluxes(sens=hs,lat=hl,timeid='time',forc_wind='z0',z0=0.01)

# Constant SST [K]
ts=274.01
case.add_surface_temp(ts,timeid='time')

################################################
# 4. Writing file
################################################

case.write('MPACE_E3SM_DEF_driver.nc')

if lverbose:
    case.info()

################################################
# 5. Ploting, if asked
################################################

if lplot:
    case.plot(rep_images='./images/driver_DEF/',timeunits='hours')
