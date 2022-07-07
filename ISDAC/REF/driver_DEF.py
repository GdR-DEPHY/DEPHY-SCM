#!/usr/bin/env python

# -*- coding: utf-8 -*-
"""
Created on 25 January 2021

@author: Etienne Vignon

Modifications:
   2022/07/07, Romain Roehrig, update and bugfixes for nudging coefficients

ISDAC original case definition
From Ovchinnikov, M., et al. (2014), JAMES. Initial profiles of temperature, moisture, and horizontal wind components are based on aircraft observations
in the mixed layer and idealization of a sounding at Barrow, AK. For simplicity, both sensible and latent surface heat fluxes are set to zero in the
simulations
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

case = Case('ISDAC/REF',
        lat=71.32,
        lon=-156.61,
        startDate="20080426180000",
        endDate="20080427020000",
        surfaceType="seaice",
        zorog=0.)

case.set_title("Forcing and initial conditions for ISDAC case - Original definition")
case.set_reference("Ovchinnikov, M., et al. (2014), JAMES")
case.set_author("E. Vignon")
case.set_script("DEPHY-SCM/ISDAC/REF/driver_DEF.py")
case.set_comment(" ")

################################################
# 2. Initial state
################################################
# Surface pressure
ps = 102000.

z=np.array(np.linspace(0,5000,51),dtype=np.float64)
thetal=z*0.
qt=z*0
u=z*0
v=z*0.
tke=z*0

thetal[z<400]=265.+0.004*(z[z<400]-400)
thetal[z>=400]=265.
thetal[z>=825]=266+(z[z>=825]-825.)**0.3
thetal[z>=2045]=271+(z[z>=2045]-2000.)**0.33

qt[z<400]=1.5-0.00075*(z[z<400]-400.)
qt[z>=400]=1.5
qt[z>=825]=1.2
qt[z>=2045]=0.5-0.000075*(z[z>=2045]-2045)

qt=qt/1000.0 # to kg/kg


u=u-7.
v=-2.+0.003*z
tke = tke+0.1

# Surface pressure
case.add_init_ps(ps)

# height
case.add_init_height(z,lev=z,levtype='altitude')

# Liquid potential Temp
case.add_init_thetal(thetal,lev=z,levtype='altitude')

# Total water 
case.add_init_qt(qt,lev=z,levtype='altitude')

# Wind
case.add_init_wind(u=u,v=v,lev=z,levtype='altitude')

# TKE
case.add_init_tke(tke,lev=z,levtype='altitude')



################################################
# 3. Forcing
################################################


# Surface Forcing. Constant  surface temperature
ts = 267.0

# 0 latent and sensible heat flux at the surface
hs=0.
hl=0.

# nudging of wind towards constant values
u_nudg=u
v_nudg=v


zuv=825.
inv_tau_nudg_uv=u*0.
index = (z>0) & (z<=zuv)
inv_tau_nudg_uv[index]=1./7200.*(1.-np.cos(np.pi*(z[index]/zuv)))/2.
inv_tau_nudg_uv[z>zuv]=1./7200.

# nudging of thetal
z1=1200.
z2=1500.

thetal_nudg=thetal
inv_tau_nudg_thetal=thetal*0.
inv_tau_nudg_thetal[z>=z1]=1./3600.*(1.-np.cos(np.pi*(z[z>=z1]-z1)/(z2-z1)))/2.
inv_tau_nudg_thetal[z>z2]=1/3600.

# nudging of qt
qt_nudg=qt
inv_tau_nudg_qt=inv_tau_nudg_thetal

# Large scale vertical pressure velocity

w=-5.0e-6*z
w[z>=825]=-0.4125e-2

# pressure levels

case.add_surface_pressure_forcing(ps,timeid='time')

# nudging of wind. Set timescale to -1 so that we read vertical profiles of 
# nudging time scale tau_nudging [s]
case.add_wind_nudging(unudg=u,vnudg=v,timeid='time',lev=z,levtype='altitude',
                      nudging_coefficient=inv_tau_nudg_uv,lev_coef=z)

# Large scale vertical  velocity
case.add_vertical_velocity(w=w,timeid='time',lev=z,levtype='altitude')

# thetalnudging. Set timescale to -1 so that we read vertical profiles of 
# nudging time scale tau_nudging  [s]
case.add_thetal_nudging(thetal_nudg,timeid='time',lev=z,levtype='altitude',
                        nudging_coefficient=inv_tau_nudg_thetal,lev_coef=z)

# qt nudging. Set timescale to -1 so that we read vertical profiles of 
# nudging time scale tau_nudging [s]
case.add_qt_nudging(qt_nudg,timeid='time',lev=z,levtype='altitude',
                    nudging_coefficient=inv_tau_nudg_qt,lev_coef=z)


# Surface Forcing
# constant surface latent and sensible fluxes [W m-2] (be careful sign convention)
case.add_surface_fluxes(sens=hs,lat=hl,timeid='time',forc_wind='z0',z0=4e-4)

# Constant SST [K]
case.add_surface_temp(ts,timeid='time')



################################################
# 4. Writing file
################################################

case.write('ISDAC_REF_DEF_driver.nc')

if lverbose:
    case.info()

################################################
# 5. Ploting, if asked
################################################

if lplot:
    case.plot(rep_images='./images/driver_DEF/',timeunits='hours')
