#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on 01 October 2020

@author: Etienne Vignon

Modifications:


MPACE original case definition
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
lverbose = True # print information about variables and case

################################################
# 1. General information about the case
################################################
# This is an idealized case so date is arbritrary
# lat/lon are fixed to a lat representative of Arctic conditions
# 8h duration with a constant surface cooling

case = Case('MPACE/REF',
        lat=71.75,
        lon=209.00,
        startDate="20041009170000",
        endDate="20041010050000",
        zorog=0.,
        z0= 0.01)

case.set_title("Forcing and initial conditions for MPACE case - Original definition")
case.set_reference("Klein et al. (2009, QJRMS)")
case.set_author("E. Vignon")
case.set_script("driver_DEF.py")
case.set_comment("Use of forcing file from E3SM, https://github.com/E3SM-Project/scmlib/wiki/E3SM-Single-Column-Model-Case-Library")

# time units are expected to be seconds since startDate
t0 = 0     # 17:00 UTC, 09 October 2004
t1 = 43200 # 05:00 UTC, 10 October 2004

################################################
# 2. Input netCDF file
################################################

fin = nc.Dataset('MPACE_forcing_dephy.nc','r')


################################################
# 3. Initial state
################################################

# Surface pressure
ps  =  fin['Ps'][0,0,0]
case.add_variable('ps',ps,timeid='time')


# Pressure
pressure  = fin['lev'][-1:0:-1]
pressure[0] = ps
#pressure = extend('p',pressure,init=True)

case.add_variable('pressure',pressure,lev=pressure,levtype='pressure',levid='lev',timeid='time')
#case.add_variable('pressure_forc',pressure,lev=pressure,levtype='pressure',levid='lev',timeid='time')

# Temperature in K and moisture (vapor and total mixing ratio)
# Note that in the ref paper, the ice-liquid potential temperature is initially set 
# but the initial state corresponds to a purely-liquid cloud

temp  = fin['T'][0,-1:0:-1,0,0]
case.add_variable('temp',temp,lev=pressure,levtype='pressure',levid='lev',timeid='time')


#  water vapor mixing ratio in kg/kg
rv =  fin['q'][0,-1:0:-1,0,0]
qv=rv/(1.+rv)
case.add_variable('qv',qv,lev=pressure,levtype='pressure',levid='lev',timeid='time')



#  total water mixing ratio in kg/kg
rl =  fin['CLDLIQ'][0,-1:0:-1,0,0]
ri =  fin['CLDICE'][0,-1:0:-1,0,0]
rt=rl+ri+rv
qt=rt/(1.+rt)
case.add_variable('qt',qt,lev=pressure,levtype='pressure',levid='lev',timeid='time')


# Zonal wind
u  = fin['u'][0,-1:0:-1,0,0]
case.add_variable('u',u,lev=pressure,levtype='pressure',levid='lev',timeid='time')

# Meridional wind
v  = fin['v'][0,-1:0:-1,0,0]
case.add_variable('v',v,lev=pressure,levtype='pressure',levid='lev',timeid='time')


################################################
# 3. Forcing
################################################

# pressure levels

ps_forc = np.zeros((2,1),dtype=np.float)
ps_forc[:]=ps
case.add_variable('ps_forc',ps_forc,time=[t0,t1],timeid='time')
pressure_forc = np.zeros((2,len(pressure)),dtype=np.float)
pressure_forc[0,:]=pressure
pressure_forc[1,:]=pressure
case.add_variable('pressure_forc',pressure_forc,time=[t0,t1],lev=pressure,levtype='pressure',levid='lev',timeid='time')

# nudging of wind towards constant values
u_nudg=np.zeros((2,len(u)),dtype=np.float)
v_nudg=np.zeros((2,len(v)),dtype=np.float)
u_nudg[0,:]=u
u_nudg[1,:]=u
v_nudg[0,:]=v
v_nudg[1,:]=v
case.add_variable('u_nudging',u_nudg,time=[t0,t1],lev=pressure,levtype='pressure',levid='lev',timeid='time')
case.add_variable('v_nudging',v_nudg,time=[t0,t1],lev=pressure,levtype='pressure',levid='lev',timeid='time')



# Large scale vertical pressure velocity

omegar  = fin['omega'][0,-1:0:-1,0,0]
omega=np.zeros((2,len(omegar)),dtype=np.float)
omega[0,:]=omegar
omega[1,:]=omegar
case.add_variable('omega',omega,time=[t0,t1],lev=pressure,levtype='pressure',levid='lev',timeid='time')

# Large-scale advection of temperature in K s-1

tadvr  = fin['divT'][0,-1:0:-1,0,0]
tadv=np.zeros((2,len(tadvr)),dtype=np.float)
tadv[0,:]=tadvr
tadv[1,:]=tadvr
case.add_variable('temp_adv',tadv,time=[t0,t1],lev=pressure,levtype='pressure',levid='lev',timeid='time') # converted in K s-1 (array type required)

# Large-scale advection of specific humidity in kg kg-1 s-1

qvadvr=fin['divq'][0,-1:0:-1,0,0]
qvadv=np.zeros((2,len(qvadvr)),dtype=np.float)
qvadv[0,:]=qvadvr
qvadv[1,:]=qvadvr
case.add_variable('qv_adv',qvadv,time=[t0,t1],lev=pressure,levtype='pressure',levid='lev',timeid='time') # converted in kg kg-1 s-1 (array type required)


# Surface Forcing
# constant surface latent and sensible fluxes [W m-2] (be careful sign convention)

hs=fin['shflx'][0,0,0]
hl=fin['lhflx'][0,0,0]

case.add_variable('sfc_sens_flx',[hs,hs],time=[t0,t1],timeid='time')
case.add_variable('sfc_lat_flx', [hl,hl],time=[t0,t1],timeid='time')

# Constant SST [K]
ts=274.01
case.add_variable('ts',[ts,ts],time=[t0,t1],timeid='time')

################################################
# 4. Attributes
################################################

# advection of theta and rt
case.set_attribute("adv_temp",1)
case.set_attribute("adv_qv",1)
# No Geostrophic wind forcing
case.set_attribute("forc_geo",0)
# Vertical pressure velocity
case.set_attribute("forc_omega",1)
# Nudged wind
case.set_attribute("nudging_u",3600.)
case.set_attribute("nudging_v",3600.)
case.set_attribute("p_nudging_u",110000.)
case.set_attribute("p_nudging_v",110000.)
# Surface type is ocean
case.set_attribute("surfaceType","ocean")
# Surface flux forcing, wind stress is computed using z0
case.set_attribute("surfaceForcingWind","z0")
#case.set_attribute("surfaceForcing","ts")
case.set_attribute("surfaceForcing","surfaceFlux")
################################################
# 5. Writing file
################################################

case.write('MPACE_REF_DEF_driver.nc',verbose=False)

if lverbose:
    case.info()

################################################
# 6. Ploting, if asked
################################################

if lplot:
    case.plot(rep_images='./images/driver_DEF/',timeunits='hours')
