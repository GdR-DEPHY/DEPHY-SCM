#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on 7 April 2020

@author: Romain Roehrig
"""

## DYNAMO/NSA3A_D1 case definition

import sys
sys.path = ['../','../../utils/',] + sys.path

import netCDF4 as nc
import numpy as np

from Case import Case

import thermo
from myfunctions import extend

################################################
# 0. General configuration of the present script
################################################

lplot = True # plot all the variables
lverbose = False # print information about variables and case

################################################
# 1. General information about the case
################################################

case = Case('DYNAMO/NSA3A_D1',
        lat=2.5,
        lon=77.5,
        startDate="20111015000000",
        endDate=  "20111105000000",
        surfaceType="ocean",
        zorog=0.)

case.set_title("Forcing and initial conditions for DYNAMO/NSA3A case - MJO event #1")
case.set_reference("Abdel-Lathif et al. (2018, JAMES); Johnson et al. (2015, JAS); Input file dynamo_nsa_v3a.nc downloaded at http://johnson.atmos.colostate.edu/dynamo/products/array_averages/")
case.set_author("R. Roehrig")
case.set_script("DEPHY-SCM/DYNAMO/NSA3A_MJO1/driver_DEF.py")
case.set_comment("Data from CSU (http://johnson.atmos.colostate.edu/dynamo/products/array_averages/) complemented with ERA-Interim above 50 hPa and the standard atmosphere above 1 hPa. Specific humidity is computed from CSU relative humidity below 50 hPa.")

################################################
# 2. Input netCDF file
################################################

fin = nc.Dataset('../dynamo_nsa_v3a.nc')

################################################
# 2. Initial state
################################################

t0 = 14*8 # 00:00 UTC, 15 October 2011

# Surface pressure
ps = fin['ps'][t0]*100. # from hPa to Pa
case.add_init_ps(ps)

# Pressure
pressure = fin['level'][:]*100. # From hPa to Pa
# Modify first surface level value (1025 hPa is used for convenience in CSU file)
pressure[0] = ps
pressure = extend('p',pressure,init=True)

case.add_init_pressure(pressure,lev=pressure,levtype='pressure',levid='lev')

# Height
height  = fin['z'][t0,:]
height = extend('zg',height,init=True)

case.add_init_height(height,lev=pressure,levtype='pressure',levid='lev')

# Zonal and meridional wind
u  = fin['u'][t0,:]
u = extend('ua',u,init=True)

v  = fin['v'][t0,:]
v = extend('va',v,init=True)

case.add_init_wind(u=u,v=v,lev=pressure,levtype='pressure',levid='lev')

# Temperature
temp = fin['T'][t0,:] + 273.15 # from °C to K
temp = extend('ta',temp,init=True)

case.add_init_temp(temp,lev=pressure,levtype='pressure',levid='lev')

# Specific humidity
#rv = fin['wmr'][0,:]/1000. # from g kg-1 to kg kg-1
#qv = rv/(1.+rv) # from rv to qv

# Specific humidity has not the appropriate precision in original CSU ASCII file, as well as in the derived netcdf files
# It leads to zero humidity in the upper troposphere, which is unrealistic
# Therefore specific humidity is computed from relative humidity, temperature and pressure from the surface to 50 hPa

rh = fin['rh'][t0,:]
t = fin['T'][t0,:] + 273.15 # from °C to K
pres = fin['level'][:]*100. # from hPa to Pa
pres[0] = ps

qv = thermo.rh2qv(rh,t,pres)
qv = extend('hus',qv,init=True)

case.add_init_qv(qv,lev=pressure,levtype='pressure',levid='lev')

################################################
# 3. Forcing
################################################

t0 = 0         # 00:00 UTC, 15 October 2011
t1 = 20*24.*3600 # 00:00 UTC, 05 November 2011 - (20*3)-hour long simulation

timeForc = np.arange(t0,t1,3600.*3.) # Forcing data are provided every 3 hours
levForc = fin['level'][:]*100. # first level pressure (1025 hPa) is arbitrary and corresponds to surface values
levForc = extend('p',levForc,init=True)

nt, = timeForc.shape
nlev, = levForc.shape

nt0 = 14*8
nt1 = nt0 + 20*8

# Surface pressure
ps_forc = fin['ps'][nt0:nt1]*100. # from hPa to Pa
case.add_surface_pressure_forcing(ps_forc,time=timeForc,timeid='time')

# Forcing pressure level
pressure_forc = np.tile(fin['level'][:]*100.,(nt,1))
pressure_forc[:,0] = ps_forc[:]
pressure_forc = extend('p',pressure_forc,init=False)
case.add_pressure_forcing(pressure_forc,time=timeForc,timeid='time',lev=levForc,levtype='pressure',levid='lev')

# Forcing height
height_forc = fin['z'][nt0:nt1,:]
height_forc = extend('zg',height_forc,init=False,time=timeForc)
case.add_height_forcing(height_forc,time=timeForc,timeid='time',lev=levForc,levtype='pressure',levid='lev')

# Temperature horizontal advection
hT = fin['hT'][nt0:nt1,:]*-1. # from advective tendency to advective forcing
hT = extend('advT',hT,init=False,time=timeForc)
case.add_temp_advection(hT,time=timeForc,timeid='time',lev=levForc,levtype='pressure',levid='lev')

# Specific humidity horizontal advection
hq = fin['hq'][nt0:nt1,:]/1000.*-1. # from advective tendency to advective forcing and from g kg-1 s-1 to kg kg-1 s-1
hq = extend('advq',hq,init=False,time=timeForc)
case.add_qv_advection(hq,time=timeForc,timeid='time',lev=levForc,levtype='pressure',levid='lev')

# Large-scale velocity
omega = fin['omega'][nt0:nt1,:]/3600.*100. # from hPa hour-1 to Pa s-1
omega = extend('wap',omega,init=False,time=timeForc)
case.add_vertical_velocity(omega=omega,time=timeForc,timeid='time',lev=levForc,levtype='pressure',levid='lev')

# Wind nudging
unudg = fin['u'][nt0:nt1,:]
vnudg = fin['v'][nt0:nt1,:]

unudg = extend('ua',unudg,init=False,time=timeForc)
vnudg = extend('va',vnudg,init=False,time=timeForc)

case.add_wind_nudging(unudg=unudg,vnudg=vnudg,timescale=3600.*3.,p_nudging=110000.,time=timeForc,timeid='time',lev=levForc,levtype='pressure',levid='lev')

# Nudging of temperature and humidity at high levels
tnudg = fin['T'][nt0:nt1,:] + 273.15 # from °C to K
rvnudg = fin['wmr'][nt0:nt1,:]/1000. # from g kg-1 s-1 to kg kg-1 s-1
qvnudg = rvnudg/(1.+rvnudg)

tnudg = extend('ta',tnudg,init=False,time=timeForc)
qvnudg = extend('hus',qvnudg,init=False,time=timeForc)

case.add_temp_nudging(tnudg,timescale=3600.*3.,p_nudging=5000.,time=timeForc,timeid='time',lev=levForc,levtype='pressure',levid='lev')
case.add_qv_nudging(qvnudg,timescale=3600.*3.,p_nudging=5000.,time=timeForc,timeid='time',lev=levForc,levtype='pressure',levid='lev')

# Sea surface temperature
ts = fin['T'][nt0:nt1,0] + 273.15 # from °C to K
case.add_forcing_ts(ts,time=timeForc,timeid='time')

fin.close()

################################################
# 4. Writing file
################################################

case.write('DYNAMO_NSA3A_MJO1_DEF_driver.nc')

if lverbose:
    case.info()

################################################
# 5. Ploting, if asked
################################################

if lplot:
    case.plot(rep_images='./images/driver_DEF/',timeunits='hours',levunits='hPa')
