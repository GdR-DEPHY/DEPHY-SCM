#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 6 September 2023

@author: Romain Roehrig

Modification
"""

## ARPEGE/SODANKYLA_2018031512 case definition
## No large-scale forcing, no radiation

CASE = 'ARPEGE'
SUBCASE = 'SODANKYLA_2018031512_NOFORC_NORAD'

import sys
sys.path.append('..')

import netCDF4 as nc
import numpy as np

from dephycf.Case import Case
from dephycf import thermo
from dephycf import constants as CC

################################################
# 0. General configuration of the present script
################################################

lplot = True # plot all the variables
lverbose = False # print information about variables and case

################################################
# 1. General information about the case
################################################

case = Case(f'{CASE}/{SUBCASE}',
        startDate="20180315120000",
        endDate=  "20180318180000",
        surfaceType="land")

case.set_title(f"Forcing and initial conditions for {CASE}/{SUBCASE} case")
case.set_reference("")
case.set_author("R. Roehrig")
case.set_script(f"DEPHY-SCM/{CASE}/{SUBCASE}/driver_DEF.py")
case.set_comment("Case implemented from an ARPEGE forecast. No large-scale forcing and no radiation")

################################################
# 2. Input netCDF file
################################################

fin = nc.Dataset('Arpege-oper-L105_Sodankyla_2018031512.nc')

################################################
# 2. Initial state
################################################

# Orography
orog = fin['height_h'][0,-1]
case.add_orography(orog)

# Surface pressure
ps = fin['ps'][0]
case.add_init_ps(ps)

# Pressure
pressure = fin['pressure_f'][0,:]
case.add_init_pressure(pressure,lev=pressure,levtype='pressure',levid='lev')

# Height
height  = fin['height_f'][0,:]
case.add_init_height(height,lev=pressure,levtype='pressure',levid='lev')

# Zonal and meridional wind
u  = fin['u'][0,:]
v  = fin['v'][0,:]
case.add_init_wind(u=u,v=v,lev=pressure,levtype='pressure',levid='lev')

# Temperature
temp = fin['t'][0,:]
case.add_init_temp(temp,lev=pressure,levtype='pressure',levid='lev')

# Specific humidity
qv = fin['qv'][0,:]
case.add_init_qv(qv,lev=pressure,levtype='pressure',levid='lev')

# Liquid water
ql = fin['ql'][0,:]
case.add_init_variable('ql',ql,lev=pressure,levtype='pressure',levid='lev')

# Ice water
qi = fin['qi'][0,:]
case.add_init_variable('qi',qi,lev=pressure,levtype='pressure',levid='lev')

################################################
# 3. Forcing
################################################

t0 = 0       # 12:00 UTC, 15 March 2018
t1 = 79*3600 # 18:00 UTC, 18 March 2018 - (24*3+6+1)-hour long simulation

timeForc = np.arange(t0,t1,3600.*1.) # Forcing data are provided every hour
levForc = fin['pressure_f'][0,:]

nt, = timeForc.shape
nlev, = levForc.shape

# latitude/longitude
lat = fin['lat'][:nt]
lon = fin['lon'][:nt]
case.add_latitude(lat,time=timeForc,timeid='time')
case.add_longitude(lon,time=timeForc,timeid='time')

# Surface pressure
ps_forc = fin['ps'][:nt]*0 + fin['ps'][0] # Constant surface pressure
case.add_surface_pressure_forcing(ps_forc,time=timeForc,timeid='time')

# Surface forcing
sens = fin['sfc_sens_flx'][:nt]*-1
lat = fin['sfc_lat_flx'][:nt]*-1
z0 = fin['mom_rough'][:nt]
case.add_surface_fluxes(sens=sens,lat=lat,time=timeForc[:-1],timeid='timef',forc_wind='z0',z0=z0,time_z0=timeForc[:-1])

# Surface temperature
ts = fin['t_soil'][:nt,0]
case.add_surface_skin_temp(ts,time=timeForc,timeid='time')

# Activate radiation
case.deactivate_radiation()

fin.close()

################################################
# 4. Writing file
################################################

case.write(f'{CASE}_{SUBCASE}_DEF_driver.nc')

if lverbose:
    case.info()

################################################
# 5. Ploting, if asked
################################################

if lplot:
    case.plot(rep_images='./images/driver_DEF/',timeunits='days',levunits='hPa')
