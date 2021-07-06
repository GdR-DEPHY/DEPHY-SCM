#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on 09 January 2021

@author: Romain Roehrig
"""

## ARM-Cumulus case definition from E3SM database

import os
import sys
sys.path = ['../../utils/',] + sys.path

import datetime

import netCDF4 as nc
import numpy as np

import constants as CC

from Case import Case

################################################
# 0. General configuration of the present script
################################################

lplot = True # plot all the variables
lverbose = False # print information about variables and case

case_name = 'ARMCU'
subcase_name = 'E3SM'
title = "Forcing and initial conditions for ARM-Cumulus case - ES3M database, provided by P. Bogenschutz"
reference = "Brown et al. (2002, QJRMS)"
author = " R. Roehrig"

fin = nc.Dataset('ARM_shallow_iopfile.nc')

################################################
# 1. General information about the case
################################################

orog = float(fin['phis'][0,0])/CC.g
# to be consistent with other flavors of ARMCU
orog = 314.

lat = float(fin['lat'][0])
lon = float(fin['lon'][0])

bdate = str(fin['bdate'][0])
tsec = fin['tsec'][:]
seconds0 = tsec[0]
seconds1 = tsec[-1]

swoff = False
lwoff = False

# Update from running script
with open('run_e3sm_scm_ARM_shallow.csh', 'r') as f:
    for line in f.readlines():
        if 'set lat' in line:
            lat = float(line.split()[3])
        if 'set lon' in line:
            lon = float(line.split()[3])
        if 'set startdate' in line:
            date = line.split()[3]
            bdate = ''.join(date.split('-'))
        if 'set start_in_sec' in line:
            seconds0 = float(line.split()[3])
        if 'set stop_option' in line:
            stop_option = line.split()[3]
        if 'set stop_n' in line:
            ntime = float(line.split()[3])
            if stop_option == 'nhours':
                seconds1 = ntime*3600
            else:
                raise ValueError
        if 'set do_turnoff_swrad' in line:
            swoff = line.split()[3] == '.false'
        if 'set do_turnoff_lwrad' in line:
            lwoff = line.split()[3] == '.false'
           
h0 = int(seconds0/3600)
m0 = int((seconds0-h0*3600)/60)
s0 = int(seconds0-h0*3600-m0*60)

d0 = datetime.datetime(int(bdate[0:4]),int(bdate[4:6]),int(bdate[6:8]))
t0 = datetime.timedelta(seconds=seconds0)

t1 = datetime.timedelta(seconds=int(seconds1/3600)*3600)

d1 = d0+t1

startDate = d0.strftime('%Y%m%d%H%M%S')
endDate = d1.strftime('%Y%m%d%H%M%S')

# to be consistent with other flavors of ARMCU
startDate="19970621113000"
endDate="19970622020000"

case = Case('{0}/{1}'.format(case_name,subcase_name),
        lat=lat,
        lon=lon,
        startDate=startDate,
        endDate=endDate,
        surfaceType='land',
        zorog=orog)

case.set_title(title)
case.set_reference(reference)
case.set_author(author)
case.set_script("DEPHY-SCM/{0}/{1}/driver_DEF.py".format(case_name,subcase_name))

################################################
# 2. Initial state
################################################

# Surface pressure
ps = fin['Ps'][0,0,0]
case.add_init_ps(ps)

# Pressure
plev = fin['lev'][:]
case.add_init_pressure(plev,lev=plev,levtype='pressure')

# Zonal and meridional wind
u = fin['u'][0,:,0,0]
v = fin['v'][0,:,0,0]
case.add_init_wind(u=u,v=v,lev=plev,levtype='pressure')

# Temperature
temp = fin['T'][0,:,0,0]
case.add_init_temp(temp,lev=plev,levtype='pressure')

# Specific humidity
qv = fin['q'][0,:,0,0]
case.add_init_qv(qv,lev=plev,levtype='pressure')

################################################
# 3. Forcing
################################################

# number of forcing time step
nforc, = tsec.shape

# surface pressure forcing
ps_forc = fin['Ps'][:,0,0] 
case.add_surface_pressure_forcing(ps_forc,time=tsec)

# Pressure forcing
pressure_forc = np.tile(plev,(nforc,1))
case.add_pressure_forcing(pressure_forc,time=tsec,lev=plev,levtype='pressure')

# Geostrophic wind
ug = fin['u'][:,:,0,0]
vg = fin['v'][:,:,0,0]
case.add_geostrophic_wind(ug=ug,vg=vg,time=tsec,lev=plev,levtype='pressure')

# Advections
temp_adv = fin['divT'][:,:,0,0]
qv_adv = fin['divq'][:,:,0,0]
case.add_temp_advection(temp_adv,time=tsec,lev=plev,levtype='pressure',include_rad=not(swoff or lwoff))
case.add_qv_advection(qv_adv,time=tsec,lev=plev,levtype='pressure')

# Vertical velocity
omega = fin['omega'][:,:,0,0]
case.add_vertical_velocity(omega=omega,time=tsec,lev=plev,levtype='pressure')

# Surface forcing
ts = fin['Tg'][:,0,0]
shf = fin['shflx'][:,0,0]
lhf = fin['lhflx'][:,0,0]
case.add_surface_fluxes(sens=shf,lat=lhf,time=tsec,forc_wind='z0',z0=0.035)

case.add_surface_temp(ts,time=tsec)

################################################
# 4. Writing file
################################################

case.write('{0}_{1}_DEF_driver.nc'.format(case_name,subcase_name))

if lverbose:
    case.info()

################################################
# 5. Ploting, if asked
################################################

if lplot:
    case.plot(rep_images='./images/driver_DEF/',timeunits='hours',levunits='hPa')

fin.close()
