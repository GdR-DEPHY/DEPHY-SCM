#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on 17 June 2020

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

lat = float(fin['lat'][0])
lon = float(fin['lon'][0])
orog = float(fin['phis'][0,0])/CC.g

bdate = str(fin['bdate'][0])
tsec = fin['tsec'][:]

h0 = int(tsec[0]/3600)
m0 = int((tsec[0]-h0*3600)/60)
s0 = int(tsec[0]-h0*3600-m0*60)

d0 = datetime.datetime(int(bdate[0:4]),int(bdate[4:6]),int(bdate[6:8]))
t0 = datetime.timedelta(seconds=int(tsec[0]))

t1 = datetime.timedelta(seconds=int(tsec[-1]/3600)*3600)

d1 = d0+t1

startDate = d0.strftime('%Y%m%d%H%M%S')
endDate = d1.strftime('%Y%m%d%H%M%S')

case = Case('{0}/{1}'.format(case_name,subcase_name),
        lat=lat,
        lon=lon,
        startDate=startDate,
        endDate=endDate,
        zorog=orog,
        z0=0.035)

case.set_title(title)
case.set_reference(reference)
case.set_author(author)
case.set_script("DEPHY-SCM/{0}/{1}/driver_DEF.py".format(case_name,subcase_name))

################################################
# 2. Initial state
################################################

# Surface pressure
ps = fin['Ps'][0,0,0]
case.add_variable('ps',[ps,])

plev = fin['lev'][:]
case.add_variable('pressure',plev,      lev=plev,levtype='pressure')

case.add_variable('temp',fin['T'][0,:,0,0],      lev=plev,levtype='pressure')
case.add_variable('qv',  fin['q'][0,:,0,0],      lev=plev,levtype='pressure')
case.add_variable('u',   fin['u'][0,:,0,0],      lev=plev,levtype='pressure')
case.add_variable('v',   fin['v'][0,:,0,0],      lev=plev,levtype='pressure')

################################################
# 3. Forcing
################################################

# number of forcing time step
nforc, = tsec.shape

pressure_forc = np.tile(plev,(nforc,1))
case.add_variable('pressure_forc',pressure_forc,time=tsec,lev=plev,levtype='pressure')

# Constant Geostrophic wind across the simulation
case.add_variable('ug',fin['u'][:,:,0,0],time=tsec,lev=plev,levtype='pressure')
case.add_variable('vg',fin['v'][:,:,0,0],time=tsec,lev=plev,levtype='pressure')

case.add_variable('ts',fin['Tg'][:,0,0]*0+310,time=tsec)
case.add_variable('sfc_sens_flx',fin['lhflx'][:,0,0],time=tsec)
case.add_variable('sfc_lat_flx', fin['shflx'][:,0,0],time=tsec)

case.add_variable('temp_adv',fin['divT'][:,:,0,0],time=tsec,lev=plev,levtype='pressure') 
case.add_variable('qv_adv',  fin['divq'][:,:,0,0],time=tsec,lev=plev,levtype='pressure') 


################################################
# 4. Attributes
################################################

# advection of theta and rt
case.set_attribute("adv_temp",1)
case.set_attribute("adv_qv",1)
# potential temperature radiative tendency is included in advection
case.set_attribute("rad_temp","adv")
# Geostrophic wind forcing
case.set_attribute("forc_geo",1)
# Surface flux forcing, wind stress is computed using z0
case.set_attribute("surfaceType","land")
case.set_attribute("surfaceForcing","surfaceFlux")
case.set_attribute("surfaceForcingWind","z0")

################################################
# 5. Writing file
################################################

case.write('{0}_{1}_DEF_driver.nc'.format(case_name,subcase_name),verbose=False)

if lverbose:
    case.info()

################################################
# 6. Ploting, if asked
################################################

if lplot:
    case.plot(rep_images='./images/driver_DEF/',timeunits='hours',levunits='hPa')

fin.close()
