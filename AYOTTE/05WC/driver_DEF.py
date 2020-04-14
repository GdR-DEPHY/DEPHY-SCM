#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on 10 April 2020

@author: Fleur Couvreux

Modifications:
  14-04-2020, R. Roehrig: some cleaning
"""

## AYOTTE/05WC

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
# This is an idealized case so date and lat/lon are arbitrary

case = Case('AYOTTE/05WC',
        lat=45,
        lon=123.3,
        startDate="20091211100000",
        endDate="20091211170000",
        zorog=0.,
        z0=0.16)

case.set_title("Forcing and initial conditions for AYOTTE-05WC case - Original definition")
case.set_reference("Ayotte et al. (1996, BLM) and files provided by F Hourdin initially from Ayotte")
case.set_author("F. couvreux")
case.set_script("driver_DEF.py")


# time units are expected to be seconds since startDate
t0 = 0     # 10:00 UTC, 11 December 2009
t1 = 25200 # 17:00 UTC, 11 December 2009


################################################
# 2. Initial state
################################################

# Surface pressure
ps = 100000.
case.add_variable('ps',[ps,])

#         z (m) theta (K) rt (g kg-1) u (m s-1) v (m s-1)
init = [   0.0,   300.2,   0.0,       7.5,     0.4,\
         100.0,   300.2,   0.0,      11.8,     0.5,\
         200.0,   300.2,   0.0,      12.5,     0.5,\
         812.0,   300.2,   0.0,      13.5,     0.33,\
         907.0,   300.3,   0.0,      13.71,    0.29,\
         946.0,   300.5,   0.0,      14.04,    0.21,\
        1003.0,   300.8,   0.0,      14.53,    0.1,\
        1155.0,   301.9,   0.0,      15.0,     0.0,\
        1346.0,   303.0,   0.0,      15.0,     0.0,\
        1518.0,   304.0,   0.0,      15.0,     0.0,\
        1671.0,   305.0,   0.0,      15.0,     0.0,\
        1690.0,   305.1,   0.0,      15.0,     0.0,\
        1709.0,   305.2,   0.0,      15.0,     0.0]

init = np.array(init,dtype=np.float64)

z = init[0::5]

case.add_variable('theta',init[1::5],      lev=z,levtype='altitude')
case.add_variable('rt',   init[2::5]/1000.,lev=z,levtype='altitude') # converted in kg kg-1
case.add_variable('u',    init[3::5],      lev=z,levtype='altitude')
case.add_variable('v',    init[4::5],      lev=z,levtype='altitude')

# Turbulent Kinetic Energy
# Initialized to 0 as in the LES
ztke = range(0,3000+1,10)
nztke = len(ztke)
tke = np.zeros(nztke,dtype=np.float64)

case.add_variable('tke',tke,lev=ztke,levtype='altitude')

################################################
# 3. Forcing
################################################

# Constant Geostrophic wind across the simulation
ug = np.zeros((2,13),dtype=np.float64)
ug[0,:] = 15.
ug[1,:] = 15.

vg = np.zeros((2,13),dtype=np.float64)
vg[0,:] = 0.
vg[1,:] = 0.

case.add_variable('ug',ug,time=[t0,t1],lev=z,levtype='altitude')
case.add_variable('vg',vg,time=[t0,t1],lev=z,levtype='altitude')

# Surface Forcing
#            t (s) H (W m-2) LE (W m-2)
sfcForc= [     0,  56.27,    0,\
           25200,  56.27,    0]

sfcForc = np.array(sfcForc,dtype=np.float64)

timeSfc = sfcForc[0::3]

case.add_variable('sfc_sens_flx',sfcForc[1::3],time=timeSfc)
case.add_variable('sfc_lat_flx', sfcForc[2::3],time=timeSfc)

################################################
# 4. Attributes
################################################

# advection of theta and rt
case.set_attribute("adv_theta",0)
case.set_attribute("adv_rt",0)
# no radiation (for now, equivalent to potential temperature radiative tendency included in advection)
case.set_attribute("rad_theta",'adv')
# Geostrophic wind forcing
case.set_attribute("forc_geo",1)
# Surface flux forcing, wind stress is computed using z0
case.set_attribute("surfaceType","land")
case.set_attribute("surfaceForcing","surfaceFlux")
case.set_attribute("surfaceForcingWind","z0")

################################################
# 5. Writing file
################################################

case.write('AYOTTE_05WC_DEF_driver.nc',verbose=False)

if lverbose:
    case.info()

################################################
# 6. Ploting, if asked
################################################

if lplot:
    case.plot(rep_images='./images/driver_DEF/',timeunits='hours')
