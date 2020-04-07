#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on 01 April 2020

@author: Fleur Couvreux

Modifications:
  07-04-2020, R. Roehrig: correction to deactivate radiation and for surface flux time axis
"""

## ARM-Cumulus original case definition
## From http://projects.knmi.nl/eurocs/ARM/case_ARM_html/

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

case = Case('AYOTTE/24SC',
        lat=45,
        lon=123.3,
        startDate="20091211100000",
        endDate="20091211170000",
        zorog=0.,
        z0=0.16)

case.set_title("Forcing and initial conditions for AYOTTE-24SC case - Original definition")
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
init = [   0.0,   301.1,   0.0,       8.0,     0.4,\
         130.0,   301.1,   0.0,      12.0,     0.6,\
         829.0,   301.1,   0.0,      12.0,     0.6,\
         848.0,   301.2,   0.0,      12.0,     0.6,\
         900.0,   301.29,  0.0,      12.0,     0.6,\
         908.0,   301.3,   0.0,      12.1,     0.576,\
         928.0,   301.4,   0.0,      12.34,    0.516,\
         968.0,   301.8,   0.0,      12.82,    0.396,\
        1000.0,   303.16,  0.0,      13.2,     0.3,\
        1008.0,   303.5,   0.0,      13.34,    0.276,\
        1048.0,   308.2,   0.0,      14.06,    0.156,\
        1100.0,   308.32,  0.0,      15.0,     0.0,\
        1388.0,   309.00,  0.0,      15.0,     0.0,\
        1750.0,   310.09,  0.0,      15.0,     0.0,\
        1787.0,   310.20,  0.0,      15.0,     0.0,\
        2000.0,   310.84,  0.0,      15.0,     0.0,\
        3000.0,   313.85,  0.0,      15.0,     0.0]

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
ug = np.zeros((2,17),dtype=np.float64)
ug[0,:] = 15.
ug[1,:] = 15.

vg = np.zeros((2,17),dtype=np.float64)
vg[0,:] = 0.
vg[1,:] = 0.


case.add_variable('ug',ug,time=[t0,t1],lev=z,levtype='altitude')
case.add_variable('vg',vg,time=[t0,t1],lev=z,levtype='altitude')

# Surface Forcing
#            t (s) H (W m-2) LE (W m-2)
sfcForc= [     0,  270.096,    0,\
           25200,  270.096,    0]

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
# potential temperature radiative tendency is included in advection
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

case.write('AYOTTE_24SC_DEF_driver.nc',verbose=False)

if lverbose:
    case.info()

################################################
# 6. Ploting, if asked
################################################

if lplot:
    case.plot(rep_images='./images/driver_DEF/',timeunits='hours')
