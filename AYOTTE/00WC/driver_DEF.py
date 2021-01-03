#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on 10 April 2020

@author: Fleur Couvreux

Modifications:
  2020/04/14, R. Roehrig: some cleaning
  2021/01/03, R. Roehrig: update for improved case definition interface.
"""

## AYOTTE/00WC

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

case = Case('AYOTTE/00WC',
        lat=45,
        lon=123.3,
        startDate="20091211100000",
        endDate="20091211170000",
        surfaceType="land",
        zorog=0.)

case.set_title("Forcing and initial conditions for AYOTTE-00WC case - Original definition")
case.set_reference("Ayotte et al. (1996, BLM) and files provided by F Hourdin initially from Ayotte")
case.set_author("F. couvreux")
case.set_script("DEPHY-SCM/AYOTTE/00WC/driver_DEF.py")

################################################
# 2. Initial state
################################################

# Surface pressure
ps = 100000.
case.add_init_ps(ps)

#         z (m) theta (K) rt (g kg-1) u (m s-1) v (m s-1)
init = [   0.0,   300.3,   0.0,       5.0,     2.4,\
         130.0,   300.3,   0.0,       9.6,     3.9,\
         245.0,   300.3,   0.0,      10.7,     3.9,\
         670.0,   300.3,   0.0,      12.5,     3.6,\
         900.0,   300.3,   0.0,      13.5,     2.9,\
         904.0,   300.3,   0.0,      13.53,    2.85,\
        1000.0,   300.58,  0.0,      14.3,     1.7,\
        1010.0,   300.6,   0.0,      14.37,    1.53,\
        1100.0,   301.45,  0.0,      15.0,     0.0,\
        1116.0,   301.6,   0.0,      15.0,     0.0,\
        1138.0,   301.8,   0.0,      15.0,     0.0,\
        1159.0,   301.9,   0.0,      15.0,     0.0,\
        1350.0,   303.0,   0.0,      15.0,     0.0,\
        1670.0,   305.0,   0.0,      15.0,     0.0,\
        1750.0,   305.48,  0.0,      15.0,     0.0,\
        1903.0,   306.4,   0.0,      15.0,     0.0,\
        2400.0,   309.39,  0.0,      15.0,     0.0]

init = np.array(init,dtype=np.float64)

z = init[0::5]

case.add_init_theta(init[1::5], lev=z, levtype='altitude')
case.add_init_rt(init[2::5]/1000., lev=z, levtype='altitude') # converted in kg kg-1
case.add_init_wind(u=init[3::5],v=init[4::5], lev=z, levtype='altitude')

# Turbulent Kinetic Energy
# Initialized to 0 as in the LES
ztke = range(0,3000+1,10)
nztke = len(ztke)
tke = np.zeros(nztke,dtype=np.float64)

case.add_init_tke(tke, lev=ztke, levtype='altitude')

################################################
# 3. Forcing
################################################

# Constant Geostrophic wind across the simulation
nz = len(z)
ug = np.zeros(nz,dtype=np.float64) + 15.
vg = np.zeros(nz,dtype=np.float64) + 0.

case.add_geostrophic_wind(ug=ug,vg=vg,lev=z,levtype='altitude')

# No radiation
case.deactivate_radiation()

# Surface Forcing
#            t (s) H (W m-2) LE (W m-2)
sfcForc= [     0,  0.0,    0,\
           25200,  0.0,    0]

timeSfc = sfcForc[0::3]

case.add_surface_fluxes(sens=sfcForc[1::3],lat=sfcForc[2::3],time=timeSfc,forc_wind='z0',z0=0.16)

################################################
# 4. Writing file
################################################

case.write('AYOTTE_00WC_DEF_driver.nc',verbose=False)

if lverbose:
    case.info()

################################################
# 5. Ploting, if asked
################################################

if lplot:
    case.plot(rep_images='./images/driver_DEF/',timeunits='hours')
