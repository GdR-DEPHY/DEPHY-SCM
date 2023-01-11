#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 01 April 2020

@author: Fleur Couvreux

Modifications:
  2020/04/07, R. Roehrig: some cleaning
  2021/01/03, R. Roehrig: update for improved case definition interface.
"""

## AYOTTE/24SC

import os

import netCDF4 as nc
import numpy as np

from dephycf.Case import Case

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
        surfaceType="land",
        zorog=0.)

case.set_title("Forcing and initial conditions for AYOTTE-24SC case - Original definition")
case.set_reference("Ayotte et al. (1996, BLM) and files provided by F Hourdin initially from Ayotte")
case.set_author("F. couvreux")
case.set_script("DEPHY-SCM/AYOTTE/24SC/driver_DEF.py")

################################################
# 2. Initial state
################################################

# Surface pressure
ps = 100000.
case.add_init_ps(ps)

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
sfcForc= [     0,  270.096,    0,\
           25200,  270.096,    0]

timeSfc = sfcForc[0::3]

case.add_surface_fluxes(sens=sfcForc[1::3],lat=sfcForc[2::3],time=timeSfc,forc_wind='z0',z0=0.16)

################################################
# 4. Writing file
################################################

case.write('AYOTTE_24SC_DEF_driver.nc')

if lverbose:
    case.info()

################################################
# 5. Ploting, if asked
################################################

if lplot:
    case.plot(rep_images='./images/driver_DEF/',timeunits='hours')
