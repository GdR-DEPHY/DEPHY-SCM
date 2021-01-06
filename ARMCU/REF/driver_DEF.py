#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on 27 November 2019

@author: Romain Roehrig

Modification
  2020/11/11, R. Roehrig: update for improved case definition interface.
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

case = Case('ARMCU/REF',
        lat=36,
        lon=-97.5,
        startDate="19970621113000",
        endDate="19970622020000",
        surfaceType='land',
        zorog=314.)

case.set_title("Forcing and initial conditions for ARM-Cumulus case - Original definition")
case.set_reference("http://projects.knmi.nl/eurocs/ARM/case_ARM_html ; Brown et al. (2002, QJRMS)")
case.set_author("R. Roehrig")
case.set_script("DEPHY-SCM/ARMCU/REF/driver_DEF.py")

################################################
# 2. Initial state
################################################

# Surface pressure
ps = 97000.
case.add_init_ps(ps)

#         z (m) theta (K) rt (g kg-1) u (m s-1) v (m s-1)
init = [  0.0,   299.00,   15.20,      10.0,     0.0,\
         50.0,   301.50,   15.17,      10.0,     0.0,\
        350.0,   302.50,   14.98,      10.0,     0.0,\
        650.0,   303.53,   14.80,      10.0,     0.0,\
        700.0,   303.70,   14.70,      10.0,     0.0,\
       1300.0,   307.13,   13.50,      10.0,     0.0,\
       2500.0,   314.00,    3.00,      10.0,     0.0,\
       5500.0,   343.20,    3.00,      10.0,     0.0]

init = np.array(init,dtype=np.float64)

z = init[0::5]

case.add_init_theta(init[1::5], lev=z, levtype='altitude')
case.add_init_rt(init[2::5]/1000., lev=z, levtype='altitude')
case.add_init_wind(u=init[3::5],v=init[4::5], lev=z, levtype='altitude')

# Turbulent Kinetic Energy
ztke = range(0,6000+1,10)
nztke = len(ztke)
tke = np.zeros(nztke,dtype=np.float64)

for iz in range(0,nztke):
    if ztke[iz] < 150:
      tke[iz] = 0.15*(1.-ztke[iz]/150.)
    else:
      tke[iz] = 0.

case.add_init_tke(tke, lev=ztke, levtype='altitude')

################################################
# 3. Forcing
################################################

# Constant Geostrophic wind across the simulation
case.add_geostrophic_wind(ug=init[3::5],vg=init[4::5],lev=z,levtype='altitude')

# Advection forcing (+ radiative tendency)
#       t (s), A_theta (K hour-1) R_theta (K hour-1) A_rt (g kg-1 hour-1)
advF = [ 41400,      0.000,            -0.125,           0.080,\
         52200,      0.000,             0.000,           0.020,\
         63000,      0.000,             0.000,          -0.040,\
         73800,     -0.080,             0.000,          -0.100,\
         84600,     -0.160,             0.000,          -0.160,\
         93600,     -0.160,            -0.100,          -0.300]

advF = np.array(advF,dtype=np.float64)

timeF = advF[0::4] - 41400
ntf = len(timeF)
A_theta = advF[1::4]
R_theta = advF[2::4]
A_rt = advF[3::4]

zforc = range(0,6000+1,10)
nzf = len(zforc)
forc_theta = np.zeros((ntf,nzf),dtype=np.float64)
forc_rt = np.zeros((ntf,nzf),dtype=np.float64)

# 2000 (Brown et al. 2002, QJRMS) ou 3000 m (http://projects.knmi.nl/eurocs/ARM/case_ARM_html/ and used in Meso-NH)
# We take 3000 m here.

for it in range(0,ntf):
    for iz in range(0,nzf):
        if zforc[iz] < 1000.:
          forc_theta[it,iz] = A_theta[it]+R_theta[it]
          forc_rt[it,iz] = A_rt[it]
        elif zforc[iz] <= 3000. :
          forc_theta[it,iz] = (A_theta[it]+R_theta[it])*(1.-(zforc[iz]-1000.)/2000.)
          forc_rt[it,iz] = A_rt[it]*(1.-(zforc[iz]-1000.)/2000.)
        else:
          forc_theta[it,iz] = 0.
          forc_rt[it,iz] = 0.

case.add_theta_advection(forc_theta/3600.,time=timeF,lev=zforc,levtype='altitude',include_rad=True) # converted in K s-1
case.add_rt_advection(forc_rt/3600./1000.,time=timeF,lev=zforc,levtype='altitude') # converted in kg kg-1 s-1

# Surface Forcing
#            t (s) H (W m-2) LE (W m-2)
sfcForc= [  41400,  -30,       5,\
            55800,   90,     250,\
            64800,  140,     450,\
            68400,  140,     500,\
            77400,  100,     420,\
            86400,  -10,     180,\
            93600,  -10,       0]

sfcForc = np.array(sfcForc,dtype=np.float64)

timeSfc = sfcForc[0::3] - 41400 # Forcing time axis should be counted from starting date

case.add_surface_fluxes(sens=sfcForc[1::3],lat=sfcForc[2::3],time=timeSfc,forc_wind='z0',z0=0.035)

################################################
# 4. Writing file
################################################

case.write('ARMCU_REF_DEF_driver.nc',verbose=False)

if lverbose:
    case.info()

################################################
# 6. Ploting, if asked
################################################

if lplot:
    case.plot(rep_images='./images/driver_DEF/',timeunits='hours')
