#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on 27 November 2019

@author: Romain Roehrig
"""

## ARM-Cumulus original case definition
## From http://projects.knmi.nl/eurocs/ARM/case_ARM_html/

import os
import sys
sys.path = ['../utils/',] + sys.path

import netCDF4 as nc
import numpy as np

import SCM_utils as utils

################################################
# 0. General configuration of the present script
################################################

lplot = True

################################################
# 1. General information about the case
################################################

case = utils.Case('ARMCU/REF',
        lat=36,
        lon=-97.5,
        startDate="19970621113000",
        endDate="19970622020000",
        zorog=314.,
        z0=0.035)

case.set_comment("Forcing and initial conditions for ARM-Cumulus case - Original definition")
case.set_reference("http://projects.knmi.nl/eurocs/ARM/case_ARM_html ; Brown et al. (2002, QJRMS)")
case.set_author("R. Roehrig")
case.set_script("DEPHY-SCM/ARMCU/setup_orig.py")


# time units are expected to be seconds since startDate
t0 = 0 # 11:30 UTC, 21 June 1997
t1 = 86400 + 2*3600 - 41400 # 02:00 UTC, 22 June 1997


################################################
# 2. Initial state
################################################

ps = 97000.
case.add_variable('ps',[ps,])

#         z (m) theta (K) rt (g kg-1) u (m s-1) v (m s-1)
init = [   0.0,   299.00,   15.20,      10.0,     0.0,\
         50.0,   301.50,   15.17,      10.0,     0.0,\
        350.0,   302.50,   14.98,      10.0,     0.0,\
        650.0,   303.53,   14.80,      10.0,     0.0,\
        700.0,   303.70,   14.70,      10.0,     0.0,\
       1300.0,   307.13,   13.50,      10.0,     0.0,\
       2500.0,   314.00,    3.00,      10.0,     0.0,\
       5500.0,   343.20,    3.00,      10.0,     0.0]

init = np.array(init,dtype=np.float64)

z = init[0::5]

case.add_variable('theta',init[1::5],      lev=z,levtype='altitude')
case.add_variable('rt',   init[2::5]/1000.,lev=z,levtype='altitude') # converted in kg kg-1
case.add_variable('u',    init[3::5],      lev=z,levtype='altitude')
case.add_variable('v',    init[4::5],      lev=z,levtype='altitude')

# Turbulent Kinetic Energy
ztke = range(0,6000+1,10)
nztke = len(ztke)
tke = np.zeros(nztke,dtype=np.float64)

for iz in range(0,nztke):
    if ztke[iz] < 150:
      tke[iz] = 0.15*(1.-ztke[iz]/150.)
    else:
      tke[iz] = 0.

case.add_variable('tke',tke,lev=ztke,levtype='altitude')

################################################
# 3. Forcing
################################################

# Constant Geostrophic wind across the simulation
ug = np.zeros((2,8),dtype=np.float64)
ug[0,:] = init[3::5]
ug[1,:] = init[3::5]

vg = np.zeros((2,8),dtype=np.float64)
vg[0] = init[4::5]
vg[1] = init[4::5]

case.add_variable('ug',ug,time=[t0,t1],lev=z,levtype='altitude')
case.add_variable('vg',vg,time=[t0,t1],lev=z,levtype='altitude')

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

timeSfc = sfcForc[0::3] - 41400

case.add_variable('sfc_sens_flx',sfcForc[1::3],time=timeSfc)
case.add_variable('sfc_lat_flx', sfcForc[2::3],time=timeSfc)

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

case.add_variable('thadv',forc_theta/3600.,   time=timeF,lev=zforc,levtype='altitude') # converted in K s-1
case.add_variable('rtadv',forc_rt/3600./1000.,time=timeF,lev=zforc,levtype='altitude') # converted in kg kg-1 s-1


################################################
# 4. Attributes
################################################

# advection of theta and rt
case.set_attribute("thadv",1)
case.set_attribute("rtadv",1)
# potential temperature radiative tendency is included in advection
case.set_attribute("thrad","adv")
# Geostrophic wind forcing
case.set_attribute("forc_geo",1)
# Surface flux forcing, wind stress is computed using z0
case.set_attribute("surfaceType","land")
case.set_attribute("surfaceForcing","surfaceFlux")
case.set_attribute("surfaceForcingWind","z0")

################################################
# 5. Writing file
################################################

case.write('ARMCU_REF_orig.nc',verbose=False)

case.info()

################################################
# 6. Ploting, if asked
################################################

if lplot:
    case.plot(rep_images='./images/setup_orig/',timeunits='hours')
