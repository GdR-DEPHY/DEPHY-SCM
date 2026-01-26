#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 27 November 2019

@author: Romain Roehrig

Modification
  2020/11/11, R. Roehrig: update for improved case definition interface.
  2025/01/26, N. Villefranque: merge with MESONH and clean for publication.
"""

## ARM-Cumulus original case definition
## From http://projects.knmi.nl/eurocs/ARM/case_ARM_html/ [dead link]

import os

import netCDF4 as nc
import numpy as np

from datetime import datetime, timedelta
from dephycf.Case import Case
from dephycf import constants

################################################
# 0. General configuration of the present script
################################################

lplot    = True  # plot all the variables
lverbose = False # print information about variables and case

################################################
# 1. General information about the case
################################################

duration=14.5
tmin = datetime(1997, 6, 21, 11, 30)
tmax = tmin + timedelta(hours=duration)

case = Case('ARMCU/REF',
        lat=36,
        lon=-97.5,
        startDate=tmin,
        endDate=tmax,
        surfaceType='land',
        zorog=314.)

case.set_title("Forcing and initial conditions for ARM-Cumulus case - Original definition")
case.set_reference("Brown et al. (2002, QJRMS)")
case.set_author("R. Roehrig")
case.set_script("DEPHY-SCM/ARMCU/REF/driver_DEF.py")
case.set_comment("LES is initialized with white noise and TKE=0")
case.set_modifications("Forcings decrease linearly up to 3 km instead of 2 km in the reference paper.\nForcings were converted back and forth, causing minor differences in values.")

################################################
# 2. Initial state
################################################

# Surface pressure
ps = 97000.
case.add_init_ps(ps)

# Zonal and meridional wind
zu = zv = [   0, 5500.]
u  = [10, 10]
v  = [0,  0]

case.add_init_wind(u=u,ulev=zu,v=v,vlev=zv,levtype='altitude')

# Potential temperature
ztheta = [  0.0,   50.0,  350.0,  650.0,  700.0, 1300.0, 2500.0, 5500.0  ]
theta  = [299.00, 301.50, 302.50, 303.53, 303.70, 307.13, 314.00, 343.20 ] 

case.add_init_theta(theta,lev=ztheta,levtype='altitude')

# Total humidity
zrt = [  0.0,   50.0,  350.0,  650.0,  700.0, 1300.0, 2500.0, 5500.0  ]
rt =  [  15.20, 15.17, 14.98,  14.80,  14.70, 13.50,  3.00,   3.00    ] # g/kg

case.add_init_rt(np.array(rt)/1000.,lev=zrt,levtype='altitude') 

# Turbulent Kinetic Energy
ztke = [0, 150, 5500.]
nztke = len(ztke)
tke = np.zeros(nztke,dtype=np.float64)

for iz in range(0,nztke):
    if ztke[iz] < 150.:
      tke[iz] = 0.15*(1.-ztke[iz]/150.)
    else:
      tke[iz] = 0.

case.add_init_tke(tke, lev=ztke, levtype='altitude')

################################################
# 3. Forcing
################################################

zforc = [0., 1000., 3000., 5000.]
tforc = [41400., 52200., 63000., 73800., 84600., 86400+9000.]

tforc = np.array(tforc) - 41400
ntf = len(tforc)
nzf = len(zforc)

# Constant Geostrophic wind across the simulation
ug = np.zeros((ntf,nzf),dtype=float) + 10.
vg = np.zeros((ntf,nzf),dtype=float) + 0.

case.add_geostrophic_wind(ug=ug,vg=vg,lev=zforc,levtype='altitude',time=tforc)

###
# The paper gives the following forcings, constant over first 1000 m then
# linear decrease to zero up to 2000 m in Brown et al. 2002, and up to 3000 m :

# Advection forcing (+ radiative tendency)
#       t (s), A_theta (K hour-1) R_theta (K hour-1) A_rt (g kg-1 hour-1)
#advF = [ 41400,      0.000,            -0.125,           0.080,\
#         52200,      0.000,             0.000,           0.020,\
#         63000,      0.000,             0.000,          -0.040,\
#         73800,     -0.080,             0.000,          -0.100,\
#         84600,     -0.160,             0.000,          -0.160,\
#         93600,     -0.160,            -0.100,          -0.300]

###
# Below, forcings are converted from K/hour to K/day, each line corresponds to
# one instant :

# Potential temperature advection (ntf, nzf)
thadv = [[-3.47222E-05, -3.47222E-05, 0.,    0.],
         [ 0.,           0.,          0.,    0.],
         [ 0.,           0.,          0.,    0.],
         [-2.22222E-05, -2.22222E-05, 0.,    0.],
         [-4.44444E-05, -4.44444E-05, 0.,    0.],
         [-7.77777E-05, -7.77777E-05, 0.,    0.]]

case.add_theta_advection(thadv,time=tforc,lev=zforc,levtype='altitude',include_rad=True)

# Water vapor mixing ratio advection (ntf, nzf)
rtadv = [[ 2.22222E-08,  2.22222E-08, 0.,    0.],
         [ 5.55555E-09,  5.55555E-09, 0.,    0.],
         [-1.11111E-08, -1.11111E-08, 0.,    0.],
         [-2.77778E-08, -2.77778E-08, 0.,    0.],
         [-4.44444E-08, -4.44444E-08, 0.,    0.],
         [-9.11111E-08, -9.11111E-08, 0.,    0.]]

case.add_rt_advection(rtadv,time=tforc,lev=zforc,levtype='altitude')

###
# The paper gives the following surface forcings :

# Surface Forcing
#            t (s) H (W m-2) LE (W m-2)
#sfcForc= [  41400,  -30,       5,\
#            55800,   90,     250,\
#            64800,  140,     450,\
#            68400,  140,     500,\
#            77400,  100,     420,\
#            86400,  -10,     180,\
#            93600,  -10,       0]

###
# Below, surface forcings are read from a txt file, each line corresponds to
# one instant (every 30 min), LE is in kg/m2/s and will be converted to W/m2 in
# the add_surface_fluxes function:

# Surface Forcing
XTIMEF, XSFTQ, XSFTH = np.genfromtxt('../aux/surface_flux_forcings.txt', dtype=float, skip_header=0, usecols=[0,1,2]).transpose()
case.add_surface_fluxes(sens=XSFTH,lat=XSFTQ*constants.Lv,time=XTIMEF,forc_wind='z0',z0=0.035)

################################################
# 4. Writing file
################################################

case.write('ARMCU_REF_DEF_driver.nc')

if lverbose:
    case.info()

################################################
# 6. Ploting, if asked
################################################

if lplot:
    case.plot(rep_images='./images/driver_DEF/',timeunits='hours')
