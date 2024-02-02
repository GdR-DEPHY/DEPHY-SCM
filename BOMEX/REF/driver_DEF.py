#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 27 November 2019

@author: Romain Roehrig

Modification
  2020/01/03, R. Roehrig: update for improved case definition interface.
"""

## BOMEX original case definition
## From ?? 

import os

import netCDF4 as nc
import numpy as np

from dephycf.Case import Case
from dephycf import constants as cc

################################################
# 0. General configuration of the present script
################################################

lplot = True # plot all the variables
lverbose = False # print information about variables and case

################################################
# 1. General information about the case
################################################

case = Case('BOMEX/REF',
        lat=15,
        lon=-56.5,
        startDate="19690624000000",
        endDate="19690625000000",
        surfaceType='ocean',
        zorog=0.)

case.set_title("Forcing and initial conditions for BOMEX case - Original definition")
case.set_reference("Siebesma and Cuijpers (JAS, 1995) and Grant and Brown (QJ, 1999)")
case.set_author("MP. Lefebvre")
case.set_script("DEPHY-SCM/BOMEX/REF/driver_DEF.py")

################################################
# 2. Initial state
################################################

# Surface pressure
ps = 101500.
case.add_init_ps(ps)

# Zonal and meridional wind
zu = [  0.,    700.,  3000.  ]
u  = [ -8.75,   -8.75,  -4.61]

zv = [ 0.,   700., 3000. ]
v  = [ 0.,    0.,  0.]

case.add_init_wind(u=u,v=v, ulev=zu, vlev=zv, levtype='altitude')

# Liquid-water Potential Temperature
zthetal = [  0.,     520.,  1480.,  2000., 3000.]
thetal  = [298.7,   298.7, 302.4,  308.2,   311.85]

case.add_init_thetal(thetal, lev=zthetal, levtype='altitude')

# Total water
zqt =[ 0.,  520., 1480., 2000., 3000. ] 
qt = [17.,   16.3,  10.7,   4.2,   3.0] # in g kg-1

case.add_init_qt(np.array(qt)/1000., lev=zqt, levtype='altitude') # converted in kg kg-1

# Turbulent Kinetic Energy
ztke = [0, 3000.]
nztke = len(ztke)
tke = np.zeros(nztke,dtype=np.float64)

for iz in range(0,nztke):
    if ztke[iz] <= 3000:
      tke[iz] = 1.-ztke[iz]/3000.
    else:
      tke[iz] = 0.

case.add_init_tke(tke, lev=ztke, levtype='altitude')

################################################
# 3. Forcing
################################################

# Constant geostrophic wind across the simulation
# Siebesma et Cuijpers donnent ug=-10.+0.0018*zz  vg=0.

zug = [0.,300.,500., 1500., 2100.,3000.]
nzug = len(zug)
ug = np.zeros(nzug,dtype=np.float64)
for iz in range(0,nzug):
    ug[iz] = -10.+1.8e-3**zug[iz]
    ug[iz] = -10.+1.8e-3*zug[iz]

vg = np.zeros(nzug,dtype=np.float64)

case.add_geostrophic_wind(ug=ug,vg=vg,lev=zug,levtype='altitude')

# Constant large-scale velocity - constant
zw = [0.,    1500, 2100.]
w  = [0., -0.0065,    0.]

case.add_vertical_velocity(w=w,lev=zw,levtype='altitude')

# Constant large-scale advection of potential temperature + radiative tendency 
zthetal_rad = [ 0.,  1500.,  3000.]
thetal_rad  = [-2.0,   -2.0,    0.] # in K day-1

case.add_thetal_radiation_tendency(np.array(thetal_rad)/86400.,
        lev=zthetal_rad, levtype='altitude') # converted in K s-1

# Constant large-scale advection of specific humidity
zqt_adv = [ 0.,    300.,     500.]
qt_adv  = [-1.2e-8, -1.2e-8,   0.] # in kg kg-1 s-1

case.add_qt_advection(qt_adv, lev=zqt_adv, levtype='altitude') 

# Surface Forcing
#            t(s)    H (W m-2)     LE (W m-2)  ustar (m s-1)
sfcForc = [    0., 8e-3*cc.Cpd,   5.2e-5*cc.Lv,  0.28,\
           86400., 8e-3*cc.Cpd,   5.3e-5*cc.Lv,  0.28]

timeSfc = sfcForc[0::3]

case.add_surface_fluxes(sens=8e-3*cc.Cpd, lat=5.2e-5*cc.Lv, 
        forc_wind='ustar', ustar=0.28)

case.add_surface_skin_temp(300.4)

################################################
# 4. Writing file
################################################

case.write('BOMEX_REF_DEF_driver.nc')

if lverbose:
    case.info()

################################################
# 5. Ploting, if asked
################################################

if lplot:
    case.plot(rep_images='./images/driver_DEF/',timeunits='hours')
