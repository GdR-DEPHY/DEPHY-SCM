#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on 06 Avril 2022

@author: Fleur Couvreux

from Kuang and Bretherton, 2006
"""

## BOMEX original case definition
## From ?? 

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

case = Case('KB2006',
        lat=15,
        lon=-56.5,
        startDate="20010927000000",
        endDate="20011002000000",
        surfaceType='ocean',
        zorog=0.)

case.set_title("Forcing and initial conditions for Kuang-Bretherton case initial state=BOMEX- Original definition")
case.set_reference("Kuang and Bretherton (JAS, 2006)")
case.set_author("F. Couvreux")
case.set_script("DEPHY-SCM/KUANG_BRETHERTON/REF/driver_DEF.py")

################################################
# 2. Initial state
################################################

# Surface pressure
ps = 101500.
case.add_init_ps(ps)

# Zonal and meridional wind
zu = [   0,    700.,  3000., 5600. ,20000.]
u  = [-8.75,  -8.75,   -4.61, 0.00, 0.00]

zv = [ 0.,   700., 3000., 5600., 20000. ]
v  = [ 0.,     0.,    0.,    0.,     0. ]

case.add_init_wind(u=u,v=v, ulev=zu, vlev=zv, levtype='altitude')

# Potential Temperature
zthetal = [  0.,     520.,  1480.,  2000., 3000., 4000., 20000.]
thetal  = [298.7,   298.7, 302.4,  308.2, 311.85, 315.5,  373.9]

case.add_init_thetal(thetal, lev=zthetal, levtype='altitude')

# Specific humidity
zqt =[ 0.,  520., 1480., 2000., 3000., 4000., 20000.] 
qt = [17.,  16.57,  10.82, 4.22, 3.01,    0., 0.] # in g kg-1

case.add_init_qt(np.array(qt)/1000., lev=zqt, levtype='altitude') # converted in kg kg-1

# Turbulent Kinetic Energy
ztke = range(0,6000+1,10)
nztke = len(ztke)
tke = np.zeros(nztke,dtype=np.float64)

#for iz in range(0,nztke):
#    if ztke[iz] < 3000:
#      tke[iz] = 1.-ztke[iz]/3000.
#    else:
#      tke[iz] = 0.

case.add_init_tke(tke, lev=ztke, levtype='altitude')

################################################
# 3. Forcing
################################################

# Constant geostrophic wind across the simulation
# Siebesma et Cuijpers donnent ug=-10.+0.0018*zz  vg=0.

zug = [0., 300., 500., 1500., 2100., 2500.,4000., 5600., 20000.]
ug =  [-10.,-9.46, -9.10, -7.30, -6.22, -5.50, -2.80, 0.0, 0.0]
nzug=len(zug)
vg = np.zeros(nzug,dtype=np.float64)

case.add_geostrophic_wind(ug=ug,vg=vg,lev=zug,levtype='altitude')

# Constant large-scale velocity - constant
zw = zug
w  = [0., -0.0013, -0.00217, -0.0065, 0., 0., 0., 0.,    0.]

case.add_vertical_velocity(w=w,lev=zw,levtype='altitude')

# Constant large-scale advection of potential temperature + radiative tendency 
zthetaladv = zug 
thetaladv  = [-2.315e-5, -2.315e-5, -2.315e-5, -2.315e-5, -0.926e-5, 0., 0., 0., 0.] # in K s-1

case.add_thetal_advection(np.array(thetaladv),lev=zthetaladv,levtype='altitude',include_rad=True) # converted in K s-1

# Constant large-scale advection of specific humidity
zqtadv = zug 
qtadv  = [-1.2e-5, -1.2e-5, 0., 0., 0., 0., 0., 0.,  0.] # in g kg-1 s-1

case.add_qt_advection(np.array(qtadv)/1000.,lev=zqtadv,levtype='altitude') # converted in kg kg-1 s-1

# Surface Forcing
#            t(s)      H (W m-2)         LE (W m-2)
sfcForc = [   0.,         10.,            150.0,\
          21600.,         10.,            150.0,\
          43200.,         10.,            150.0,\
          64800.,         12.5,           187.5,\
          86400.,         15.,            225.0,\
         108000.,         17.5,           262.5,\
         129600.,         20.,            300.0,\
         151200.,         20.,            300.0,\
         172800.,         20.,            300.0,\
         194400.,         20.,            300.0,\
         216000.,         20.,            300.0,\
         237600.,         20.,            300.0,\
         259200.,         20.,            300.0,\
         280800.,         20.,            300.0,\
         302400.,         20.,            300.0,\
         324000.,         20.,            300.0,\
         345600.,         20.,            300.0,\
         367200.,         20.,            300.0,\
         388800.,         20.,            300.0,\
         410400.,         20.,            300.0,\
         432000.,         20.,            300.0]

timeSfc = sfcForc[0::3]

case.add_surface_fluxes(sens=sfcForc[1::3],lat=sfcForc[2::3],time=timeSfc,forc_wind='z0',z0=0.001019)

################################################
# 4. Writing file
################################################

case.write('KB2006_REF_DEF_driver.nc')

if lverbose:
    case.info()

################################################
# 5. Ploting, if asked
################################################

if lplot:
    case.plot(rep_images='./images/driver_DEF/',timeunits='hours')
