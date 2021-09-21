#!/usr/bin/env python
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
zu = [   0,    700.,  20000.]
u  = [-8.75,  -8.75,   25.99]

zv = [ 0.,   20000. ]
v  = [ 0.,     0.]

case.add_init_wind(u=u,v=v, ulev=zu, vlev=zv, levtype='altitude')

# Potential Temperature
ztheta = [  0.,     520.,  1480.,  2000., 20000.]
theta  = [298.7,   298.7, 302.4,  308.2,   373.9]

case.add_init_theta(theta, lev=ztheta, levtype='altitude')

# Specific humidity
zqv =[ 0.,  520., 1480., 2000., 20000.] 
qv = [17.,  16.3,  10.7,   4.2,   0.] # in g kg-1

case.add_init_qv(np.array(qv)/1000., lev=zqv, levtype='altitude') # converted in kg kg-1

# Turbulent Kinetic Energy
ztke = range(0,6000+1,10)
nztke = len(ztke)
tke = np.zeros(nztke,dtype=np.float64)

for iz in range(0,nztke):
    if ztke[iz] < 3000:
      tke[iz] = 1.-ztke[iz]/3000.
    else:
      tke[iz] = 0.

case.add_init_tke(tke, lev=ztke, levtype='altitude')

################################################
# 3. Forcing
################################################

# Constant geostrophic wind across the simulation
# Siebesma et Cuijpers donnent ug=-10.+0.0018*zz  vg=0.

zug = range(0,6000+1,10)
nzug = len(zug)
ug = np.zeros(nzug,dtype=np.float64)
for iz in range(0,nzug):
    ug[iz] = -10.+0.0018*zug[iz]
    ug[iz] = -10.+0.0018*zug[iz]

vg = np.zeros(nzug,dtype=np.float64)

case.add_geostrophic_wind(ug=ug,vg=vg,lev=zug,levtype='altitude')

# Constant large-scale velocity - constant
zw = [0.,    1500, 2100.]
w  = [0., -0.0065,    0.]

case.add_vertical_velocity(w=w,lev=zw,levtype='altitude')

# Constant large-scale advection of potential temperature + radiative tendency 
zthetaladv = [ 0.,        1500.,  2500.]
thetaladv  = [-2.00016, -2.00016,   0.] # in K day-1

case.add_thetal_advection(np.array(thetaladv)/86400.,lev=zthetaladv,levtype='altitude',include_rad=True) # converted in K s-1

# Constant large-scale advection of specific humidity
zqtadv = [ 0.,      300.,  500.]
qtadv  = [-1.0368, -1.0368,  0.] # in g kg-1 day-1

case.add_qt_advection(np.array(qtadv)/86400./1000.,lev=zqtadv,levtype='altitude') # converted in kg kg-1 s-1

# Surface Forcing
#            t(s)      H (W m-2)             LE (W m-2)
sfcForc = [   0., 9.4589217426112491,   153.03603975594794,\
          86400., 9.4589217426112491,   153.03603975594794]

timeSfc = sfcForc[0::3]

case.add_surface_fluxes(sens=sfcForc[1::3],lat=sfcForc[2::3],time=timeSfc,forc_wind='z0',z0=0.001019)

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
