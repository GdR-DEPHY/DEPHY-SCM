#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on 27 November 2019

@author: Romain Roehrig
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
        zorog=0.,
        z0= 0.001019)

case.set_title("Forcing and initial conditions for BOMEX case - Original definition")
case.set_reference("Siebesma and Cuijpers (JAS, 1995) and Grant and Brown (QJ, 1999)")
case.set_author("MP. Lefebvre")
case.set_script("DEPHY-SCM/BOMEX/REF/driver_DEF.py")


# time units are expected to be seconds since startDate
t0 = 0 # 00:00 UTC, 24 June 1969
t1 = 86400 # 00:00 UTC, 25 June 1969


################################################
# 2. Initial state
################################################

# Surface pressure
ps = 101500.
case.add_variable('ps',[ps,])

# Zonal wind
zu = [   0,    700.,  20000.]
u  = [-8.75,  -8.75,   25.99]

case.add_variable('u',u,lev=zu,levtype='altitude')

# Meridional wind
zv = [ 0.,   20000. ]
v  = [ 0.,     0.]

case.add_variable('v',v,lev=zv,levtype='altitude')

# Potential Temperature
ztheta = [  0.,     520.,  1480.,  2000., 20000.]
theta  = [298.7,   298.7, 302.4,  308.2,   373.9]

case.add_variable('theta',theta,lev=ztheta,levtype='altitude')

# Specific humidity
# put 9000 m instead of 10000 m which seems more relevant 
# according to the formula at http://projects.knmi.nl/rico/setup1d_composite.html
zqv =[ 0.,  520., 1480., 2000., 20000.] 
qv = [17.,  16.3,  10.7,   4.2,   0.] # in g kg-1

case.add_variable('qv',np.array(qv)/1000.,lev=zqv,levtype='altitude') # converted in kg kg-1 (array type required)

# Turbulent Kinetic Energy
ztke = range(0,6000+1,10)
nztke = len(ztke)
tke = np.zeros(nztke,dtype=np.float64)

for iz in range(0,nztke):
    if ztke[iz] < 3000:
      tke[iz] = 1.-ztke[iz]/3000.
    else:
      tke[iz] = 0.

case.add_variable('tke',tke,lev=ztke,levtype='altitude')

################################################
# 3. Forcing
################################################

# Geostrophic wind across the simulation
# Siebesma et Cuijpers donnent ug=-10.+0.0018*zz  vg=0.

zug = range(0,6000+1,10)
nzug = len(zug)
ug = np.zeros(nzug,dtype=np.float64)
for iz in range(0,nzug):
    ug[iz] = -10.+0.0018*zug[iz]
    ug[iz] = -10.+0.0018*zug[iz]

vg = np.zeros(nzug,dtype=np.float64)
vg[iz] = 0.

case.add_variable('ug',[ug,ug],lev=zug,time=[t0,t1],levtype='altitude')
case.add_variable('vg',[vg,vg],lev=zug,time=[t0,t1],levtype='altitude')

# Surface Forcing
#            t(s)      H (W m-2)             LE (W m-2)
sfcForc = [   0., 9.4589217426112491,   153.03603975594794,\
          86400., 9.4589217426112491,   153.03603975594794]

sfcForc = np.array(sfcForc,dtype=np.float64)

timeSfc = sfcForc[0::3]

case.add_variable('sfc_sens_flx',sfcForc[1::3],time=timeSfc)
case.add_variable('sfc_lat_flx', sfcForc[2::3],time=timeSfc)

# Advection forcing (+ radiative tendency)

# Large-scale velocity - constant
zw = [0.,    1500, 2100.]
w  = [0., -0.0065,    0.]

case.add_variable('w',[w,w],time=[t0,t1],lev=zw,levtype='altitude')

# Large-scale advection of potential temperature + radiative tendency - constant 
zthetaadv = [ 0.,        1500.,  2500.]
thetaadv  = [-2.00016, -2.00016,   0.] # in K day-1

case.add_variable('thetal_adv',np.array([thetaadv,thetaadv])/86400.,time=[t0,t1],lev=zthetaadv,levtype='altitude') # converted in K s-1 (array type required)

# Large-scale advection of specific humidity - constant
zqvadv = [ 0.,      300.,  500.]
qvadv  = [-1.0368, -1.0368,  0.] # in g kg-1 day-1

case.add_variable('qt_adv',np.array([qvadv,qvadv])/86400./1000.,time=[t0,t1],lev=zqvadv,levtype='altitude') # converted in kg kg-1 s-1 (array type required)


################################################
# 4. Attributes
################################################

# advection of theta and rt
case.set_attribute("adv_thetal",1)
case.set_attribute("adv_qt",1)
# potential temperature radiative tendency is included in advection
case.set_attribute("rad_thetal","adv")
# Geostrophic wind forcing
case.set_attribute("forc_geo",1)
# Surface flux forcing, wind stress is computed using z0
case.set_attribute("surfaceType","land")
case.set_attribute("surfaceForcing","surfaceFlux")
case.set_attribute("surfaceForcingWind","z0")

################################################
# 5. Writing file
################################################

case.write('BOMEX_REF_DEF_driver.nc',verbose=False)

if lverbose:
    case.info()

################################################
# 6. Ploting, if asked
################################################

if lplot:
    case.plot(rep_images='./images/driver_DEF/',timeunits='hours')
