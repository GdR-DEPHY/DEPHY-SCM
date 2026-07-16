#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 06 Avril 2022

@author: Fleur Couvreux

Modification
  2026/07/15, N. Villefranque: clean for publication.
"""

## KB2006 original case definition
# from Kuang and Bretherton, 2006

import numpy as np

from datetime import datetime, timedelta
from dephycf.Case import Case

################################################
# 0. General configuration of the present script
################################################

lplot = True # plot all the variables
lverbose = False # print information about variables and case

################################################
# 1. General information about the case
################################################

duration=5*24
tmin = datetime(2001, 9, 27)
tmax = tmin + timedelta(hours=duration)

case = Case('KB2006',
        lat=15,
        lon=-56.5,
        startDate=tmin,
        endDate=tmax,
        surfaceType='ocean',
        zorog=0.)

case.set_title("Forcing and initial conditions for Kuang-Bretherton case - Original definition")
case.set_reference("Kuang and Bretherton (JAS, 2006)")
case.set_author("F. Couvreux")
case.set_script("DEPHY-SCM/KB2006/REF/driver_DEF.py")
case.set_comment("Initial state as in BOMEX")
case.set_modifications("Wind profiles as in BOMEX instead of null as in Kuang and Bretherton, and with geostrophic forcings.\nModified initial theta profile with an inversion at tropopause.")

################################################
# 2. Initial state
################################################

# Surface pressure
ps = 101500.
case.add_init_ps(ps)

# Zonal and meridional wind
zu = zv = [  0.,    700.,  3000., 5600, 10000]
u  = [ -8.75, -8.75, -4.61, 0., 0.]
v  = [  0.,    0.,    0.,   0., 0.]

case.add_init_wind(u=u,v=v, ulev=zu, vlev=zv, levtype='altitude')

# Liquid-water Potential Temperature
zthetal = [  0.,     520.,  1480.,  2000., 3000., 4000., 15000., 17500., 20000.]
thetal  = [298.7,   298.7, 302.4,  308.2, 311.85, 316.93, 362.95, 392.25, 470.56]

case.add_init_thetal(thetal, lev=zthetal, levtype='altitude')

# Total water
zqt =[ 0., 520., 1480., 2000., 3000., 4000., 20000.] 
qt = [17., 16.3, 10.7,   4.2,   3.0,  0.,     0] # in g kg-1

case.add_init_qt(np.array(qt)/1000., lev=zqt, levtype='altitude') # converted in kg kg-1

# Turbulent Kinetic Energy
ztke = [0, 3000., 6000.]
nztke = len(ztke)
tke = np.zeros(nztke,dtype=float)

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
# add points at z = 5600 and z = 10000 m with u = 0 m/s

zug = [0., 300., 500., 1500., 2100., 3000., 5600, 10000]
nzug = len(zug)
ug = np.zeros(nzug,dtype=float)
for iz in range(0,nzug):
  ug[iz] = -10.+1.8e-3*zug[iz]
  if zug[iz] > 5000: ug[iz] = 0.

vg = np.zeros(nzug,dtype=float)

case.add_geostrophic_wind(ug=ug,vg=vg,lev=zug,levtype='altitude')


# Constant large-scale velocity - constant
zw = [0.,  300.,    500.,     1500., 2100.]
w  = [0., -0.0013, -0.00217, -0.0065, 0.  ]

case.add_vertical_velocity(w=w,lev=zw,levtype='altitude')

# Constant large-scale advection of potential temperature + radiative tendency 
zthetaladv = [0.,         1500.,     2100.,    2500.]
thetaladv  = [-2.315e-5, -2.315e-5, -0.926e-5, 0.   ] # in K s-1

case.add_thetal_advection(np.array(thetaladv),lev=zthetaladv,levtype='altitude',include_rad=True) # converted in K s-1

# Constant large-scale advection of specific humidity
zqtadv = [0,        300,    500]
qtadv  = [-1.2e-8, -1.2e-8, 0. ] # in kg kg-1 s-1

case.add_qt_advection(qtadv, lev=zqtadv, levtype='altitude')

# Surface Forcing
ustar  = 0.28
time, sensib, latent = np.genfromtxt('../aux/surface_flux_forcings.txt',
                                     dtype=float).transpose()

case.add_surface_fluxes(sens=sensib, lat=latent, time=time,
                        forc_wind='ustar',ustar=ustar)

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
