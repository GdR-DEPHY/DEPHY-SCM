#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 01 April 2020

@author: Fleur Couvreux

Modifications:
  2020/04/07, R. Roehrig: some cleaning
  2021/01/03, R. Roehrig: update for improved case definition interface.
  2026/07/20, N. Villefranque: clean for publication.
"""

## AYOTTE/24SC original case definition
## From Ayotte et al., 1996, BML

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
# This is an idealized case so date and lat/lon are arbitrary

scase = "24SC"
fsens = 270.096

duration=7
tmin = datetime(2009, 12, 11, 10, 00)
tmax = tmin + timedelta(hours=duration)

case = Case('AYOTTE/%s'%scase,
        lat=45,
        lon=123.3,
        startDate=tmin,
        endDate=tmax,
        surfaceType="land",
        zorog=0.)

case.set_title("Forcing and initial conditions for AYOTTE/%s case - Original definition"%scase)
case.set_reference("Ayotte et al. (1996, BLM)")
case.set_author("F. Couvreux")
case.set_script("DEPHY-SCM/AYOTTE/%s/driver_DEF.py"%scase)
case.set_comment("Files provided by F. Hourdin initially from Ayotte")

################################################
# 2. Initial state
################################################

# Surface pressure
ps = 100000.
case.add_init_ps(ps)

# Initial profiles : read from aux file
z, theta, rt, u, v = np.genfromtxt("../aux/init_%s.txt"%scase, dtype=float).transpose()

# Potential temperature
case.add_init_theta(theta, lev=z, levtype='altitude')

# Total water mixing ratio
case.add_init_rt(rt/1000., lev=z, levtype='altitude') # converted in kg kg-1

# Zonal and meridional wind
case.add_init_wind(u=u, v=v, lev=z, levtype='altitude')

################################################
# 3. Forcing
################################################

# Constant Geostrophic wind across the simulation
ug = z*0. + 15.
vg = z*0. + 0.

case.add_geostrophic_wind(ug=ug,vg=vg,lev=z,levtype='altitude')

# No radiation
case.deactivate_radiation()

# Surface Forcing
timeSfc = [0, duration*3600]
sensib  = [fsens, fsens] # W/m2
latent  = [0,     0    ] # W/m2 
z0=0.16

case.add_surface_fluxes(sens=sensib, lat=latent, time=timeSfc,forc_wind='z0',z0=z0)

################################################
# 4. Writing file
################################################

case.write('AYOTTE_%s_DEF_driver.nc'%scase)

if lverbose:
    case.info()

################################################
# 5. Ploting, if asked
################################################

if lplot:
    case.plot(rep_images='./images/driver_DEF/',timeunits='hours')
