#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 27 November 20199 June 2026

@author: Benoît Vié
"""

## DEPHY cirrus case definition

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

case = Case('CIRRUS/GPT',
        lat=50,
        lon=0,
        startDate="20260609000000",
        endDate="20260609120000",
        surfaceType='ocean',
            zorog=0.01)

case.set_title("Forcing and initial conditions for cirrus case - Original definition")
case.set_reference("")
case.set_author("B. Vié")
case.set_script("DEPHY-SCM/CIRRUS/GPT/driver_DEF.py")

################################################
# 2. Initial state
################################################

# Surface pressure
ps = 100000.
case.add_init_ps(ps)

# Altitude (m)
z = np.array([
      0,  500, 1000, 1500, 2000, 2500, 3000,
   4000, 5000, 6000, 7000, 8000, 9000,
    10000,11000,12000, 13000
])

# Température potentielle θ (K)
theta = np.array([
    290.0, 292.0, 294.0, 296.0, 298.0, 300.0, 302.0,
    306.0, 311.0, 317.0, 323.0, 330.0, 338.0,
    347.0, 356.0, 365.0, 370.0
])

# Vent horizontal (module) (m/s)
u = np.array([
     4,  6,  8, 10, 12, 14, 16,
    20, 24, 28, 34, 40, 48,
    55, 60, 55, 50
])

# Rapport de mélange vapeur d'eau (g/kg)
rv = np.array([
    7.0, 6.2, 5.5, 4.8, 4.0, 3.2, 2.55,
    1.75, 1.2, 0.9, 0.65, 0.5, 0.45,
    0.85, 0.25, 0.05, 0.
])

case.add_init_theta(theta, lev=z, levtype='altitude')
case.add_init_rt(rv/1000., lev=z, levtype='altitude')
case.add_init_wind(u,u*0., lev=z, levtype='altitude')

# Turbulent Kinetic Energy
ztke = z
nztke = len(ztke)
tke = np.zeros(nztke,dtype=np.float64)

case.add_init_tke(tke, lev=ztke, levtype='altitude')

################################################
# 3. Forcing
################################################


# Surface Forcing
#            t (s) H (W m-2) LE (W m-2)
case.add_surface_fluxes(sens=0.,lat=0.,forc_wind='z0',z0=0.01)


# add zero tendencies 
case.add_theta_advection([0]*len(z), lev=z, levtype="altitude", include_rad=False)


################################################
# 4. Writing file
################################################

case.write('CIRRUS_GPT_DEF_driver.nc')

if lverbose:
    case.info()

################################################
# 6. Ploting, if asked
################################################

if lplot:
    case.plot(rep_images='./images/driver_DEF/',timeunits='hours')
