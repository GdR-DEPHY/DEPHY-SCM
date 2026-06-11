#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 10 December 2019

@author: Romain Roehrig

Comment: based on Meso-NH namelists ARMCU_EXSEG1.nam and ARMCU_PRE_IDEA1.nam

Modification
  2025/11/13, N. Villefranque: clean for publication
  2021/01/03, R. Roehrig: update for improved case definition interface.
"""

## ARM-Cumulus case definition corresponding to MESONH LES setup, with TKE

import os

import netCDF4 as nc
import numpy as np

from dephycf import constants
from dephycf.Case import Case

################################################
# 0. General configuration of the present script
################################################

lplot = True # plot all the variables
lverbose = False # print information about variables and case

################################################
# 1. General information about the case
################################################

case = Case('ARMCU/MESONH-withTKE',
        lat=36.6,
        lon=-97.5,
        startDate="19970621113000",
        endDate="19970622023000",
        surfaceType="land",
        zorog=314.)

case.set_title("Forcing and initial conditions for ARM-Cumulus case - Meso-NH definition")
case.set_reference("Brown et al. (2002, QJRMS)")
case.set_author("R. Roehrig, F. Couvreux")
case.set_script("DEPHY-SCM/ARMCU/MESONH-withTKE/driver_DEF.py")
#case.set_comment("http://projects.knmi.nl/eurocs/ARM/case_ARM_html") # CE LIEN EST MORT
case.set_modifications("TKE is initialized to 0; forcings are applied up to 3 km")

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

# Specific humidity
zrv = [  0.0,   50.0,  350.0,  650.0,  700.0, 1300.0, 2500.0, 5500.0  ]
rv =  [  15.20e-3, 15.17e-3, 14.98e-3, 14.80e-3, 14.70e-3, 13.50e-3, 3.00e-3, 3.00e-3] # converted in kg kg-1 

case.add_init_rv(rv,lev=zrv,levtype='altitude') 

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

# Potential temperature advection (ntf, nzf)
thadv = [[-3.47222E-05, -3.47222E-05, 0.,    0.],
         [ 0.,           0.,          0.,    0.],
         [ 0.,           0.,          0.,    0.],
         [-2.22222E-05, -2.22222E-05, 0.,    0.],
         [-4.44444E-05, -4.44444E-05, 0.,    0.],
         [-7.77777E-05, -7.77777E-05, 0.,    0.]]

case.add_theta_advection(thadv,time=tforc,lev=zforc,levtype='altitude',include_rad=True)

# Water vapor mixing ratio advection (ntf, nzf)
rvadv = [[ 2.22222E-08,  2.22222E-08, 0.,    0.],
         [ 5.55555E-09,  5.55555E-09, 0.,    0.],
         [-1.11111E-08, -1.11111E-08, 0.,    0.],
         [-2.77778E-08, -2.77778E-08, 0.,    0.],
         [-4.44444E-08, -4.44444E-08, 0.,    0.],
         [-9.11111E-08, -9.11111E-08, 0.,    0.]]

case.add_rv_advection(rvadv,time=tforc,lev=zforc,levtype='altitude')

# Surface Forcing
XTIMEF, XSFTQ, XSFTH = np.genfromtxt('surface_flux_forcings.txt', dtype=float, skip_header=0, usecols=[0,1,2]).transpose()

case.add_surface_fluxes(sens=XSFTH,lat=XSFTQ*constants.Lv,time=XTIMEF,forc_wind='z0',z0=0.035)

################################################
# 4. Writing file
################################################

case.write('ARMCU_MESONH-withTKE_DEF_driver.nc')

if lverbose:
    case.info()

################################################
# 5. Ploting, if asked
################################################

if lplot:
    case.plot(rep_images='./images/driver_DEF/',timeunits='hours')
