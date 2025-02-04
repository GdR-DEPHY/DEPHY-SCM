#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 04 February 2025

@author: Najda Villefranque

Comment: based on ARMCU/MESONH but with surface fluxes forced from ARMCU/MNHRADSURF

Modification
"""

## ARM-Cumulus original case definition
## From http://projects.knmi.nl/eurocs/ARM/case_ARM_html/

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

case = Case('ARMCU/MNHRADFRCSURF',
        lat=36.6,
        lon=-97.5,
        startDate="19970621113000",
        endDate="19970622023000",
        surfaceType="land",
        zorog=314.)

case.set_title("Forcing and initial conditions for ARM-Cumulus case - modified surface fluxes")
case.set_reference("http://projects.knmi.nl/eurocs/ARM/case_ARM_html ; Brown et al. (2002, QJRMS)")
case.set_author("R. Roehrig, F. Couvreux, N. Villefranque")
case.set_script("DEPHY-SCM/ARMCU/MNHRADFRCSURF/driver_DEF.py")

################################################
# 2. Initial state
################################################

# Surface pressure
ps = 97000.
case.add_init_ps(ps)

#         z (m) theta (K) rv (kg kg-1)
init = [  0.0,   299.00,   15.20e-3,\
         50.0,   301.50,   15.17e-3,\
        350.0,   302.50,   14.98e-3,\
        650.0,   303.53,   14.80e-3,\
        700.0,   303.70,   14.70e-3,\
       1300.0,   307.13,   13.50e-3,\
       2500.0,   314.00,    3.00e-3,\
       5500.0,   343.20,    3.00e-3]

init = np.array(init,dtype=np.float64)

z = init[0::3]

case.add_init_theta(init[1::3], lev=z, levtype='altitude')
case.add_init_rv(init[2::3], lev=z, levtype='altitude')

#         z (m) u (m s-1) v (m s-1)
init = [  0.0,      10.0,     0.0,\
       5500.0,      10.0,     0.0]

init = np.array(init,dtype=np.float64)

z = init[0::3]

case.add_init_wind(u=init[1::3],v=init[2::3], lev=z, levtype='altitude')

################################################
# 3. Forcing
################################################

# Constant Geostrophic wind across the simulation
zforc = [0., 1000., 3000., 5000.]
timeF = [41400., 52200., 63000., 73800., 84600., 86400+9000.]

timeF = np.array(timeF) - 41400
ntf, = timeF.shape
nzf = len(zforc)

ug = np.zeros((ntf,nzf),dtype=np.float64) + 10.
vg = np.zeros((ntf,nzf),dtype=np.float64) + 0.

case.add_geostrophic_wind(ug=ug,vg=vg,lev=zforc,levtype='altitude',time=timeF)

# Potential temperature advection
zadv = [                  0.,        1000.,       3000., 5000.]
tmp =  [41400.,       -3.47222E-05, -3.47222E-05, 0.,    0.,\
        52200.,        0.,           0.,          0.,    0.,\
        63000.,        0.,           0.,          0.,    0.,\
        73800.,       -2.22222E-05, -2.22222E-05, 0.,    0.,\
        84600.,       -4.44444E-05, -4.44444E-05, 0.,    0.,\
        86400.+9000., -7.77777E-05, -7.77777E-05, 0.,    0.]


timeF = np.array(tmp[0::5]) - 41400
ntf, = timeF.shape
nzf = len(zadv)

thadv = np.zeros((ntf,nzf),dtype=np.float64)
for it in range(0,ntf):
    thadv[it,:] = tmp[5*it+1:5*it+5]

case.add_theta_advection(np.array(thadv)+1/(3600*24),time=timeF,lev=zadv,levtype='altitude',include_rad=False)
#case.add_theta_advection(thadv,time=timeF,lev=zadv,levtype='altitude',include_rad=True)

# Water vapor mixing ratio advection
zadv = [                  0.,        1000.,       3000., 5000.]
tmp =  [41400.,        2.22222E-08,  2.22222E-08, 0.,    0.,\
        52200.,        5.55555E-09,  5.55555E-09, 0.,    0.,\
        63000.,       -1.11111E-08, -1.11111E-08, 0.,    0.,\
        73800.,       -2.77778E-08, -2.77778E-08, 0.,    0.,\
        84600.,       -4.44444E-08, -4.44444E-08, 0.,    0.,\
        86400.+9000., -9.11111E-08, -9.11111E-08, 0.,    0.]


timeF = np.array(tmp[0::5]) - 41400
ntf, = timeF.shape
nzf = len(zadv)

rvadv = np.zeros((ntf,nzf),dtype=np.float64)
for it in range(0,ntf):
    rvadv[it,:] = tmp[5*it+1:5*it+5]

case.add_rv_advection(rvadv,time=timeF,lev=zadv,levtype='altitude')

# Surface Forcing

XTIMEF, XSFTQ, XSFTH = np.genfromtxt('surface_flux_forcings.txt', dtype=float, skip_header=0, usecols=[0,1,2]).transpose()
XTIMEF = XTIMEF-XTIMEF[0]
XSFTH = XSFTH*constants.Cpd

case.add_surface_fluxes(sens=XSFTH,lat=XSFTQ*constants.Lv,time=XTIMEF,forc_wind='z0',z0=0.035)

################################################
# 4. Writing file
################################################

case.write('ARMCU_MNHRADFRCSURF_DEF_driver.nc')

if lverbose:
    case.info()

################################################
# 5. Ploting, if asked
################################################

if lplot:
    case.plot(rep_images='./images/driver_DEF/',timeunits='hours')
