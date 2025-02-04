#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 10 December 2019

@author: Najda Villefranque

Comment: based on MESONH case but with active radiation
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

case = Case('ARMCU/MNHRAD',
        lat=36.6,
        lon=-97.5,
        startDate="19970621113000",
        endDate="19970622023000",
        surfaceType="land",
        zorog=314.)

case.set_title("Forcing and initial conditions for ARM-Cumulus case - Meso-NH definition")
case.set_reference("http://projects.knmi.nl/eurocs/ARM/case_ARM_html ; Brown et al. (2002, QJRMS)")
case.set_author("R. Roehrig, F. Couvreux, N. Villefranque")
case.set_script("DEPHY-SCM/ARMCU/MNHRAD/driver_DEF.py")

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

# +1 K/d to remove radiative cooling as it will be interactive
case.add_theta_advection(np.array(thadv)+1/(3600*24),time=timeF,lev=zadv,levtype='altitude',include_rad=False)

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

XTIMEF = np.zeros(31+1)
XSFTH = np.zeros(31+1)
XSFTQ = np.zeros(31+1)

XTIMEF[1] = 0.
XTIMEF[2] = 1800.
XTIMEF[3] = 3600.
XTIMEF[4] = 5400.
XTIMEF[5] = 7200.
XTIMEF[6] = 9000.
XTIMEF[7] = 10800.
XTIMEF[8] = 12600.
XTIMEF[9] = 14400.
XTIMEF[10] = 16200.
XTIMEF[11] = 18000.
XTIMEF[12] = 19800.
XTIMEF[13] = 21600.
XTIMEF[14] = 23400.
XTIMEF[15] = 25200.
XTIMEF[16] = 27000.
XTIMEF[17] = 28800.
XTIMEF[18] = 30600.
XTIMEF[19] = 32400.
XTIMEF[20] = 34200.
XTIMEF[21] = 36000.
XTIMEF[22] = 37800.
XTIMEF[23] = 39600.
XTIMEF[24] = 41400.
XTIMEF[25] = 43200.
XTIMEF[26] = 45000.
XTIMEF[27] = 46800.
XTIMEF[28] = 48600.
XTIMEF[29] = 50400.
XTIMEF[30] = 52200.
XTIMEF[31] = 54000.

XSFTH[1] = -30.
XSFTH[2] = -15.0
XSFTH[3] = 0.
XSFTH[4] = 15.0
XSFTH[5] = 30.0
XSFTH[6] = 45.0
XSFTH[7] = 60.0
XSFTH[8] = 75.0
XSFTH[9] = 90.0
XSFTH[10] = 100.0
XSFTH[11] = 110.0
XSFTH[12] = 120.0
XSFTH[13] = 130.0
XSFTH[14] = 140.0
XSFTH[15] = 140.0
XSFTH[16] = 140.0
XSFTH[17] = 132.0
XSFTH[18] = 124.0
XSFTH[19] = 116.0
XSFTH[20] = 108.0
XSFTH[21] = 100.0
XSFTH[22] = 78.0
XSFTH[23] = 56.0
XSFTH[24] = 34.0
XSFTH[25] = 12.0
XSFTH[26] = -10.0
XSFTH[27] = -10.0
XSFTH[28] = -10.0
XSFTH[29] = -10.0
XSFTH[30] = -10.0
XSFTH[31] = -10.0

XSFTQ[1] = 1.99936020473448506E-006  
XSFTQ[2] = 1.42454414587332047E-005  
XSFTQ[3] = 2.64915227127319252E-005  
XSFTQ[4] = 3.87376039667306440E-005  
XSFTQ[5] = 5.09836852207293695E-005  
XSFTQ[6] = 6.32297664747280883E-005  
XSFTQ[7] = 7.54758477287268071E-005  
XSFTQ[8] = 8.77219289827255259E-005  
XSFTQ[9] = 9.99680102367242446E-005  
XSFTQ[10] = 1.15962891874600132E-004  
XSFTQ[11] = 1.31957773512476006E-004  
XSFTQ[12] = 1.47952655150351879E-004  
XSFTQ[13] = 1.63947536788227780E-004  
XSFTQ[14] = 1.79942418426103654E-004  
XSFTQ[15] = 1.89939219449776085E-004  
XSFTQ[16] = 1.99936020473448489E-004  
XSFTQ[17] = 1.93538067818298140E-004  
XSFTQ[18] = 1.87140115163147790E-004  
XSFTQ[19] = 1.80742162507997441E-004  
XSFTQ[20] = 1.74344209852847091E-004  
XSFTQ[21] = 1.67946257197696742E-004  
XSFTQ[22] = 1.48752399232245693E-004  
XSFTQ[23] = 1.29558541266794618E-004  
XSFTQ[24] = 1.10364683301343569E-004  
XSFTQ[25] = 9.11708253358925209E-005  
XSFTQ[26] = 7.19769673704414589E-005
XSFTQ[27] = 5.39827255278310908E-005  
XSFTQ[28] = 3.59884836852207294E-005  
XSFTQ[29] = 1.79942418426103647E-005   
XSFTQ[30] = 0.0000000000000000      
XSFTQ[31] = -1.79942418426103647E-005

case.add_surface_fluxes(sens=XSFTH[1:],lat=XSFTQ[1:]*constants.Lv,time=XTIMEF[1:],forc_wind='z0',z0=0.035)

time_ts_rad = [0, 3600*8, 3600*15]
ts_rad = [295, 312.5, 300]
ts_rad = [295, 305, 295]
case.add_rad_ts(ts_rad, time=time_ts_rad)

################################################
# 4. Writing file
################################################

case.write('ARMCU_MNHRAD_DEF_driver.nc')

if lverbose:
    case.info()

################################################
# 5. Ploting, if asked
################################################

if lplot:
    case.plot(rep_images='./images/driver_DEF/',timeunits='hours')
