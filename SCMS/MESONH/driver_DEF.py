#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on 5 November 2020

@author: Fleur Couvreux
"""

## SCMS original case definition
## From Neggers et al 2003

import os
import sys
sys.path = ['../../utils/',] + sys.path

import netCDF4 as nc
import numpy as np

import constants
from Case import Case

################################################
# 0. General configuration of the present script
################################################

lplot = True # plot all the variables
lverbose = False # print information about variables and case

################################################
# 1. General information about the case
################################################

case = Case('SCMS/MESONH',
        lat=28.7,
        lon=-81.,
        startDate="19950805120000",
        endDate="19950806000000",
        zorog=0.,
        z0=0.035)

case.set_title("Forcing and initial conditions for SCMS-Cumulus case - Meso-NH definition")
case.set_reference(" Neggers et al. (2003, QJRMS)")
case.set_author("F. Couvreux, N. Villefranque")
case.set_script("DEPHY-SCM/SCMS/MESONH/driver_MESONH.py")


# time units are expected to be seconds since startDate
t0 = 0 # 12:00 UTC, 05 August 1995
t1 = 43200 # 00:00 UTC, 06 August 1995


################################################
# 2. Initial state
################################################

# Surface pressure
ps = 102100.
case.add_variable('ps',[ps,])
nz=126
z = np.zeros((nz),dtype=np.float64)
th = np.zeros((nz),dtype=np.float64)
rv = np.zeros((nz),dtype=np.float64)

z[:]=[0.0,20.0,60.0,100.0,140.0,180.0,220.0,260.0,300.0,340.0,380.0,420.0,460.0,500.0,540.0,580.0,620.0,660.0,700.0,740.0,780.0,820.0,860.0,900.0,940.0,980.0,1020.0,1060.0,1100.0,1140.0,1180.0,1220.0,1260.0,1300.0,1340.0,1380.0,1420.0,1460.0,1500.0,1540.0,1580.0,1620.0,1660.0,1700.0,1740.0,1780.0,1820.0,1860.0,1900.0,1940.0,1980.0,2020.0,2060.0,2100.0,2140.0,2180.0,2220.0,2260.0,2300.0,2340.0,2380.0,2420.0,2460.0,2500.0,2540.0,2580.0,2620.0,2660.0,2700.0,2740.0,2780.0,2820.0,2860.0,2900.0,2940.0,2980.0,3020.0,3060.0,3100.0,3140.0,3180.0,3220.0,3260.0,3300.0,3340.0,3380.0,3420.0,3460.0,3500.0,3540.0,3580.0,3620.0,3660.0,3700.0,3740.0,3780.0,3820.0,3860.0,3900.0,3940.0,3980.0,4020.0,4060.0,4100.0,4140.0,4180.0,4220.0,4260.0,4300.0,4340.0,4380.0,4420.0,4460.0,4500.0,4540.0,4580.0,4620.0,4660.0,4700.0,4740.0,4780.0,4820.0,4860.0,4900.0,4940.0,4980.0]
th[:]=[297.2,297.2,297.4,297.6,297.7,297.9,298.1,298.3,298.4,298.6,298.8,298.9,299.1,299.3,299.4,299.6,299.8,299.9,300.1,300.3,300.4,300.6,300.8,300.9,301.1,301.3,301.4,301.6,301.8,301.9,302.1,302.3,302.5,302.6,302.8,303.0,303.1,303.3,303.5,303.9,304.3,304.7,305.1,305.5,305.9,306.3,306.7,307.1,307.5,307.9,308.3,308.7,309.1,309.5,309.9,310.3,310.6,310.7,310.8,310.8,310.9,311.0,311.1,311.2,311.2,311.3,311.4,311.5,311.6,311.6,311.7,311.8,311.9,312.0,312.0,312.1,312.2,312.3,312.4,312.4,312.5,312.6,312.7,312.8,312.8,312.9,313.0,313.1,313.2,313.2,313.3,313.4,313.5,313.6,313.6,313.7,313.8,313.9,314.0,314.0,314.1,314.2,314.3,314.4,314.4,314.5,314.6,314.7,314.8,314.8,314.9,315.0,315.1,315.2,315.2,315.3,315.4,315.5,315.6,315.6,315.7,315.8,315.9,316.0,316.0,316.1]
rv[:]=[0.0178,0.0178,0.0177,0.0177,0.0176,0.0176,0.0175,0.0175,0.0174,0.0174,0.0173,0.0173,0.0172,0.0172,0.0171,0.0171,0.0170,0.0170,0.0169,0.0169,0.0168,0.0167,0.0165,0.0163,0.0161,0.0158,0.0156,0.0154,0.0152,0.0150,0.0148,0.0146,0.0144,0.0142,0.0140,0.0138,0.0136,0.0134,0.0132,0.0125,0.0118,0.0111,0.0104,0.0098,0.0091,0.0084,0.0077,0.0070,0.0064,0.0057,0.0050,0.0044,0.0037,0.0030,0.0030,0.0030,0.0030,0.0030,0.0030,0.0030,0.0030,0.0030,0.0030,0.0030,0.0030,0.0030,0.0030,0.0030,0.0030,0.0030,0.0030,0.0030,0.0030,0.0030,0.0030,0.0030,0.0030,0.0030,0.0030,0.0030,0.0030,0.0030,0.0030,0.0030,0.0030,0.0030,0.0030,0.0030,0.0030,0.0030,0.0030,0.0030,0.0030,0.0030,0.0030,0.0030,0.0030,0.0030,0.0030,0.0030,0.0030,0.0030,0.0030,0.0030,0.0030,0.0030,0.0030,0.0030,0.0030,0.0030,0.0030,0.0030,0.0030,0.0030,0.0030,0.0030,0.0030,0.0030,0.0030,0.0030,0.0030,0.0030,0.0030,0.0030,0.0030,0.0030]

case.add_variable('theta',th,lev=z,levtype='altitude')
case.add_variable('rv',   rv,lev=z,levtype='altitude')

#         z (m) u (m s-1) v (m s-1)
init = [  0.0,      -4.0,     4.0,\
       4980.0,      -4.0,     4.0]

init = np.array(init,dtype=np.float64)

z = init[0::3]

case.add_variable('u',    init[1::3],      lev=z,levtype='altitude')
case.add_variable('v',    init[2::3],      lev=z,levtype='altitude')

################################################
# 3. Forcing
################################################

# Constant Geostrophic wind across the simulation
zforc = [0., 1500., 5000.]
timeF = [0., 43200.]

ntf = len(timeF)
nzf = len(zforc)

ug = np.zeros((ntf,nzf),dtype=np.float64)
vg = np.zeros((ntf,nzf),dtype=np.float64)
for it in range(0,ntf):
    ug[it,:] = [-4., -4., -4.]
    vg[it,:] = [ 4.,  4.,  4.]

print(ug)
print(vg)
case.add_variable('ug',ug,time=timeF,lev=zforc,levtype='altitude')
case.add_variable('vg',vg,time=timeF,lev=zforc,levtype='altitude')

# Surface Forcing

XTIMEF = np.zeros(13+1)
XSFTH = np.zeros(13+1)
XSFTQ = np.zeros(13+1)

XTIMEF[1]=0. 
XTIMEF[2]=3600. 
XTIMEF[3]=7200. 
XTIMEF[4]=10800. 
XTIMEF[5]=14400. 
XTIMEF[6]=18000. 
XTIMEF[7]=21600. 
XTIMEF[8]=25200. 
XTIMEF[9]=28800. 
XTIMEF[10]=32400. 
XTIMEF[11]=36000. 
XTIMEF[12]=39600. 
XTIMEF[13]=43200. 
XSFTH[1]=0.0000000 
XSFTH[2]=25.881905 
XSFTH[3]=50.000000 
XSFTH[4]=70.710678 
XSFTH[5]=86.602539 
XSFTH[6]=96.592590 
XSFTH[7]=100.000000 
XSFTH[8]=96.592583  
XSFTH[9]=86.602531 
XSFTH[10]=70.710678 
XSFTH[11]=50.000008 
XSFTH[12]=25.881893 
XSFTH[13]=-0.000009 
XSFTQ[1]=0.000000 
XSFTQ[2]=77.645714 
XSFTQ[3]=150.000000 
XSFTQ[4]=212.132034 
XSFTQ[5]=259.807617 
XSFTQ[6]=289.777771 
XSFTQ[7]=300.000000 
XSFTQ[8]=289.777740 
XSFTQ[9]=259.807587 
XSFTQ[10]=212.132034  
XSFTQ[11]=150.000015 
XSFTQ[12]=77.645676 
XSFTQ[13]=-0.000026

case.add_variable('sfc_sens_flx',XSFTH[1:],time=XTIMEF[1:])
case.add_variable('sfc_lat_flx', XSFTQ[1:],time=XTIMEF[1:])

# Potential temperature advection
zadv = [0.,        1500., 3000., 5000.]
thadv = [-0.0000347, 0., 0., 0.]
rvadv = [0., 0., 0., 0.]
timeF = [0.,43200.]
nzf = len(zadv)

thadv = np.zeros((ntf,nzf),dtype=np.float64)
rvadv = np.zeros((ntf,nzf),dtype=np.float64)
for it in range(0,ntf):
    thadv[it,:] = [-0.0000347, 0., 0., 0.]
    rvadv[it,:] = [ 0., 0., 0., 0.]

case.add_variable('theta_adv',thadv,time=timeF,lev=zadv,levtype='altitude')
case.add_variable('rv_adv',rvadv,time=timeF,lev=zadv,levtype='altitude')


################################################
# 4. Attributes
################################################

# advection of theta and rv
case.set_attribute("adv_theta",1)
case.set_attribute("adv_rv",1)
# potential temperature radiative tendency is included in advection
case.set_attribute("rad_theta","adv")
# Geostrophic wind forcing
case.set_attribute("forc_geo",1)
# Surface flux forcing, wind stress is computed using z0
case.set_attribute("surfaceType","land")
case.set_attribute("surfaceForcing","surfaceFlux")
case.set_attribute("surfaceForcingWind","z0")

################################################
# 5. Writing file
################################################

case.write('SCMS_MESONH_DEF_driver.nc',verbose=False)

if lverbose:
    case.info()

################################################
# 6. Ploting, if asked
################################################

if lplot:
    case.plot(rep_images='./images/driver_DEF/',timeunits='hours')
