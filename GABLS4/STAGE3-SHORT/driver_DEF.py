#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on 09 June 2020

@author: Romain Roehrig
"""

## GABLS4/STAGE3-SHORT original case definition

import sys
sys.path = ['../../utils/',] + sys.path

import netCDF4 as nc
import numpy as np

import thermo

from Case import Case

################################################
# 0. General configuration of the present script
################################################

lplot = True # plot all the variables
lverbose = False # print information about variables and case

################################################
# 1. General information about the case
################################################

case_name = 'GABLS4'
subcase_name = 'STAGE3-SHORT'

case = Case('{0}/{1}'.format(case_name,subcase_name),
        lat=-75.1,
        lon=-123.3,
        startDate="20091211100000",
        endDate="20091211220000",
        zorog=3233.)

case.set_title("Forcing and initial conditions for the {0}/{1} case - Original definition".format(case_name,subcase_name))
case.set_reference("Couvreux et al. (2020, BLM); Input file GABLS4_SCM_LES_STAGE3_10h.nc downloaded at http://www.umr-cnrm.fr/aladin/meshtml/GABLS4/description_October2016.html")
case.set_author("R. Roehrig")
case.set_script("DEPHY-SCM/{0}/{1}/driver_DEF.py".format(case_name,subcase_name))
case.set_comment("Use of file GABLS4_SCM_LES_STAGE3_10h.nc downloaded at  http://www.umr-cnrm.fr/aladin/meshtml/GABLS4/description_October2016.html")

################################################
# 2. Input netCDF file
################################################

fin = nc.Dataset('GABLS4_SCM_LES_STAGE3_10h.nc','r')

################################################
# 2. Initial state
################################################

# Surface pressure
ps = fin['psurf'][:]
case.add_variable('ps',ps)

nlev, = fin['height'].shape

# Height
height = np.zeros(nlev+1,dtype=np.float)
height[1:] = fin['height'][::-1] # reverse altitude order
case.add_variable('height',height,lev=height,levtype='altitude')

# Pressure
pressure = np.zeros(nlev+1,dtype=np.float)
pressure[0] = ps # add surface level value 
pressure[1:]  = fin['pf'][::-1]
case.add_variable('pressure',pressure,lev=height,levtype='altitude')

# Zonal wind (same value for the first two levels)
u = np.zeros(nlev+1,dtype=np.float)
u[1:]  = fin['u'][::-1]
u[0] = u[1]
case.add_variable('u',u,lev=height,levtype='altitude')

# Meridional wind (same value for the first two levels)
v = np.zeros(nlev+1,dtype=np.float)
v[1:] = fin['v'][::-1]
v[0] = v[1] 
case.add_variable('v',v,lev=height,levtype='altitude')

# Potential temperature
theta = np.zeros(nlev+1,dtype=np.float)
theta[0] = thermo.t2theta(p=ps,temp=fin['Tg'][0])
theta[1:] = fin['theta'][::-1]
case.add_variable('theta',theta,lev=height,levtype='altitude')

# Temperature
temp = np.zeros(nlev+1,dtype=np.float)
temp[0] = fin['Tg'][0]
temp[1:] = fin['t'][::-1]
case.add_variable('temp',temp,lev=height,levtype='altitude')

# Total water content (0 everywhere)
qv = np.zeros(nlev+1,dtype=np.float)
qv[1:] = fin['qv'][::-1]
case.add_variable('qv',qv,lev=height,levtype='altitude')

################################################
# 3. Forcing
################################################

timeForc = fin['time'][0:12]
nt, = timeForc.shape

# Constant geostrophic wind (same value for the first two levels)
ug = np.zeros((nt,nlev+1),dtype=np.float)
vg = np.zeros((nt,nlev+1),dtype=np.float)
ug[:,1:] = fin['Ug'][0:12,::-1]
vg[:,1:] = fin['Vg'][0:12,::-1]
ug[:,0] = ug[:,1]
vg[:,0] = vg[:,1]
case.add_variable('ug',ug,time=timeForc,lev=height,levtype='altitude')
case.add_variable('vg',vg,time=timeForc,lev=height,levtype='altitude')

# Surface temperature
ts = fin['Tg'][0:12]
case.add_variable('ts',ts,time=timeForc)

################################################
# 4. Attributes
################################################

# Radiation should be activated
case.set_attribute("rad_temp",'adv')

# Geostrophic wind 
case.set_attribute("forc_geo",1)

# Surface flux forcing, wind stress is computed using z0
case.set_attribute("surfaceType","land")
case.set_attribute("surfaceForcing","ts")
case.set_attribute("surfaceForcingWind","z0")
case.set_attribute("z0",0.001) # The last version in Couvreux et al. (2020) recommend 0.001 

################################################
# 5. Writing file
################################################

case.write('{0}_{1}_DEF_driver.nc'.format(case_name,subcase_name),verbose=False)

fin.close()

if lverbose:
    case.info()

################################################
# 6. Ploting, if asked
################################################

if lplot:
    case.plot(rep_images='./images/driver_DEF/',timeunits='hours')
