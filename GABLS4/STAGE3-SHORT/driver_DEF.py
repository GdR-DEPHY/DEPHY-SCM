#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on 09 June 2020

@author: Romain Roehrig

Modification
  2020/11/12, E. Vignon:  bugfix for ps + land->landice + no evaporation option
  2021/01/03, R. Roehrig: update for improved case definition interface.
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
        lon=123.3,
        startDate="20091211100000",
        endDate="20091211220000",
        surfaceType="land",
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
case.add_init_ps(ps)

nlev, = fin['height'].shape

# Height
height = np.zeros(nlev+1,dtype=np.float)
height[1:] = fin['height'][::-1] # reverse altitude order
case.add_init_height(height,lev=height,levtype='altitude')

# Pressure
pressure = np.zeros(nlev+1,dtype=np.float)
pressure[0] = ps # add surface level value 
pressure[1:]  = fin['pf'][::-1]
case.add_init_pressure(pressure,lev=height,levtype='altitude')

# Zonal and meridional wind (same value for the first two levels)
u = np.zeros(nlev+1,dtype=np.float)
u[1:]  = fin['u'][::-1]
u[0] = u[1]

v = np.zeros(nlev+1,dtype=np.float)
v[1:] = fin['v'][::-1]
v[0] = v[1] 

case.add_init_wind(u=u,v=v, lev=height, levtype='altitude')

# Potential temperature
theta = np.zeros(nlev+1,dtype=np.float)
theta[0] = thermo.t2theta(p=ps,temp=fin['Tg'][0])
theta[1:] = fin['theta'][::-1]

case.add_init_theta(theta, lev=height, levtype='altitude')

# Temperature
temp = np.zeros(nlev+1,dtype=np.float)
temp[0] = fin['Tg'][0]
temp[1:] = fin['t'][::-1]

case.add_init_temp(temp, lev=height, levtype='altitude')

# Total water content (0 everywhere)
qv = np.zeros(nlev+1,dtype=np.float)
qv[1:] = fin['qv'][::-1]

case.add_init_qv(qv, lev=height, levtype='altitude')

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

case.add_geostrophic_wind(ug=ug,vg=vg,time=timeForc,lev=height,levtype='altitude')

# Surface temperature
ts = fin['Tg'][0:12]

case.add_forcing_ts(ts,time=timeForc,z0=0.001) # The last version in Couvreux et al. (2020) recommend 0.001

# No surface evaporation (in fact no moisture at all)
case.deactivate_surface_evaporation()

# No radiation
case.deactivate_radiation()

################################################
# 4. Writing file
################################################

case.write('{0}_{1}_DEF_driver.nc'.format(case_name,subcase_name),verbose=False)

fin.close()

if lverbose:
    case.info()

################################################
# 5. Ploting, if asked
################################################

if lplot:
    case.plot(rep_images='./images/driver_DEF/',timeunits='hours')
