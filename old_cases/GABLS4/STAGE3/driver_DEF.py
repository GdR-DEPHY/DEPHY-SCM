#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 09 June 2020

@author: Romain Roehrig

Modification
  2020/11/12, E. Vignon:  bugfix for ps + land->landice + no evaporation option
  2021/01/03, R. Roehrig: update for improved case definition interface.
  2024/02/27, R. Roehrig: make the case more consistent with Couvreux et al. (2020) wrt first level (z=0)
"""

## GABLS4/STAGE3 original case definition

import netCDF4 as nc
import numpy as np

from dephycf import thermo
from dephycf.Case import Case

################################################
# 0. General configuration of the present script
################################################

lplot = True # plot all the variables
lverbose = False # print information about variables and case

################################################
# 1. General information about the case
################################################

case_name = 'GABLS4'
subcase_name = 'STAGE3'

case = Case('{0}/{1}'.format(case_name,subcase_name),
        lat=-75.1,
        lon=123.3,
        startDate="20091211000000",
        endDate="20091212120000",
        surfaceType="landice",
        zorog=3233.)

case.set_title("Forcing and initial conditions for the {0}/{1} case - Original definition".format(case_name,subcase_name))
case.set_reference("Couvreux et al. (2020, BLM); Input file GABLS4_SCM_LES_STAGE3.nc downloaded at http://www.umr-cnrm.fr/aladin/meshtml/GABLS4/GABLS4.html")
case.set_author("R. Roehrig")
case.set_script("DEPHY-SCM/{0}/{1}/driver_DEF.py".format(case_name,subcase_name))
case.set_comment("Use of file GABLS4_SCM_LES_STAGE3.nc downloaded at  http://www.umr-cnrm.fr/aladin/meshtml/GABLS4/GABLS4.html")

################################################
# 2. Input netCDF file
################################################

fin = nc.Dataset('GABLS4_SCM_LES_STAGE3.nc','r')

################################################
# 2. Initial state
################################################

# Surface pressure
ps = fin['psurf'][:]
case.add_init_ps(ps)

nlev, = fin['height'].shape

# Height
height = np.zeros(nlev+1,dtype=np.float64)
height[1:] = fin['height'][::-1] # reverse altitude order
case.add_init_height(height,lev=height,levtype='altitude')

# Pressure
pressure = np.zeros(nlev+1,dtype=np.float64)
pressure[0] = ps # add surface level value 
pressure[1:]  = fin['pf'][::-1]
case.add_init_pressure(pressure,lev=height,levtype='altitude')

# Zonal and meridional wind (same value for the first two levels)
u = np.zeros(nlev+1,dtype=np.float64)
u[1:]  = fin['u'][::-1]
u[0] = 0 # from Couvreux et al. (2020), Table 5

v = np.zeros(nlev+1,dtype=np.float64)
v[1:] = fin['v'][::-1]
v[0] = 0 # from Couvreux et al. (2020), Table 5

case.add_init_wind(u=u,v=v, lev=height, levtype='altitude')

# Potential temperature
theta = np.zeros(nlev+1,dtype=np.float64)
#theta[0] = thermo.t2theta(p=ps,temp=fin['Tg'][0])
theta[0] = 271.3 # from Couvreux et al. (2020), Table 5
theta[1:] = fin['theta'][::-1]

case.add_init_theta(theta, lev=height, levtype='altitude')

# Temperature
temp = np.zeros(nlev+1,dtype=np.float64)
temp[0] = fin['Tg'][0]
temp[1:] = fin['t'][::-1]

case.add_init_temp(temp, lev=height, levtype='altitude')

# Total water content (0 everywhere)
qv = np.zeros(nlev+1,dtype=np.float64)
qv[1:] = fin['qv'][::-1]

case.add_init_qv(qv, lev=height, levtype='altitude')

################################################
# 3. Forcing
################################################

timeForc = fin['time'][:]
nt, = timeForc.shape

# Constant geostrophic wind (same value for the first two levels)
ug = np.zeros((nt,nlev+1),dtype=np.float64)
vg = np.zeros((nt,nlev+1),dtype=np.float64)
ug[:,1:] = fin['Ug'][:,::-1]
vg[:,1:] = fin['Vg'][:,::-1]
ug[:,0] = 0 # from Couvreux et al. (2020), Table 5
vg[:,0] = 0 # from Couvreux et al. (2020), Table 5

case.add_geostrophic_wind(ug=ug,vg=vg,time=timeForc,lev=height,levtype='altitude')

# Surface temperature
ts = fin['Tg'][:]

case.add_forcing_ts(ts,time=timeForc,z0=0.001,z0h=1.e-4) # The last version in Couvreux et al. (2020) recommend 0.001

# No surface evaporation (in fact no moisture at all)
case.deactivate_surface_evaporation()

# No radiation
case.deactivate_radiation()

################################################
# 4. Writing file
################################################

case.write('{0}_{1}_DEF_driver.nc'.format(case_name,subcase_name))

fin.close()

if lverbose:
    case.info()

################################################
# 5. Ploting, if asked
################################################

if lplot:
    case.plot(rep_images='./images/driver_DEF/',timeunits='hours')
