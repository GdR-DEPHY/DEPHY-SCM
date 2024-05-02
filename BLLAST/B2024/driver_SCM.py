#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 6 March 2024

@author: Romain Roehrig 

Modifications:

"""

## BLLAST/B2024 SCM-enabled case definition

import netCDF4 as nc
import numpy as np

from dephycf.Case import Case
from dephycf import thermo

################################################
# 0. General configuration of the present script
################################################

lplot = True     # plot the new version of the case
lcompare = True  # plot comparisons between original and new versions
lverbose = False # print information on variables and case

################################################
# 1. Get the original version of the case
################################################

# initialize the case structure for the original version
case = Case('BLLAST/B2024')

# read case information in file
case.read('BLLAST_B2024_DEF_driver.nc')

# display some information about the case
if lverbose:
    case.info()

################################################
# 2. Interpolate onto a new grid, same for all the variables
#    and add new variables if needed
################################################

# Extend profiles constantly towards the surface
#case.extend_init_wind(height=0)
#case.extend_init_theta(height=0)
#case.extend_init_rv(height=0)

# Extend profiles above what 17500 m
hera5 = 18000 # continue with ERA5 above hera5 (in m)
with nc.Dataset('../aux/ERA5/ERA5_P2OA_20110620000000-20110620230000.nc') as f:
    temp = np.average(np.squeeze(f['ta'][:,::-1]), axis=0)
    pa = f['plev'][::-1]
    theta = thermo.t2theta(p=pa, temp=temp)
    height0 = np.average(np.squeeze(f['zg'][:,::-1]), axis=0)
    theta = theta[height0 >= hera5]
    height = height0[height0 >= hera5]
    case.extend_init_theta(theta=theta, height=height)

    hur = np.average(np.squeeze(f['qv'][:,::-1]), axis=0)*0
    hur = hur[height0 >= hera5]
    case.extend_init_hur(hur=hur, height=height)

    u = np.average(np.squeeze(f['ua'][:,::-1]), axis=0)
    u = u[height0 >= hera5]
    v = np.average(np.squeeze(f['va'][:,::-1]), axis=0)
    v = v[height0 >= hera5]
    case.extend_init_wind(u=u, v=v, height=height)

# add a surface temperature
with nc.Dataset('../aux/ERA5/ERA5_P2OA_20110620000000-20110620230000.nc') as f:
    dates = nc.num2date(f['time'][:], units=f['time'].units)
    index = [d.day == 20 and d.hour >= 5 for d in dates ]
    ts = np.squeeze(f['tskin'][index])
    times = f['time'][index]
    case.add_surface_skin_temp(ts, time=times-times[0])

# Grids onto which interpolate the input data

# New vertical grid, 10-m resolution from surface to 3000 m and 100-m resolution above, up to htop
htop = 50000
levout = np.array(list(range(0,3000,10)) + list(range(3100,int(htop)+1,100)),dtype=np.float64)

#  New temporal grid, from 05:00 UTC to 18:00 UTC, 20 June 2011, 30-min timestep
timeout = np.array(range(0,46800+1,1800),dtype=np.float64) 

# conversion
newcase = case.convert2SCM(time=timeout,lev=levout,levtype='altitude')

# update some attributes
newcase.set_title("Forcing and initial conditions for BLLAST/B2024 case - SCM-enabled version")
newcase.set_script("DEPHY-SCM/BLLAST/B2024/driver_SCM.py")

# display some information about the new version of the case
if lverbose:
    newcase.info()

################################################
# 3. Save new version of the case in netcdf file
################################################

# save the new version of the case in netcdf file 
newcase.write('BLLAST_B2024_SCM_driver.nc')

################################################
# 4. Plots if asked
################################################

if lplot:
    newcase.plot(rep_images='./images/driver_SCM/',timeunits='hours')

if lcompare:
    newcase.plot_compare(case,rep_images='./images/compare/',label1="SCM-enabled",label2="Original",timeunits='hours')
