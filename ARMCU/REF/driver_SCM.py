#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 29 November 2019

@author: Romain Roehrig

Modification
  2020/11/11, R. Roehrig: update for improved case definition interface.
  2021/07/05, R. Roehrig: update tskin
"""

## ARM-Cumulus SCM-enabled case definition

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
case = Case('ARMCU/REF')

# read case information in file
case.read('ARMCU_REF_DEF_driver.nc')

# display some information about the case
if lverbose:
    case.info()

################################################
# 2. Interpolate onto a new grid, same for all the variables
#    and add new variables if needed
################################################

# Extend profiles 
# Wind: just constant above what is presently defined
case.extend_init_wind(height=20000.)

# Total water: 0 above what is presently defined
case.extend_init_rt(rt=[0,0], height=[5600., 20000.])

# Temperature using ERA5
# Only above 12 km to have a nicer and more stable extrapolated profile.
with nc.Dataset('../aux/ERA5/ERA5_SGP_19970621000000-19970622230000.nc') as f:
    temp = np.average(np.squeeze(f['ta'][:,::-1]), axis=0)
    pa = f['plev'][::-1]
    theta = thermo.t2theta(p=pa, temp=temp)
    height = np.average(np.squeeze(f['zg'][:,::-1]), axis=0)
    theta = theta[height >= 12000]
    height = height[height >= 12000]
    case.extend_init_theta(theta=theta, height=height)

# Add a surface temperature based on SGP observations
with nc.Dataset('../aux/tskin/tskin_SGP_C1_irt10m_19970621003000-19970622233000.nc') as f:
    dates = nc.num2date(f['time'][:], units=f['time'].units)
    index = [d.day == 21 and d.hour >= 11 or d.day == 22 and d.hour <= 2 for d in dates ]
    times = f['time'][index]
    tskin = f['tskin'][index]

    case.add_surface_skin_temp(tskin,time=times-times[0])

# grid onto which interpolate the input data

# New vertical grid, 10-m resolution from surface to 6000 m (above the surface)
levout = np.array(range(0,20001,10),dtype=np.float64) 

# New temporal grid, from 11:30 UTC, 21 June 1997 to 02:00 UTC, 22 June 1997, 30-min timestep
timeout = np.array(range(0,86400+2*3600+1-41400,1800),dtype=np.float64) 

# conversion
newcase = case.convert2SCM(time=timeout,lev=levout,levtype='altitude')

# update some attributes
newcase.set_title("Forcing and initial conditions for ARM-Cumulus case - SCM-enabled version")
newcase.set_script("DEPHY-SCM/ARMCU/REF/driver_SCM.py")

# display some information about the new version of the case
if lverbose:
    newcase.info()

################################################
# 3. Save new version of the case in netcdf file
################################################

# save the new version of the case in netcdf file 
newcase.write('ARMCU_REF_SCM_driver.nc')

################################################
# 4. Plots if asked
################################################

if lplot:
    newcase.plot(rep_images='./images/driver_SCM/',timeunits='hours')

if lcompare:
    newcase.plot_compare(case,rep_images='./images/compare/',label1="SCM-enabled",label2="Original",timeunits='hours')
