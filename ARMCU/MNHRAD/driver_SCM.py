#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 10 December 2019

@author: Romain Roehrig

Modification
  2020/11/11, R. Roehrig: update for improved case definition interface.
  2021/07/05, R. Roehrig: update tskin
"""

## ARM-Cumulus SCM-enabled case definition

import netCDF4 as nc
import numpy as np

from dephycf.Case import Case

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
case = Case('ARMCU/MNHRAD')

# read case information in file
case.read('ARMCU_MNHRAD_DEF_driver.nc')

# display some information about the case
if lverbose:
    case.info()

################################################
# 2. Interpolate onto a new grid, same for all the variables
#    and add new variables if needed
################################################

# add a surface temperature. To be improved...
with nc.Dataset('../aux/tskin/tskin_SGP_C1_irt10m_19970621003000-19970622233000.nc') as f:
    dates = nc.num2date(f['time'][:], units=f['time'].units)
    index = [d.day == 21 and d.hour >= 11 or d.day == 22 and d.hour <= 2 for d in dates ]
    times = f['time'][index]
    tskin = f['tskin'][index]

    case.add_surface_skin_temp(tskin,time=times-times[0])

# grid onto which interpolate the input data

# New vertical grid, 10-m resolution from surface to 6000 m (above the surface)
levout = np.array(range(0,6001,10),dtype=np.float64) 

# New temporal grid, from 11:30 UTC, 21 June 1997 to 02:30 UTC, 22 June 1997, 30-min timestep
timeout = np.array(range(0,86400+2*3600+1800+1-41400,1800),dtype=np.float64) 

# conversion
newcase = case.convert2SCM(time=timeout,lev=levout,levtype='altitude')

# update some attributes
newcase.set_title("Forcing and initial conditions for ARM-Cumulus case - SCM-enabled Meso-NH version")
newcase.set_script("DEPHY-SCM/ARMCU/MNHRAD/driver_SCM.py")

# display some information about the new version of the case
if lverbose:
    newcase.info()

################################################
# 3. Save new version of the case in netcdf file
################################################

# save the new version of the case in netcdf file 
newcase.write('ARMCU_MNHRAD_SCM_driver.nc')

################################################
# 4. Plots if asked
################################################

if lplot:
    newcase.plot(rep_images='./images/driver_SCM/',timeunits='hours')

if lcompare:
    newcase.plot_compare(case,rep_images='./images/compare/',label1="SCM-enabled",label2="Original",timeunits='hours')
