#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on 29 November 2019

@author: Romain Roehrig

Modification
  2020/11/11, R. Roehrig: update for improved case definition interface.
"""

## ARM-Cumulus SCM-enabled case definition

import sys
sys.path = ['../../utils/',] + sys.path

import netCDF4 as nc
import numpy as np

from Case import Case

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
case = Case('AMMA/REF')

# read case information in file
case.read('AMMA_REF_DEF_driver.nc')

fin = nc.Dataset('AMMA_REF_DEF_driver.nc','r')

# display some information about the case
if lverbose:
    case.info()

################################################
# 2. Interpolate onto a new grid, same for all the variables
#    and add new variables if needed
################################################

# grid onto which interpolate the input data

# New vertical grid, 10-m resolution from surface to 6000 m (above the surface)
#levout = np.array(range(0,6001,10),dtype=np.float64) 
levout = np.array(fin['lev_zh'][:])
# New temporal grid, from 11:30 UTC, 21 June 1997 to 02:00 UTC, 22 June 1997, 30-min timestep
#timeout = np.array(range(0,86400+2*3600+1-41400,1800),dtype=np.float64) 
timeout = np.array(fin['time_hfss'][:])
# conversion
newcase = case.convert2SCM(time=timeout,lev=levout,levtype='altitude')

# add a surface temperature. To be improved...
#ts = timeout*0. + 300 # same shape as timeout

#newcase.add_surface_temp(ts,time=timeout,timeid='time')

# update some attributes
newcase.set_title("Forcing and initial conditions for AMMA-Cumulus case - SCM-enabled version")
newcase.set_script("DEPHY-SCM/AMMA/REF/driver_SCM.py")

# display some information about the new version of the case
if lverbose:
    newcase.info()

################################################
# 3. Save new version of the case in netcdf file
################################################

# save the new version of the case in netcdf file 
newcase.write('AMMA_REF_SCM_driver.nc',verbose=False)

################################################
# 4. Plots if asked
################################################

if lplot:
    newcase.plot(rep_images='./images/driver_SCM/',timeunits='hours')

if lcompare:
    newcase.plot_compare(case,rep_images='./images/compare/',label1="SCM-enabled",label2="Original",timeunits='hours')
