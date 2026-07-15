#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 9 November 2022

@author: Catherine Rio

Modification
  2022/11/08, C. Rio: EUROCS case
  2024/02/05, F. Couvreux: adaptation au dernier format
  2026/07/15, N. Villefranque: clean for publication
"""

## EUROCS SCM-enabled case definition

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
case = Case('EUROCS/REF')

# read case information in file
case.read('EUROCS_REF_DEF_driver.nc')

fin = nc.Dataset('EUROCS_REF_DEF_driver.nc','r')

# display some information about the case
if lverbose:
    case.info()

################################################
# 2. Interpolate onto a new grid, same for all the variables
#    and add new variables if needed
################################################

# Grid onto which interpolate the input data

# Same grid as in driver DEF
levout = np.array(fin['lev_theta'][:])

# Same grid as in driver DEF
timeout = np.array(fin['time_hfss'][:])

# Conversion
newcase = case.convert2SCM(time=timeout,lev=levout,levtype='pressure')

# Update some attributes
newcase.set_title("Forcing and initial conditions for EUROCS case - SCM-enabled version")
newcase.set_script("DEPHY-SCM/EUROCS/REF/driver_SCM.py")

# display some information about the new version of the case
if lverbose:
    newcase.info()

################################################
# 3. Save new version of the case in netcdf file
################################################

# save the new version of the case in netcdf file 
newcase.write('EUROCS_REF_SCM_driver.nc')

################################################
# 4. Plots if asked
################################################

if lplot:
    newcase.plot(rep_images='./images/driver_SCM/',timeunits='hours')

if lcompare:
    newcase.plot_compare(case,rep_images='./images/compare/',label1="SCM-enabled",label2="Original",timeunits='hours')
