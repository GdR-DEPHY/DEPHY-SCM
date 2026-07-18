#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 10 December 2019

@author: Romain Roehrig
"""

## ARM-Cumulus SCM-enabled case definition

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
case = Case('IHOP/REF')

# read case information in file
case.read('IHOP_REF_DEF_driver.nc')

# display some information about the case
if lverbose:
    case.info()

################################################
# 2. Interpolate onto a new grid, same for all the variables
#    and add new variables if needed
################################################

# grid onto which interpolate the input data

# New vertical grid, 10-m resolution from surface to 4975 m (above the surface)
levout = np.array(range(0,4976,10), dtype=float)

# New temporal grid, 7h with 30min time step
timeout = np.array(range(0,7*3600+1,1800), dtype=float)

# conversion
newcase = case.convert2SCM(time=timeout,lev=levout,levtype='altitude')

# add a surface temperature. To be improved...
ts = timeout*0. + 300 # same shape as timeout

newcase.add_variable('ts',ts,time=timeout,timeid='time')

# update some attributes
newcase.set_title("Forcing and initial conditions for IHOP case - SCM-enabled version")
newcase.set_script("DEPHY-SCM/IHOP/driver_SCM.py")

# display some information about the new version of the case
if lverbose:
    newcase.info()

################################################
# 3. Save new version of the case in netcdf file
################################################

# save the new version of the case in netcdf file 
newcase.write('IHOP_REF_SCM_driver.nc')

################################################
# 4. Plots if asked
################################################

if lplot:
    newcase.plot(rep_images='./images/driver_SCM/',timeunits='hours')

if lcompare:
    newcase.plot_compare(case,rep_images='./images/compare/',label1="SCM-enabled",label2="Original",timeunits='hours')
