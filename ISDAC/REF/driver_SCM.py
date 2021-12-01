#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on 05 January 2021

@author: Etienne Vignon
"""

## ISDAC SCM-enabled case definition

import sys
sys.path = ['../../utils/',] + sys.path

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
case = Case('ISDAC/REF')

# read case information in file
case.read('ISDAC_REF_DEF_driver.nc')

# display some information about the case
if lverbose:
    case.info()

################################################
# 2. Interpolate onto a new grid, same for all the variables
#    and add new variables if needed
################################################

# grid onto which interpolate the input data

# New vertical grid, 10-m resolution from surface to 5000 m (above the surface)
levout=np.array(np.linspace(0,5000,501),dtype=np.float64)

# New temporal grid, from 17:00 UTC, 09 October 2004 to 05:00 UTC, 10 October 2004, 1-hour timestep
timeout = np.array([0,28800],dtype=np.float64) 

# conversion
newcase = case.convert2SCM(time=timeout,lev=levout,levtype='altitude')


# update some attributes
newcase.set_title("Forcing and initial conditions for ISDAC/REF case - SCM-enabled version")
newcase.set_script("driver_SCM.py")

# display some information about the new version of the case
if lverbose:
    newcase.info()

################################################
# 3. Save new version of the case in netcdf file
################################################

# save the new version of the case in netcdf file 
newcase.write('ISDAC_REF_SCM_driver.nc')

################################################
# 4. Plots if asked
################################################

if lplot:
    newcase.plot(rep_images='./images/driver_SCM/',timeunits='hours')

if lcompare:
    newcase.plot_compare(case,rep_images='./images/compare/',label1="SCM-enabled",label2="Original",timeunits='hours')
