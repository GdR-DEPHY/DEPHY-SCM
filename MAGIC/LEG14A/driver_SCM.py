#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on 21 January 2022

@author: Maike Ahlgrimm

Modification

"""

## MAGIC Leg14A SCM-enabled case definition

import sys
sys.path = ['../../utils/',] + sys.path

import netCDF4 as nc
import numpy as np

from datetime import datetime, timedelta

from Case import Case
import thermo

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
case = Case('MAGIC/LEG14A')

# read case information in file
case.read('MAGIC_LEG14A_DEF_driver.nc')

# display some information about the case
if lverbose:
    case.info()

################################################
# 2. Interpolate onto a new grid, same for all the variables
#    and add new variables if needed
################################################

# New vertical grid, 10-m resolution from surface to 6000 m (above the surface)
levout = np.array(range(0,20001,10),dtype=np.float64) 

# New temporal grid, from 5:30 UTC, 8 July 2013 to 5:30 UTC, 12 August 2013, 30-min timestep
timeout = np.array(range(0,345600,1800),dtype=np.float64)

# conversion
newcase = case.convert2SCM(time=timeout,lev=levout,usetemp=True,levtype='altitude')

# update some attributes
newcase.set_title("Forcing and initial conditions for MAGIC Leg14A case - SCM-enabled version")
newcase.set_script("DEPHY-SCM/MAGIC/LEG14A/driver_SCM.py")

# display some information about the new version of the case
if lverbose:
    newcase.info()

################################################
# 3. Save new version of the case in netcdf file
################################################

# save the new version of the case in netcdf file 
newcase.write('MAGIC_LEG14A_SCM_driver.nc')

################################################
# 4. Plots if asked
################################################

if lplot:
    newcase.plot(rep_images='./images/driver_SCM/',timeunits='hours')

if lcompare:
    newcase.plot_compare(case,rep_images='./images/compare/',label1="SCM-enabled",label2="Original",timeunits='hours')
