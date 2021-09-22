#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on 21 September 2021

@author: Romain Roehrig
"""

## FIRE Reference case SCM-enabled definition

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
case = Case('FIRE/REF')

# read case information in file
case.read('FIRE_REF_DEF_driver.nc')

# display some information about the case
if lverbose:
    case.info()

################################################
# 2. Interpolate onto a new grid, same for all the variables
#    and add new variables if needed
################################################

# Extend profiles
# Wind: just make it constant above what is presently defined
case.extend_init_wind(height=80000.)

# Total water content: linear to 0 between 1200 and 1300m, 0 above
case.extend_init_qt(qt=[0,0], height=[1300,80000.])

# Geostrophic wind: just make it constant above what is presently defined
case.extend_geostrophic_wind(height=80000.)

# Thetal and qt large-scale advection: linear to 0 between 1200 and 1300m, 0 above
case.extend_thetal_advection(thetal_adv=[0,0], height=[1300.,80000.])
case.extend_qt_advection(qt_adv=[0,0], height=[1300.,80000.])

# grid onto which interpolate the input data

# New vertical grid, 10-m resolution from surface to 60000 m (above the surface)
levout = np.array(list(range(0,2001,10)) + list(range(2100,80001,100)),dtype=np.float64) 

# New temporal grid, from 00:00 UTC, 16 December 2004 to 00:00 UTC 19 December 2004, 1-hour timestep
timeout = np.arange(0.,(37+1)*3600.,3600.)

# conversion
newcase = case.convert2SCM(time=timeout,lev=levout,levtype='altitude')

# update some attributes
newcase.set_title("Forcing and initial conditions for FIRE Reference case - SCM-enabled version")
newcase.set_script("DEPHY-SCM/FIRE/REF/driver_SCM.py")

# display some information about the new version of the case
if lverbose:
    newcase.info()

################################################
# 3. Save new version of the case in netcdf file
################################################

# save the new version of the case in netcdf file 
newcase.write('FIRE_REF_SCM_driver.nc')

################################################
# 4. Plots if asked
################################################

if lplot:
    newcase.plot(rep_images='./images/driver_SCM/',timeunits='hours')

if lcompare:
    newcase.plot_compare(case,rep_images='./images/compare/',label1="SCM-enabled",label2="Original",timeunits='hours')
