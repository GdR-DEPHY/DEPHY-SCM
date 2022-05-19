#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on 29 November 2019

@author: Romain Roehrig
"""

## Bomex SCM-enabled case definition

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
case = Case('KB2006/REF')

# read case information in file
case.read('KB2006_REF_DEF_driver.nc')

# display some information about the case
if lverbose:
    case.info()

################################################
# 2. Interpolate onto a new grid, same for all the variables
#    and add new variables if needed
################################################

htop = 20000.

# Extend profiles
# Total water following same gradient as below, until reaching 0
#case.extend_init_qt([0, 0], height=[5500.,htop])

# Same for thetal
#case.extend_init_thetal(308.2+3.65e-3*htop, height=htop)

# grid onto which interpolate the input data

# New vertical grid, 10-m resolution from surface to 6000 m (above the surface)
levout = np.array(list(range(0,3000,10)) + list(range(3100,int(htop)+1,100)),dtype=np.float64) 

# New temporal grid, from 11:30 UTC, 21 June 1997 to 02:00 UTC, 22 June 1997, 30-min timestep
timeout = np.array(range(0,86400*7+3600,3600),dtype=np.float64) 

# conversion
newcase = case.convert2SCM(time=timeout,lev=levout,levtype='altitude')

# update some attributes
newcase.set_title("Forcing and initial conditions for BOMEX case - SCM-enabled version")
newcase.set_script("DEPHY-SCM/BOMEX/REF/driver_SCM.py")

# display some information about the new version of the case
if lverbose:
    newcase.info()

################################################
# 3. Save new version of the case in netcdf file
################################################

# save the new version of the case in netcdf file 
newcase.write('KB2006_REF_SCM_driver.nc')

################################################
# 4. Plots if asked
################################################

if lplot:
    newcase.plot(rep_images='./images/driver_SCM/',timeunits='hours')

if lcompare:
    newcase.plot_compare(case,rep_images='./images/compare/',label1="SCM-enabled",label2="Original",timeunits='hours')
