#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 10 Avril 2020

@author: Fleur Couvreux
"""

## AYOTTE/00SC

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
case = Case('AYOTTE/00SC')

# read case information in file
case.read('AYOTTE_00SC_DEF_driver.nc')

# display some information about the case
if lverbose:
    case.info()

################################################
# 2. Interpolate onto a new grid, same for all the variables
#    and add new variables if needed
################################################

# grid onto which interpolate the input data

# New vertical grid, 10-m resolution from surface to 6000 m (above the surface)
levout = np.array(range(0,6001,10),dtype=np.float64) 

# New temporal grid, from 10:00 UTC to 17:00 UTC, 11 December 2009, 30-min timestep
timeout = np.array(range(0,25201,1800),dtype=np.float64) 

# conversion
newcase = case.convert2SCM(time=timeout,lev=levout,levtype='altitude')

# add a surface temperature. To be improved...
ts = timeout*0. + 310. # same shape as timeout

newcase.add_variable('ts',ts,time=timeout,timeid='time')

# update some attributes
newcase.set_title("Forcing and initial conditions for AYOTTE-00SC case - SCM-enabled version")
newcase.set_script("driver_SCM.py")

# display some information about the new version of the case
if lverbose:
    newcase.info()

################################################
# 3. Save new version of the case in netcdf file
################################################

# save the new version of the case in netcdf file 
newcase.write('AYOTTE_00SC_SCM_driver.nc')

################################################
# 4. Plots if asked
################################################

if lplot:
    newcase.plot(rep_images='./images/driver_SCM/',timeunits='hours')

if lcompare:
    newcase.plot_compare(case,rep_images='./images/compare/',label1="SCM-enabled",label2="Original",timeunits='hours')
