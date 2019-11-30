#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on 29 November 2019

@author: Romain Roehrig
"""

## ARM-Cumulus SCM-enabled case definition

import sys
sys.path = ['../../utils/',] + sys.path

import numpy as np

import SCM_utils as utils

################################################
# 0. General configuration of the present script
################################################

lplot = True    # plot the new version of the case
lcompare = True # plot comparisons between original and new versions

################################################
# 1. Get the original version of the case
################################################

# initialize the case structure for the original version
case = utils.Case('RICO/SHORT')

# read case information in file
case.read('RICO_SHORT_orig.nc')

# display some information about the case
case.info()

################################################
# 2. Interpolate onto a new grid, same for all the variables
#    and add new variables if needed
################################################

# grid onto which interpolate the input data

# New vertical grid, 10-m resolution from surface to 60000 m (above the surface)
levout = np.array(range(0,60001,10),dtype=np.float64) 

# New temporal grid, from 00:00 UTC, 16 December 2004 to 00:00 UTC 19 December 2004, 1-hour timestep
timeout = np.arange(0.,(72+1)*3600.,3600.)

# conversion
newcase = case.convert2SCM(time=timeout,lev=levout,levtype='altitude')

# update some attributes
newcase.set_comment("Forcing and initial conditions for RICO composite short case - SCM-enabled version")
newcase.set_script("DEPHY-SCM/RICO/SHORT/convert_SCM.py")

# display some information about the new version of the case
newcase.info()

################################################
# 3. Save new version of the case in netcdf file
################################################

# save the new version of the case in netcdf file 
newcase.write('RICO_SHORT_1D.nc',verbose=False)

################################################
# 4. Plots if asked
################################################

if lplot:
    newcase.plot(rep_images='./images/setup_1D/',timeunits='hours')

if lcompare:
    newcase.plot_compare(case,rep_images='./images/compare/',label1="SCM-enabled",label2="Original",timeunits='hours')
