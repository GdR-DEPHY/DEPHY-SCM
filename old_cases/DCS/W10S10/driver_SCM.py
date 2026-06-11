#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 30 May 2024

@author: Gaston Bidoux, Romain Roehrig, Catherine Rio

Modification
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
case = Case('DCS/W10S10')

# read case information in file
case.read('DCS_W10S10_DEF_driver.nc')

fin = nc.Dataset('DCS_W10S10_DEF_driver.nc','r')

# display some information about the case
if lverbose:
    case.info()

################################################
# 2. Interpolate onto a new grid, same for all the variables
#    and add new variables if needed
################################################

levout = np.arange(0,20001,5)
timeout = [0.,36000.]

newcase = case.convert2SCM(time=timeout, lev=levout, levtype='altitude')

# update some attributes
newcase.set_title("Forcing and initial conditions for DCS/W10S10 case - SCM-enabled version")
newcase.set_script("DEPHY-SCM/DCS/W10S10/driver_SCM.py")

# display some information about the new version of the case
if lverbose:
    newcase.info()

################################################
# 3. Save new version of the case in netcdf file
################################################

# save the new version of the case in netcdf file 
newcase.write('DCS_W10S10_SCM_driver.nc')

################################################
# 4. Plots if asked
################################################

if lplot:
    newcase.plot(rep_images='./images/driver_SCM/',timeunits='hours')

if lcompare:
    newcase.plot_compare(case,rep_images='./images/compare/',label1="SCM-enabled",label2="Original",timeunits='hours')
