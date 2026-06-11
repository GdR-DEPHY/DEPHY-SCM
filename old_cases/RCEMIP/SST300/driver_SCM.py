#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 2 December 2024

@author: Romain Roehrig
"""

## RCEMIP case definition - SST = 300 K - SCM-enabled definition

SST = 300

import numpy as np
from datetime import datetime, timedelta

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
case = Case(f'RCEMIP/SST{SST}')

# read case information in file
case.read(f'RCEMIP_SST{SST}_DEF_driver.nc')

# display some information about the case
if lverbose:
    case.info()

################################################
# 2. Interpolate onto a new grid, same for all the variables
#    and add new variables if needed
################################################

# grid onto which interpolate the input data

levout = np.arange(0,80001,10,dtype=np.float64) 

# New temporal grid, from 00:00 UTC, 16 December 2004 to 00:00 UTC 19 December 2004, 1-hour timestep
timeout = [0,100*86400]

# conversion
newcase = case.convert2SCM(time=timeout,lev=levout,levtype='altitude')

# update some attributes
newcase.set_title("Forcing and initial conditions for Radiative-Convection Equilibrium MIP (RCEMIP) case - SST=300K - SCM-enabled version")
newcase.set_script(f"DEPHY-SCM/RCEMIP/SST{SST}/driver_SCM.py")

# display some information about the new version of the case
if lverbose:
    newcase.info()

################################################
# 3. Save new version of the case in netcdf file
################################################

# save the new version of the case in netcdf file 
newcase.write(f'RCEMIP_SST{SST}_SCM_driver.nc')

################################################
# 4. Plots if asked
################################################

if lplot:
    newcase.plot(rep_images='./images/driver_SCM/',timeunits='hours')

if lcompare:
    newcase.plot_compare(case,rep_images='./images/compare/',label1="SCM-enabled",label2="Original",timeunits='hours')
