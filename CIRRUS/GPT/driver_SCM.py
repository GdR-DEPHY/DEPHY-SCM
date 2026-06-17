#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 29 November 20199 June 2026

@author: Benoît Vié
"""

## Cirru case definition

import netCDF4 as nc
import numpy as np

from dephycf.Case import Case
from dephycf import thermo

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
case = Case('CIRRUS/GPT')

# read case information in file
case.read('CIRRUS_GPT_DEF_driver.nc')

# display some information about the case
if lverbose:
    case.info()

################################################
# 2. Interpolate onto a new grid, same for all the variables
#    and add new variables if needed
################################################

# Extend profiles 
# Wind: just constant above what is presently defined
case.extend_init_wind(height=20000.)

# Total water: 0 above what is presently defined
case.extend_init_rt(rt=[0,0], height=[12000., 20000.])

# Temperature
case.extend_init_theta(height=20000)

# Add a surface temperature based on SGP observations
case.add_surface_skin_temp(290.)

# grid onto which interpolate the input data

# New vertical grid, 10-m resolution from surface to 20000 m (above the surface)
levout = np.array(range(0,20001,10),dtype=np.float64) 

# New temporal grid, from 00 UTC, 12:00 UTC, 30-min timestep
timeout = np.array(range(0,43200+1,1800),dtype=np.float64) 

# conversion
newcase = case.convert2SCM(time=timeout,lev=levout,levtype='altitude')

# update some attributes
newcase.set_title("Forcing and initial conditions for Cirrus case - SCM-enabled version")
newcase.set_script("DEPHY-SCM/CIRRUS/GPT/driver_SCM.py")

# display some information about the new version of the case
if lverbose:
    newcase.info()

################################################
# 3. Save new version of the case in netcdf file
################################################

# save the new version of the case in netcdf file 
newcase.write('CIRRUS_GPT_SCM_driver.nc')

################################################
# 4. Plots if asked
################################################

if lplot:
    newcase.plot(rep_images='./images/driver_SCM/',timeunits='hours')

if lcompare:
    newcase.plot_compare(case,rep_images='./images/compare/',label1="SCM-enabled",label2="Original",timeunits='hours')
