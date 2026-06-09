#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 09 June 2026

@author: DEPHY team

Modification
"""

## Cirrus case (idealized) - Borella 2025

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
case = Case('CIRRUS/W01')

# read case information in file
case.read('CIRRUS_W01_DEF_driver.nc')

# display some information about the case
if lverbose:
    case.info()

################################################
# 2. Interpolate onto a new grid, same for all the variables
#    and add new variables if needed
################################################

# grid onto which interpolate the input data

# New vertical grid, 10-m resolution from surface to 12000 m (above the surface)
levout = np.array(range(0,12000,10),dtype=np.float64) 

# New temporal grid, 12h with 5-min timestep
timeout = np.arange(0.,12*3600., 300.)

# conversion
newcase = case.convert2SCM(time=timeout,lev=levout,levtype='altitude')

# update some attributes
newcase.set_title("Forcing and initial conditions for CIRRUS case - SCM-enabled original version")
newcase.set_script("DEPHY-SCM/CIRRUS/W01/driver_SCM.py")

# display some information about the new version of the case
if lverbose:
    newcase.info()

################################################
# 3. Save new version of the case in netcdf file
################################################

# save the new version of the case in netcdf file 
newcase.write('CIRRUS_W01_SCM_driver.nc')

################################################
# 4. Plots if asked
################################################

if lplot:
    newcase.plot(rep_images='./images/driver_SCM/',timeunits='hours')

if lcompare:
    newcase.plot_compare(case,rep_images='./images/compare/',label1="SCM-enabled",label2="Original",timeunits='hours')
