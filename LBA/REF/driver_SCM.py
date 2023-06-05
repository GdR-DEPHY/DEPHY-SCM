#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 24 May 2023

@author: F. Couvreux

Modification
  2023/06/05, R. Roehrig: some fixes
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
case = Case('LBA/REF')

# read case information in file
case.read('LBA_REF_DEF_driver.nc')

fin = nc.Dataset('LBA_REF_DEF_driver.nc','r')

# display some information about the case
if lverbose:
    case.info()

################################################
# 2. Interpolate onto a new grid, same for all the variables
#    and add new variables if needed
################################################

# grid onto which interpolate the input data

# New vertical grid, 10-m resolution from surface to 6000 m (above the surface)
levout = np.array(range(0,20001,10),dtype=np.float64) 
# Common temporal grid using that from sensible heat flux
timeout = np.array(fin['time_hfss'][:])
# conversion
newcase = case.convert2SCM(time=timeout,lev=levout,levtype='altitude')

# add a surface temperature. To be improved...
ts = timeout*0. + 310 # same shape as timeout
newcase.add_surface_temp(ts,time=timeout,timeid='time')

# update some attributes
newcase.set_title("Forcing and initial conditions for LBA/REF case (shallow-to-deep transition) - SCM-enabled version")
newcase.set_script("DEPHY-SCM/LBA/REF/driver_SCM.py")

# display some information about the new version of the case
if lverbose:
    newcase.info()

################################################
# 3. Save new version of the case in netcdf file
################################################

# save the new version of the case in netcdf file 
newcase.write('LBA_REF_SCM_driver.nc')

################################################
# 4. Plots if asked
################################################

if lplot:
    newcase.plot(rep_images='./images/driver_SCM/',timeunits='hours')

if lcompare:
    newcase.plot_compare(case,rep_images='./images/compare/',label1="SCM-enabled",label2="Original",timeunits='hours')
