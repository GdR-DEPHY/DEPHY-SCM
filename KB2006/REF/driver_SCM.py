#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 06 Avril 2022

@author: Fleur Couvreux

Modification
  2026/07/15, N. Villefranque: clean for publication.
"""

## KB2006 SCM-enabled case definition

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
case = Case('KB2006/REF')

# read case information in file
case.read('KB2006_REF_DEF_driver.nc')

fin = nc.Dataset('KB2006_REF_DEF_driver.nc','r')

# display some information about the case
if lverbose:
    case.info()

################################################
# 2. Interpolate onto a new grid, same for all the variables
#    and add new variables if needed
################################################

# Grid onto which interpolate the input data

htop=20000
# New vertical grid, 10-m resolution from surface to 6000 m (above the surface)
levout = np.array(list(range(0,6000,10)) + list(range(6100,int(htop)+1,100)),dtype=float)

# New temporal grid, every half hour
timeout = np.array(range(0,86400*7+3600,3600),dtype=float) 

# Conversion
newcase = case.convert2SCM(time=timeout,lev=levout,levtype='altitude')

# add a surface temperature because some models need it even when 
# the case is forced with fluxes. To be improved...
ts = timeout*0. + 310 # same shape as timeout
newcase.add_surface_temp(ts,time=timeout,timeid='time')

# Update some attributes
newcase.set_title("Forcing and initial conditions for KB2006 case - SCM-enabled version")
newcase.set_script("DEPHY-SCM/KB2006/REF/driver_SCM.py")

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
