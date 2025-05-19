#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on 5 September 2024

@author: Frederique Cheruy

Modifications:
"""

## DICE/REF SCM-enabled case definition

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
case = Case('DICE/REF')

# read case information in file
case.read('DICE_REF_DEF_driver.nc')

# display some information about the case
if lverbose:
    case.info()

################################################
# 2. Interpolate onto a new grid, same for all the variables
#    and add new variables if needed
################################################

# Extend profiles constantly towards the surface
case.extend_init_wind(height=0)
case.extend_init_theta(height=0)
case.extend_init_qv(height=0)

case.extend_geostrophic_wind(height=0)
case.extend_temperature_advection(height=0)
case.extend_qv_advection(height=0)
case.extend_wind_advection(height=0)


# Grids onto which interpolate the input data

# New vertical grid, 10-m resolution from surface to 3000 m and 100-m resolution above, up to htop
htop = 30000
levout = np.array(list(range(0,3000,10)) + list(range(3100,int(htop)+1,100)),dtype=np.float64)

#  New temporal grid, 30-min timestep
timeout = np.array(range(0,261000+1,1800),dtype=np.float64) 

# conversion
newcase = case.convert2SCM(time=timeout,lev=levout,levtype='altitude')

# update some attributes
newcase.set_title("Forcing and initial conditions for DICE case - SCM-enabled version")
newcase.set_script("DEPHY-SCM/DICE/REF/driver_SCM.py")

# display some information about the new version of the case
if lverbose:
    newcase.info()

################################################
# 3. Save new version of the case in netcdf file
################################################

# save the new version of the case in netcdf file 
newcase.write('DICE_REF_SCM_driver.nc')

################################################
# 4. Plots if asked
################################################

if lplot:
    newcase.plot(rep_images='./images/driver_SCM/',timeunits='hours')

if lcompare:
    newcase.plot_compare(case,rep_images='./images/compare/',label1="SCM-enabled",label2="Original",timeunits='hours')
