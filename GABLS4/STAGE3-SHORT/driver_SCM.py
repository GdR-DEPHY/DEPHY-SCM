#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 09 June 2020

@author: Romain Roehrig
"""

## GABLS4/STAGE3 SCM-enabled case definition

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

case_name = 'GABLS4'
subcase_name = 'STAGE3-SHORT'

# initialize the case structure for the original version
case = Case('{0}/{1}'.format(case_name,subcase_name))

# read case information in file
case.read('{0}_{1}_DEF_driver.nc'.format(case_name,subcase_name))

# display some information about the case
if lverbose:
    case.info()

################################################
# 2. Interpolate onto a new grid, same for all the variables
#    and add new variables if needed
################################################

# grid onto which interpolate the input data

# New vertical grid, 1-m resolution from the surface to 200 m and then 10-m resolution up to 30000 m (above the surface)
levout = np.array(list(range(0,201,1))+list(range(210,30001,10)),dtype=np.float64) 

# New temporal grid, from 10:00 UTC, 11 December 2009 to 22:00 UTC 11 December 2009, 1-hour timestep
timeout = np.arange(0.,(12+1)*3600.,3600.)

# conversion
newcase = case.convert2SCM(time=timeout,lev=levout,levtype='altitude')

# update some attributes
newcase.set_title("Forcing and initial conditions for {0}/{1} case - SCM-enabled version".format(case_name,subcase_name))
newcase.set_script("DEPHY-SCM/{0}/{1}/driver_SCM.py".format(case_name,subcase_name))

# display some information about the new version of the case
if lverbose:
    newcase.info()

################################################
# 3. Save new version of the case in netcdf file
################################################

# save the new version of the case in netcdf file 
newcase.write('{0}_{1}_SCM_driver.nc'.format(case_name,subcase_name))

################################################
# 4. Plots if asked
################################################

if lplot:
    newcase.plot(rep_images='./images/driver_SCM/',timeunits='hours')

if lcompare:
    newcase.plot_compare(case,rep_images='./images/compare/',label1="SCM-enabled",label2="Original",timeunits='hours')
