#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 10 Avril 2020

@author: Fleur Couvreux

Modifications:
  2026/07/20, N. Villefranque: clean for publication.
"""

## AYOTTE/00SC SCM-enabled case definition

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

scase = "00SC"

# initialize the case structure for the original version
case = Case('AYOTTE/%s'%scase)

# read case information in file
case.read('AYOTTE_%s_DEF_driver.nc'%scase)

# display some information about the case
if lverbose:
    case.info()

################################################
# 2. Interpolate onto DEF grid, same for all the variables
#    and add new variables if needed
################################################

# conversion
newcase = case.convert2SCM()

# add a surface temperature. To be improved...
ts = [310., 310.]
newcase.add_variable('ts',ts,time=[0, 7*3600],timeid='time')

# update some attributes
newcase.set_title("Forcing and initial conditions for AYOTTE/%s case - SCM-enabled version"%scase)
newcase.set_script("driver_SCM.py")

# display some information about the new version of the case
if lverbose:
    newcase.info()

################################################
# 3. Save new version of the case in netcdf file
################################################

# save the new version of the case in netcdf file 
newcase.write('AYOTTE_%s_SCM_driver.nc'%scase)

################################################
# 4. Plots if asked
################################################

if lplot:
    newcase.plot(rep_images='./images/driver_SCM/',timeunits='hours')

if lcompare:
    newcase.plot_compare(case,rep_images='./images/compare/',label1="SCM-enabled",label2="Original",timeunits='hours')
