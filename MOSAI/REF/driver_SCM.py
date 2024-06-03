#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 29 May 2024

@author: DEPHY team atelier cas 1D

Modifications:
"""

## MOSAI/REF SCM-enabled case definition

import netCDF4 as nc
import numpy as np

import argparse
parser=argparse.ArgumentParser()

from dephycf.Case import Case
from dephycf import thermo

################################################
# 0. General configuration of the present script
################################################

parser.add_argument("-c", help="cover type: MAIZE|DECIDUOUS", metavar="cover", required=True)
parser.add_argument("-n", help="name of subcase variant", metavar="name_config", default="")
args=parser.parse_args()

lplot = True     # plot the new version of the case
lcompare = True  # plot comparisons between original and new versions
lverbose = False # print information on variables and case

################################################
# 1. Get the original version of the case
################################################

# initialize the case structure for the original version
cover = args.c
scase = cover+args.n
case = Case('MOSAI/%s'%scase)

# read case information in file
case.read('MOSAI_%s_DEF_driver.nc'%scase)

# display some information about the case
if lverbose:
    case.info()

################################################
# 2. Interpolate onto a new grid, same for all the variables
#    and add new variables if needed
################################################

# Extend profiles constantly towards the surface
case.extend_init_wind(height=0)
case.extend_init_temp(height=0)
case.extend_init_rt(height=0)
case.extend_init_hur(height=0)

#case.extend_theta_advection(theta_adv=0, height=0)
#case.extend_rt_advection(rt_adv=0, height=0)
#case.extend_temperature_advection(height=0)
#case.extend_qv_advection(height=0)

####


# Grids onto which interpolate the input data

# New vertical grid, 10-m resolution from surface to 3000 m and 100-m resolution above, up to htop
htop = 50000
levout = np.array(list(range(0,3000,10)) + list(range(3100,int(htop)+1,100)),dtype=np.float64)

#  New temporal grid, from 05:00 UTC to 18:00 UTC, 20 June 2011, 30-min timestep
timeout = np.array(range(0,86400+46800+1,1800),dtype=np.float64) 

# conversion
newcase = case.convert2SCM(time=timeout,lev=levout,levtype='altitude')

# update some attributes
newcase.set_title("Forcing and initial conditions for MOSAI case - SCM-enabled version")
newcase.set_script("DEPHY-SCM/MOSAI/%s/driver_SCM.py"%scase)

# display some information about the new version of the case
if lverbose:
    newcase.info()

################################################
# 3. Save new version of the case in netcdf file
################################################

# save the new version of the case in netcdf file 
newcase.write('MOSAI_%s_SCM_driver.nc'%scase)

################################################
# 4. Plots if asked
################################################

if lplot:
    newcase.plot(rep_images='./images/driver_SCM/',timeunits='hours')

if lcompare:
    newcase.plot_compare(case,rep_images='./images/compare/',label1="SCM-enabled",label2="Original",timeunits='hours')
