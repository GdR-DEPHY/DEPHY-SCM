#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 14 Avril 2025

@author: Najda Villefranque

Modifications
"""


import netCDF4 as nc
import numpy as np

import argparse
parser=argparse.ArgumentParser()

from dephycf.Case import Case

################################################
# 0. General configuration of the present script
################################################

case_duration=60

parser.add_argument("-n", help="name of subcase variant", metavar="name_config", default="")
parser.add_argument("-p", help="plot all variables",      action="store_true")
parser.add_argument("-c", help="plot comparaisons",       action="store_true")
parser.add_argument("-v", help="verbose",                 action="store_true")
args=parser.parse_args()

scase    = args.n
lplot    = args.p
lcompare = args.c
lverbose = args.v

################################################
# 1. Get the original version of the case
################################################

# initialize the case structure for the original version
case = Case('BOTANY/%s'%scase)

# read case information in file
case.read('BOTANY_%s_DEF_driver.nc'%scase)

# display some information about the case
if lverbose:
    case.info()

################################################
# 2. Interpolate onto a new grid, same for all the variables
#    and add new variables if needed
################################################

htop = 30000
hmid = 12000
ext_height=[hmid, htop]

def get_thl_(z):
  thl=case.variables['thetal'].data[0,-2:]
  hei=case.variables['thetal'].height.data[0,-2:]
  gam = (thl[1]-thl[0])/(hei[1]-hei[0])
  return thl[1] + gam*(z-hei[1])

thl = [get_thl_(z) for z in ext_height]
case.extend_init_thetal(thetal=thl, height=ext_height)
case.extend_init_wind(u=[0,0], v=[0,0], height=ext_height)
case.extend_init_qt(qt=[0,0], height=[hmid, htop])
case.extend_vertical_velocity([0,0], height=[hmid, htop])

####


# Grids onto which interpolate the input data

# New vertical grid, 10-m resolution from surface to 4000 m and 100-m resolution above, up to htop
hmid = 4000
levout = np.array(list(range(0,int(hmid),10)) \
        + list(range(int(hmid)+100,int(htop)+1,100)),dtype=np.float64)

#  New temporal grid, case duration, 1h time step
timeout = np.arange(0.,(case_duration+1)*3600.,3600.)

# conversion
newcase = case.convert2SCM(time=timeout,lev=levout,levtype='altitude')

# update some attributes
newcase.set_title("Forcing and initial conditions for BOTANY case - SCM-enabled version")
newcase.set_script("DEPHY-SCM/BOTANY/%s/driver_SCM.py"%scase)

# display some information about the new version of the case
if lverbose:
    newcase.info()

################################################
# 3. Save new version of the case in netcdf file
################################################

# save the new version of the case in netcdf file 
newcase.write('BOTANY_%s_SCM_driver.nc'%scase)

################################################
# 4. Plots if asked
################################################

if lplot:
    newcase.plot(rep_images='./images/driver_SCM/',timeunits='hours')

if lcompare:
    newcase.plot_compare(case,rep_images='./images/compare/',label1="SCM-enabled",label2="Original",timeunits='hours')
