#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 21 September 2021

@author: Romain Roehrig
"""

## FIRE Reference case SCM-enabled definition

import numpy as np
import netCDF4 as nc

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
case = Case('FIRE/REF')

# read case information in file
case.read('FIRE_REF_DEF_driver.nc')

# display some information about the case
if lverbose:
    case.info()

################################################
# 2. Interpolate onto a new grid, same for all the variables
#    and add new variables if needed
################################################

htop = 40000.

# Extend profiles
# Wind: just make it constant above what is presently defined
case.extend_init_wind(height=htop)

# First extend thetal up to 2500m with the same thetal vertical gradient 
# (and limit discontinuity with ERA5 profile)
zthetal = [  0., 595., 605., 650., 800., 1000., 1200.,2500.]
thetal  = [ 287.5, 287.5, 299.57, 299.91, 301.03, 302.54, 304.04,313.79]

# Then use ERA5 above
with nc.Dataset('../aux/ERA5/ERA5_FIRE_19870714000000-19870715230000.nc') as f:
    temp = np.average(f['ta'][:,::-1,0,0], axis=0) # time average over the two days
    pa = f['plev'][::-1]
    thetal = thermo.t2theta(p=pa, temp=temp)
    height = np.average(f['zg'][:,::-1,0,0], axis=0)
    case.extend_init_thetal(thetal=thetal, height=height)

# First extend qt up to 2800m with the same qt vertical gradient
zqt = [  0., 595., 605., 650., 800., 1000., 1200.,2800.]
qt  = [  9.6e-3, 9.6e-3, 6.57e-3, 6.43e-3, 5.99e-3, 5.39e-3, 4.79e-3,0.e-3]

# Then put it to zero above
case.extend_init_qt(qt=[0,0], height=[3000,htop])

# Geostrophic wind: just make it constant above what is presently defined
case.extend_geostrophic_wind(height=htop)

# Large-scale forcing vertical velocity: linear to 0 between 1200 and 1300m, 0 above
case.extend_vertical_velocity(w=[0,0],height=[1205.,htop])

# Thetal and qt large-scale advection: linear to 0 between 1200 and 1300m, 0 above
case.extend_thetal_advection(thetal_adv=[0,0], height=[1205.,htop])
case.extend_qt_advection(qt_adv=[0,0], height=[1205.,htop])

# grid onto which interpolate the input data

# New vertical grid, 10-m resolution from surface to 60000 m (above the surface)
levout = np.array(list(range(0,2001,5)) + list(range(2100,int(htop)+1,100)),dtype=np.float64) 

# New temporal grid, from 00:00 UTC, 16 December 2004 to 00:00 UTC 19 December 2004, 1-hour timestep
timeout = np.arange(0.,(37+1)*3600.,3600.)

# conversion
newcase = case.convert2SCM(time=timeout,lev=levout,levtype='altitude')

# update some attributes
newcase.set_title("Forcing and initial conditions for FIRE Reference case - SCM-enabled version")
newcase.set_script("DEPHY-SCM/FIRE/REF/driver_SCM.py")

# display some information about the new version of the case
if lverbose:
    newcase.info()

################################################
# 3. Save new version of the case in netcdf file
################################################

# save the new version of the case in netcdf file 
newcase.write('FIRE_REF_SCM_driver.nc')

################################################
# 4. Plots if asked
################################################

if lplot:
    newcase.plot(rep_images='./images/driver_SCM/',timeunits='hours')

if lcompare:
    newcase.plot_compare(case,rep_images='./images/compare/',label1="SCM-enabled",label2="Original",timeunits='hours')
