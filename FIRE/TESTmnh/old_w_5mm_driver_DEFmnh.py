#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 21 Septembre 2021

@author: Romain Roehrig

Modifs :
  01/03/2025 - Najda Villefranque - 3 days instead of 37h

"""

## EUROCS FIRE straotumulus case original case definition

import os

import netCDF4 as nc
import numpy as np
import math

from dephycf.Case import Case
from dephycf import thermo

################################################
# 0. General configuration of the present script
################################################

lplot = True # plot all the variables
lverbose = False # print information about variables and case

################################################
# 1. General information about the case
################################################

case = Case('FIRE/TESTmnh',
        lat=33.3,
        lon=-119.5,
        startDate="19870714080000",
        endDate="19870724080000",
        surfaceType="ocean",
        zorog=0.)

case.set_title("Forcing and initial conditions for idealized stratocumulus case")
case.set_reference("DEPHY")
case.set_author("DEPHY")
case.set_script("DEPHY-SCM/FIRE/TEST/driver_DEF.py")

################################################
# 2. Initial state
################################################

htop = 50000.
hint = 2000.

# Surface pressure
ps = 101250.
case.add_init_ps(ps)

# Zonal and meridional wind
zwind = [   0, 1200. ]
u  =    [ 5 for z in zwind]
v  =    [ 0 for z in zwind]

case.add_init_wind(u=u,ulev=zwind,v=v,vlev=zwind,levtype='altitude')
case.extend_init_wind(height=htop)

# Liquid-water potential temperature
zthetal = [  0., 595., 605., 650., 800., 1000., 1200., hint]
thetal  = [ 287.5, 287.5, 299.57, 299.91, 301.03, 302.54, 304.04, 299.5+0.0075*(hint-595)]

case.add_init_thetal(thetal,lev=zthetal,levtype='altitude')

file_std = 'info_profil_std.txt'
zstd = np.genfromtxt(file_std,dtype=float,skip_header=1,usecols=0)*1000
pstd = np.genfromtxt(file_std,dtype=float,skip_header=1,usecols=1)*100
tstd = np.genfromtxt(file_std,dtype=float,skip_header=1,usecols=2)
rhostd = np.genfromtxt(file_std,dtype=float,skip_header=1,usecols=3)
rhovstd = np.genfromtxt(file_std,dtype=float,skip_header=1,usecols=4)

rvstd = rhovstd/rhostd
#profile used in the radiation for MNH is quite moist at 3000m => very important to use a top at least of 3000m
#for i in range(len(rvstd)):
#  print('i=',zstd[i],'rvstd',rvstd[i])
# so we devide the rvstd/3
rvstd=rvstd/3.

thetastd = thermo.t2theta(p=pstd, temp=tstd)
qtstd = thermo.rt2qt(rvstd)

indz = zstd > hint
case.extend_init_thetal(thetal=thetastd[indz], height=zstd[indz])

# Total water content
zqt = [  0., 595., 605., 650., 800., 1000., 1200., hint]
qt  = [  9.6e-3, 9.6e-3, 6.57e-3, 6.43e-3, 5.99e-3, 5.39e-3, 4.79e-3, 6.6e-3-0.003e-3*(hint-595)]

case.add_init_qt(qt,lev=zqt,levtype='altitude')
indz = zstd > hint
case.extend_init_qt(qt=qtstd[indz], height=zstd[indz])

################################################
# 3. Forcing
################################################

# Large-scale velocity - constant
zw = [0., 100., 300., 500., 595., 605., 650., 800., 900., 1050., 1100., 1200.] + list(range(1300,6000,100))
#w  = [0., -0.001, -0.003, -0.005, -0.00595, -0.00605, -0.0065, -0.008, -0.009, -0.0020, 0.,0.]

w0 = -0.005
z0 = 1000

w = []
for z in zw:
   w.append(w0*math.tanh(z/z0))

case.add_vertical_velocity(w=w,lev=zw,levtype='altitude')

# Constant Geostrophic wind across the simulation
ugeo = case.variables['ua'].data[0,:]
vgeo = case.variables['va'].data[0,:]
levForc_wind = case.variables['ua'].height.data[0,:]

case.add_geostrophic_wind(ug=ugeo,uglev=levForc_wind,vg=vgeo,vglev=levForc_wind,levtype='altitude')

# Nudging of thetal and qt
thetal_nudg = case.variables['thetal'].data[0,:]
levForc_thetal = case.variables['thetal'].height.data[0,:]
case.add_thetal_nudging(thetal_nudg,timescale=3600.*3.,z_nudging=1400.,lev=levForc_thetal,levtype='altitude')

qt_nudg = case.variables['qt'].data[0,:]
levForc_qt = case.variables['qt'].height.data[0,:]
case.add_qt_nudging(qt_nudg,timescale=3600.*3.,z_nudging=1400.,lev=levForc_qt,levtype='altitude')

# Surface Forcing. Constant sea surface temperature
ts = 289.
case.add_forcing_ts(ts)

################################################
# 4. Writing file
################################################

case.write('FIRE_TESTmnh_DEF_driver.nc')

if lverbose:
    case.info()

################################################
# 5. Ploting, if asked
################################################

if lplot:
    case.plot(rep_images='./images/driver_DEF/',timeunits='hours')
