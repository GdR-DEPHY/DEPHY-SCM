#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 9 November 2022

@author: Catherine Rio

Modification
  2022/11/08, C. Rio: EUROCS case
  2024/02/05, F. Couvreux: adaptation au dernier format
  2026/06/02, N. Villefranque: clean for publication
"""

## EUROCS original case definition
## From  https://rmets.onlinelibrary.wiley.com/doi/abs/10.1256/qj.03.145

import numpy as np

from datetime import datetime, timedelta
from dephycf.Case import Case

################################################
# 0. General configuration of the present script
################################################

lplot = True # plot all the variables
lverbose = False # print information about variables and case

################################################
# 1. General information about the case
################################################

duration=36 # forcings are defined up to 96 h
tmin = datetime(1997, 6, 27, 11, 30)
tmax = tmin + timedelta(hours=duration)

case = Case('EUROCS/REF',
        lat=36.61,
        lon=-97.49,
        startDate=tmin,
        endDate=tmax,
        surfaceType='land',
        zorog=360.)

case.set_title("Forcing and initial conditions for EUROCS case - Original definition")
case.set_reference("Guichard et al. (2004, QJRMS)")
case.set_author("C. Rio")
case.set_script("DEPHY-SCM/EUROCS/REF/driver_DEF.py")

################################################
# 2. Initial state
################################################

# Surface pressure
ps = 97285.89
case.add_init_ps(ps)

# Pressure, potential temperature, specific humidity profiles from aux file
P,th,rv = np.genfromtxt('../aux/init_thrv.txt', dtype=float, skip_header=0,
                        usecols=[0,1,2]).transpose()
case.add_init_theta(th, lev=P, levtype='pressure')
case.add_init_rv(rv, lev=P, levtype='pressure')

# Pressure, horizontal winds profiles from aux file
Pu,u,v = np.genfromtxt('../aux/init_uv.txt', dtype=None, skip_header=0,
                       usecols=[0,1,2]).transpose()
case.add_init_wind(u=u,v=v, lev=Pu, levtype='pressure')

################################################
# 3. Forcing
################################################

# time varying forcings, get forcing times 
ys,ms,ds,ss = np.genfromtxt('../aux/Forc_Tq_sfc.txt', dtype=float, skip_header=0,
                            usecols=[0,1,2,3]).transpose()

TimeForc=[]
for it, (y,m,d,t) in enumerate(zip(ys,ms,ds,ss)):
  dateforc = datetime(int(y),int(m),int(d))+timedelta(hours=t/3600.)
  difftime = dateforc-tmin
  TimeForc += [difftime.total_seconds()]

TimeForc = np.array(TimeForc)
ntimes=len(TimeForc)

# read forcing variables from file, where there is one variable per column,
# each column = ntimes x nlevels values => reshape to (ntimes, nlevels)
Pfrck                = np.genfromtxt('../aux/Forc_Tq.txt', dtype=float, skip_header=0, usecols=[0]).transpose().reshape((ntimes,-1))
ufrck, vfrck         = np.genfromtxt('../aux/Forc_Tq.txt', dtype=float, skip_header=0, usecols=[1,2]).transpose().reshape((2,ntimes,-1))
thfrck, rvfrck       = np.genfromtxt('../aux/Forc_Tq.txt', dtype=float, skip_header=0, usecols=[3,4]).transpose().reshape((2,ntimes,-1))
dthdtfrck, drvdtfrck = np.genfromtxt('../aux/Forc_Tq.txt', dtype=float, skip_header=0, usecols=[6,7]).transpose().reshape((2,ntimes,-1))

Pfrck1D = Pfrck[0,:]

case.add_wind_nudging(unudg=ufrck, vnudg=vfrck,
                      timescale=3600.*2., p_nudging=110000., 
                      time=TimeForc, timeid='time', 
                      lev=Pfrck1D, levtype='pressure', levid='lev')

case.add_theta_advection(dthdtfrck, time=TimeForc, lev=Pfrck1D,
                         levtype='pressure', include_rad=False)

case.add_rv_advection(drvdtfrck, time=TimeForc, lev=Pfrck1D,
                      levtype='pressure')

# surface forcing variables
timeSfc, sensSfc, lateSfc, tsSfc = np.genfromtxt('../aux/surface_flux_forcings.txt', dtype=float, skip_header=0).transpose()

case.add_surface_fluxes(sens=sensSfc, lat=lateSfc, time=timeSfc, forc_wind='z0', z0=0.15)
case.add_rad_ts(tsSfc, time=timeSfc)

# add emissivity_lw = 0.994 and albedo_sw = 0.17
alb=0.17
emis=0.994
case.add_forcing_variable('alb',alb)
case.add_forcing_variable('emis',emis)

################################################
# 4. Writing file
################################################

case.write('EUROCS_REF_DEF_driver.nc')

if lverbose:
    case.info()

################################################
# 6. Ploting, if asked
################################################

if lplot:
    case.plot(rep_images='./images/driver_DEF/',timeunits='hours')
