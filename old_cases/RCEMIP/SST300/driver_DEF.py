#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 2 Decembre 2024

@author: Romain Roehrig

Modification
"""

## RCEMIP case definition - SST = 300 K

SST = 300

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

dmin = datetime(1979, 1, 1, 0, 0, 0)
dmax = dmin + timedelta(days=50)

case = Case(f'RCEMIP/SST{SST}',
        lat=0,
        lon=0,
        case_type='RCE',
        startDate=dmin,
        endDate=dmax,
        surfaceType="ocean",
        zorog=0.)

case.set_title("Forcing and initial conditions for Radiative-Convection Equilibrium MIP (RCEMIP) case - SST=300K - Original definition")
case.set_reference("Wing et al. (2018, GMD)")
case.set_author("R. Roehrig")
case.set_script(f"DEPHY-SCM/RCEMIP/SST{SST}/driver_DEF.py")

################################################
# 2. Initial state
################################################

q0 = {}
q0[295] = 12.00/1000.
q0[296] = 13.8/1000. # Perso
q0[297] = 14.65/1000. # Perso
q0[298] = 15.55/1000. # Perso
q0[299] = 16.5/1000. # Perso
q0[300] = 18.65/1000.
q0[301.15] = 18.8/1000. # Perso
q0[301] = 18.6/1000. # Perso
q0[302] = 19.7/1000. # Perso
q0[303] = 20.9/1000. # Perso
q0[304] = 22.2/1000. # Perso
q0[305] = 24.00/1000.

# Surface pressure
ps = 101480.
case.add_init_ps(ps)

def q(z,SST=300,qt=1.e-14,zq1=4000.,zq2=7500.,zt=15000.):
  tmp = q0[SST]*np.exp(-z/zq1)*np.exp(-z*z/(zq2*zq2))
  if isinstance(tmp,float):
    if z > zt:
        tmp = qt
  else:
    tmp = np.where(z > zt, qt, tmp)
  return tmp

def Tv(z,SST=300.,zt=15000.,gamma=0.0067):
  Tv0 = SST*(1+0.608*q(0,SST=SST))
  tmp = Tv0 - gamma*z
  Tvt = Tv0 - gamma*zt
  if isinstance(tmp,float):
    if z > zt:
        tmp = Tvt
  else:
    tmp = np.where(z > zt, Tvt, tmp)
  return tmp

def T(z,SST=300):
  tmp = Tv(z,SST=SST)/(1.+0.608*q(z,SST=SST))
  return tmp

def p(z,SST=300,zt=15000.,p0=ps,gamma=0.0067,g=9.79764,Rd=287.04):
  Tv0 = Tv(0,SST=SST)
  Tvt = Tv(zt,SST=SST)
  tmp1 = p0*np.exp(g/(Rd*gamma)*np.log((Tv0-gamma*z)/Tv0))
  pt = p0*(Tvt/Tv0)**(g/(Rd*gamma))
  tmp2 = pt*np.exp(-g*(z-zt)/(Rd*Tvt))
  if isinstance(tmp1,float):
    if z > zt:
        tmp = tmp2
    else:
        tmp = tmp1
  else:  
    tmp = np.where(z > zt, tmp2, tmp1)
  return tmp


# Altitude axis
zlev = np.arange(0,80001,10)
nlev, = zlev.shape

u = np.zeros(nlev,np.float32)
v = np.zeros(nlev,np.float32)
case.add_init_wind(u=u,ulev=zlev,v=v,vlev=zlev,levtype='altitude')

# Temperature
temp  = T(zlev, SST=SST)
case.add_init_temp(temp,lev=zlev,levtype='altitude')

# Specific humidity
qv = q(zlev, SST=SST)
case.add_init_qv(qv,lev=zlev,levtype='altitude')

# Pressure
pres = p(zlev, SST=SST)
case.add_init_pressure(pres,lev=zlev,levtype='altitude')

################################################
# 3. Forcing
################################################

# Surface Forcing. Constant sea surface temperature
ts = SST
case.add_forcing_ts(ts)

# Activate radiation
case.activate_radiation()
case.add_forcing_variable('alb', 0.07)
case.add_forcing_variable('i0', 551.58)
case.add_forcing_variable('sza', 42.05)

################################################
# 4. Writing file
################################################

case.write(f'RCEMIP_SST{SST}_DEF_driver.nc')

if lverbose:
    case.info()

################################################
# 5. Ploting, if asked
################################################

if lplot:
    case.plot(rep_images='./images/driver_DEF/',timeunits='hours')
