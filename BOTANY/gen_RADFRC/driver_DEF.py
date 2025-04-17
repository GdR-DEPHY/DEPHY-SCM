#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 11 Avril 2025

@author: Najda Villefranque

Modifications
"""

import numpy as np

import argparse
parser=argparse.ArgumentParser()

from dephycf.Case import Case
from dephycf.thermo import theta2t

################################################
# 0. General configuration of the present script
################################################

parser.add_argument("-n", help="name of subcase variant", metavar="name_config", default="")
parser.add_argument("-p", help="plot all variables",      action="store_true")
parser.add_argument("-v", help="verbose",                 action="store_true")
args=parser.parse_args()

scase    = args.n
lplot    = args.p
lverbose = args.v

################################################
# 1. General information about the case
################################################

# 24h
start_date = "20200201000000"
end_date  =  "20200202000000"
Zorog = 0

from datetime import datetime
timeref = datetime.strptime(start_date, "%Y%m%d%H%M%S")

case = Case('BOTANY/%s'%scase,
        lat=13.1,
        lon=-52,
        startDate=start_date,
        endDate=end_date,
        surfaceType='ocean',
        zorog=Zorog)

case.set_title("Forcing and initial conditions for BOTANY case - %s case ; 24h only, rad tend forced"%scase)
case.set_reference("Cloud Botany, Jansson et al. 2023, https://doi.org/10.1029/2023MS003796")
case.set_author("N Villefranque")
case.set_script("DEPHY-SCM/BOTANY/%s/driver_DEF.py"%scase)

################################################
# 2. Initial state
################################################

# Surface pressure
ps = 101605.
case.add_init_ps(ps)

# Initial profiles
file="profiles_init.txt"

# Altitude, theta_l, qt, u, w_ls
z, thl, qt, u, w_ls   = np.genfromtxt(file,dtype=float,skip_header=0,usecols=[0,1,2,3,4]).transpose()

# Wind initial profiles
case.add_init_wind(u=u, v=u*0, lev=z, levtype='altitude')

# Temperature
case.add_init_thetal(thl, lev=z, levtype='altitude')

# Total water mixing ratio
case.add_init_qt(qt, lev=z, levtype='altitude')

################################################
# 3. Forcing
################################################

# advection tendencies 
file="profiles_tendencies.txt"

# Altitude, theta_l, qt
z, thl_adv, qt_adv, tau_h = np.genfromtxt(file,dtype=float,skip_header=0,usecols=[0,1,2,3]).transpose()

tau_s = tau_h*3600
nudging_coefficient = 1/tau_s

# radiative tendencies
file="profiles_radtend.txt"
radtend = np.genfromtxt(file,dtype=float,skip_header=0).transpose() # t,z
zz = radtend[0, 1:]
tt = radtend[1:, 0]
thl_rad = radtend[1:,1:]

thl_frc = thl_rad+thl_adv

# Temperature advection
case.add_thetal_advection(thl_frc, lev=zz, time=tt, levtype='altitude', include_rad=True)
# Temperature nudging
case.add_thetal_nudging(thl, lev=z, levtype='altitude',
        nudging_coefficient=nudging_coefficient, lev_coef=z)

# Total water mixing ratio advection
case.add_qt_advection(qt_adv, lev=z, levtype='altitude')
# Total water mixing ratio nudging
case.add_qt_nudging(qt, lev=z, levtype='altitude',
        nudging_coefficient=nudging_coefficient, lev_coef=z)

# Wind forcing
case.add_wind_nudging(unudg=u,vnudg=0*u, lev=z, levtype='altitude',
        nudging_coefficient=nudging_coefficient, lev_coef=z)

# Vertical velocity forcing
case.add_vertical_velocity(w=w_ls, lev=z, levtype='altitude')

# Surface Forcings
# Constant sea surface temperature, 1.25 warmer than first level
albedo = 0.065
emissi = 0.96
z0 = 0.001
ths = thl[0]+1.25
ts = theta2t(p=ps, theta=ths)
case.add_forcing_ts(ts, z0=z0)
case.add_forcing_variable('alb',albedo)
case.add_forcing_variable('emis',emissi)

################################################
# 4. Writing file
################################################

case.write('BOTANY_%s_DEF_driver.nc'%scase)

if lverbose:
    case.info()

################################################
# 5. Ploting, if asked
################################################

if lplot:
    case.plot(rep_images='./images/driver_DEF/',timeunits='hours')
