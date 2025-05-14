#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 29 Mai 2024

@author: DEPHY team "atelier cas 1D"

Modifications
"""

import os, sys

import netCDF4 as nc
import numpy as np

import argparse
parser=argparse.ArgumentParser()

from dephycf.Case import Case
from dephycf.other_thermo import rm_from_hu
from dephycf.utils import kernel_mean_gauss

################################################
# 0. General configuration of the present script
################################################

parser.add_argument("-c", help="cover type: MAIZE|DECIDUOUS", metavar="cover", required=True)
parser.add_argument("-n", help="name of subcase variant", metavar="name_config", default="")
parser.add_argument("-s", help="sensible heat flux threshold", metavar="thresh_sensib", default=-9999, type=float)
parser.add_argument("-a", help="flag to activate advection tendencies", action="store_true")
args=parser.parse_args()

lplot = True # plot all the variables
lverbose = False # print information about variables and case

################################################
# 1. General information about the case
################################################

start_date = "20230819050000"
end_date  =  "20230820180000"
Zorog = 594.

cover = args.c
scase = cover+args.n
sensib_thresh = args.s
advect = args.a

from datetime import datetime
timeref = datetime.strptime(start_date, "%Y%m%d%H%M%S")

case = Case('MOSAI/%s'%scase,
        lat=43.1,
        lon=0.36,
        startDate=start_date,
        endDate=end_date,
        surfaceType='land',
        zorog=Zorog)

case.set_title("Forcing and initial conditions for MOSAI case - %s case"%scase)
case.set_reference("MOSAI campaign")
case.set_author("DEPHY team")
case.set_script("DEPHY-SCM/MOSAI/%s/driver_DEF.py"%scase)

################################################
# 2. Initial state
################################################

# Surface pressure
ps = 95000. # Approximate value from ERA5
case.add_init_ps(ps)

# Thermodynamical initial profiles
# Altitude, u, v, temperature, pression, relative humidity
ds = nc.Dataset("input.nc", "r")

# utility functions to read netCDF variables, smooth the soundings and keep
# only the ascending part of the profiles
def getvar(var):
  if (var == "height"):
    return ds.variables[var][:] - Zorog
  else:
    return ds.variables[var][:]
def smoothvar(var, npts=40):
  return kernel_mean_gauss(getvar(var)[imin:imax],z,npts)

z = getvar("height")    # m
imin = 0
zdiff = 0
while zdiff < 2:
  imin+=1
  zdiff=z[imin]-z[imin-1]
imax = np.argmax(z)     # to select only the upward sounding
z = z[imin:imax]

# get all variables and smooth them
u = smoothvar("ua")     # m/s
v = smoothvar("va")     # m/s
temp = smoothvar("ta")  # °C
hur = smoothvar("hur")  # %
pre = smoothvar("pa")   # hPa
z = smoothvar("height") # m

temp = [t + 273.15 for t in temp] # => to K
# in rm_from_hu, p is in hPa, T in K and h in %
# returns g/kg
rt = [0.001*rm_from_hu(p, t, h) for (p,t,h) in zip(pre,temp,hur)]
pre = [p*100 for p in pre] # in common format, P is in Pa

# Wind initial profiles
case.add_init_wind(u=u, v=v, lev=z, levtype='altitude')

# Pressure
case.add_init_pressure(pre, lev=z, levtype='altitude')

# Temperature
case.add_init_temp(temp, lev=z, levtype='altitude')

# Relative humidity
case.add_init_hur(hur, lev=z, levtype='altitude')

# Total water mixing ratio
case.add_init_rt(rt, lev=z, levtype='altitude')

################################################
# 3. Forcing
################################################

if advect :
  # advection tendencies from ARPEGE oper files
  # I reverse altitude dimension to get height-increasing variables

  fin = nc.Dataset("tendencies.nc", "r")
  timeForc= fin['time'][:]-fin['time'][0]
  levForc = fin['height_f'][0,::-1]

  # Temperature advection
  tadv = fin['tadv'][:,::-1]
  case.add_temp_advection(tadv,time=timeForc,timeid='time',lev=levForc,levtype='altitude')

  # Specific humidity advection
  qadv = fin['qadv'][:,::-1]
  case.add_qv_advection(qadv,time=timeForc,timeid='time',lev=levForc,levtype='altitude')
else:
  # No advection forcing => set to 0
  z_tend_adv     = [0,15000]
  tadv = [0,0]
  qadv = [0,0]
  case.add_temp_advection(tadv, lev=z_tend_adv, levtype="altitude")
  case.add_qv_advection(   qadv, lev=z_tend_adv, levtype="altitude")

# No wind forcing

# No radiation
case.deactivate_radiation()

# Surface Forcings : day time latent sensible ustar tskin
flux_file = "flux.txt"
dateSfc = np.genfromtxt(flux_file,dtype=str,skip_header=1,usecols=0)
hourSfc = np.genfromtxt(flux_file,dtype=str,skip_header=1,usecols=1)
lhf     = np.genfromtxt(flux_file,dtype=float,skip_header=1,usecols=2)
shf     = np.genfromtxt(flux_file,dtype=float,skip_header=1,usecols=3)
us      = np.genfromtxt(flux_file,dtype=float,skip_header=1,usecols=4)
ts      = np.genfromtxt(flux_file,dtype=float,skip_header=1,usecols=5)

time_from_start_sec = []
latent_flux = []
sensib_flux = []
ustar = []
surface_temp = []
for it in range(0, dateSfc.shape[0]):
  datetime_str = dateSfc[it]+" "+hourSfc[it]
  datetime_object = datetime.strptime(datetime_str, '"%Y-%m-%d %H:%M:%S"')
  delta_t_sec = (datetime_object-timeref).total_seconds()
  if delta_t_sec >= 0:
    time_from_start_sec += [delta_t_sec]
    latent_flux += [lhf[it]]
    sensib_flux += [max(sensib_thresh,shf[it])]
    ustar += [us[it]]
    surface_temp += [ts[it]]

case.add_surface_fluxes(sens=sensib_flux,lat=latent_flux,time=time_from_start_sec,
        forc_wind='ustar',ustar=ustar,time_ustar=time_from_start_sec)

case.add_surface_skin_temp(surface_temp, time=time_from_start_sec)

################################################
# 4. Writing file
################################################

case.write('MOSAI_%s_DEF_driver.nc'%scase)

if lverbose:
    case.info()

################################################
# 5. Ploting, if asked
################################################

if lplot:
    case.plot(rep_images='./images/driver_DEF/',timeunits='hours')
