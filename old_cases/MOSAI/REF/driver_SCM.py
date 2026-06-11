#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 29 May 2024

@author: DEPHY team atelier cas 1D

Modifications:
27/05/2025 Alice Maison : various options added, including geostrophic wind, wind advection and forcing smoothing
"""

## MOSAI/REF SCM-enabled case definition

import netCDF4 as nc
import numpy as np

import argparse
parser=argparse.ArgumentParser()

from dephycf.Case import Case
from dephycf import thermo
from driver_DEF import str2bool

################################################
# 0. General configuration of the present script
################################################
parser.add_argument("cover", type=str, metavar="cover", choices=['MAIZE','DECIDUOUS'], help="cover type: MAIZE|DECIDUOUS")
parser.add_argument("initpf", type=str, metavar="initpf", choices=['rs_smooth','rs_idea'], help="initial profile: rs_smooth|rs_idea")
parser.add_argument("--advTq", type=str2bool, nargs='?', const=True, default=False, help="flag to activate advection tendencies of temperature and humidity")
parser.add_argument("--advuv", type=str2bool, nargs='?', const=True, default=False, help="flag to activate advection tendencies of horizontal wind")
parser.add_argument("--vertvel", type=str2bool, nargs='?', const=True, default=False, help="flag to add vertical velocity")
parser.add_argument("fadv", type=str, metavar="adv_from", choices=['AROME','ARPEGEoper','ERA5'], help="advection from: AROME|ARPEGEoper|ERA5")
parser.add_argument("pt_AROME", type=int, default=0, metavar="AROME_pt_nb", help="miniAROME point number")
parser.add_argument("sadv", type=float, default=0., metavar="smooth_adv", help="smooth advection tendencies")
parser.add_argument("zadv", type=float, default=50000., metavar="max_alt_adv", help="maximum altitude of advection tendencies")
parser.add_argument("--geo", type=str2bool, nargs='?', const=True, default=False, help="flag to activate geostrophic wind")
parser.add_argument("sgeo", type=float, default=0., metavar="smooth_geos_wind", help="smooth geostrophic wind")
parser.add_argument("zgeo", type=float, default=16000., metavar="max_alt_geos_wind", help="maximum altitude of geostrophic wind")
parser.add_argument("--rad", type=str2bool, nargs='?', const=True, default=False, help="flag to activate radiation")
parser.add_argument("--ffx", type=str2bool, nargs='?', const=True, default=False, help="flag to activate surface flux forcing")
parser.add_argument("--fts", type=str2bool, nargs='?', const=True, default=False, help="flag to activate surface temperature forcing")
parser.add_argument("-s", type=float, default=-9999, metavar="thresh_sensib", help="sensible heat flux threshold")
args=parser.parse_args()

lplot = True     # plot the new version of the case
lcompare = True  # plot comparisons between original and new versions
lverbose = False # print information on variables and case

################################################
# 1. Get the original version of the case
################################################

# initialize the case structure for the original version
cover = args.cover
initial_profile = args.initpf
advTq = args.advTq
advuv = args.advuv
vertvel = args.vertvel
adv_from = args.fadv
zmax_adv = args.zadv
geoswd = args.geo

# name of the case
add_nam = ''
if ((advTq) | (advuv) | (vertvel)):
  if ((advTq) | (advuv)): add_nam = add_nam+'_adv'
  if (advTq): add_nam = add_nam+'Tq'
  if (advuv): add_nam = add_nam+'uv'
  if (vertvel): add_nam = add_nam+'_w'
  add_nam = add_nam+'_'+adv_from[0:3]
if (geoswd): add_nam = add_nam+'_geo'
scase = cover+'_'+args.initpf+add_nam

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


# Grids onto which interpolate the input data
# New vertical grid, 10-m resolution from surface to 3000 m and 100-m resolution above, up to htop
htop = 50000

# set advection profiles to 0 towards the top of the domain
if (advTq):
  case.extend_temperature_advection(temp_adv=[0.,0.],height=[zmax_adv,htop])
  case.extend_qv_advection(qv_adv=[0.,0.],height=[zmax_adv,htop])
if (advuv):
  case.extend_wind_advection(ua_adv=[0.,0.],va_adv=[0.,0.],height=[zmax_adv,htop])
if (vertvel):
  if ((adv_from == 'ARPEGEoper') | (adv_from == 'ERA5')):
    case.extend_vertical_velocity(omega=[0.,0.],height=[zmax_adv,htop])
  elif (adv_from == 'AROME'):
    case.extend_vertical_velocity(w=[0.,0.],height=[zmax_adv,htop])

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
