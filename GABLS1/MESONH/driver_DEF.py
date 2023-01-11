#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 01 May 2020

@author: Fleur Couvreux

Modifications:
    2020/06/04, R. Roehrig: Add z0 value + some cleaning/formatting
    2020/11/12, E. Vignon:  Add no evaporation option
    2021/01/03, R. Roehrig: update for improved case definition interface.

GABLS1 MESO-NH case definition (compared to REF, the TKE profile is set to 0 m2 s-2)
From Kosovic and Curry 2000; Cuxart et al 2006; Beare et al 2006
From Kosovic and Curry 2000: The initial conditions,surface cooling rate, and the inversion strength for these simulations were based on the measurements made during BASE on 1 October 1994 from flight number 7
In  the  baseline  simulation,  the  latitude  was 73.8N, the geostrophic wind was set to 8 m s-1, thesurface cooling rate was 0.25 K h-1, the overlying inversion strength was 0.01 K m-1, and the surface roughness was 0.1 m. The baseline case roughness length is higher than the typical roughness length over sea ice in the Arctic ocean. However, due to the limitations of LES resolution, using a significantly lower roughness length would  result  in an underresolved surface layer.
"""

import os

import netCDF4 as nc
import numpy as np

from dephycf.Case import Case

################################################
# 0. General configuration of the present script
################################################

lplot = True # plot all the variables
lverbose = False # print information about variables and case

################################################
# 1. General information about the case
################################################
# This is an idealized case so date is arbritrary
# lat/lon are fixed to a lat representative of Arctic conditions
# 8h duration with a constant surface cooling

case = Case('GABLS1/MESONH',
        lat=73,
        lon=123.33,
        startDate="20000101100000",
        endDate="20000101190000",
        surfaceType='land',
        zorog=0.)

case.set_title("Forcing and initial conditions for GABLS1 case - Original definition")
case.set_reference("Beare et al. (2006, BLM), Cuxart et al (2006, BLM), Kosovic and Curry (2000) reproduced in Meso-NH by Q Rodier")
case.set_author("F. Couvreux")
case.set_script("driver_DEF.py")


# time units are expected to be seconds since startDate
t0 = 0     # 10:00 UTC, 1 January 2000
t1 = 32400 # 17:00 UTC, 1 January 2000


################################################
# 2. Initial state
################################################

# Surface pressure
ps = 101320.
case.add_init_ps(ps)

#         z (m) theta (K) rt (g kg-1) u (m s-1) v (m s-1)
init = [   0.0,   265.0,   0.0,       0.0,     0.0,\
           2.0,   265.0,   0.0,       8.0,     0.0,\
         100.0,   265.0,   0.0,       8.0,     0.0,\
         400.0,   268.0,   0.0,       8.0,     0.0,\
         700.0,   271.0,   0.0,       8.0,     0.0]

init = np.array(init,dtype=np.float64)

z = init[0::5]

case.add_init_theta(init[1::5], lev=z, levtype='altitude')
case.add_init_rt(init[2::5]/1000., lev=z, levtype='altitude') # converted in kg kg-1
case.add_init_wind(u=init[3::5],v=init[4::5], lev=z, levtype='altitude')

case.add_init_height(z,lev=z,levtype='altitude')

# Turbulent Kinetic Energy
ztke = range(0,400+1,10)
nztke = len(ztke)
tke = np.zeros(nztke,dtype=np.float64)
# tke initialized to 0 in the Meso-NH LES instead of 0.4(1-z/250)**3 m 2 s-2 ,for z<250m as proposed by Beare et al (2006)

#for iz in range(0,nztke):
#    if ztke[iz] < 250:
#      tke[iz] = 0.4*(1.-ztke[iz]/250.)**3
#    else:
#      tke[iz] = 0.

case.add_init_tke(tke, lev=ztke, levtype='altitude')

################################################
# 3. Forcing
################################################

# Constant Geostrophic wind across the simulation
nz = len(z)
ug = np.zeros(nz,dtype=np.float64) + 8.
vg = np.zeros(nz,dtype=np.float64) + 0.

case.add_geostrophic_wind(ug=ug,vg=vg,lev=z,levtype='altitude')

# Surface Forcing
# constant cooling rate 0.25K/hr from 265 K

ts=[265., 264.75, 264.5, 264.25, 264., 263.75, 263.5, 263.25, 263.0, 262.75]
timets=[0., 3600., 7200., 10800., 14400., 18000., 21600., 25200., 28800., 32400.]

case.add_forcing_ts(ts,time=timets,z0=0.1)

# The following flux-gradient relations are recommended
# du/dz= ∂v=dz u*/(Kz)*(1.+Bm(z/L))
# dtheta/dz= θ*/(Kz)*(1+Bh(z/L))
# K=0.4
# Bm=4.8
# Bh=7.8

# No surface evaporation (in fact no moisture at all)
case.deactivate_surface_evaporation()

# Radiation scheme is switched off 
case.deactivate_radiation()

################################################
# 4. Writing file
################################################

case.write('GABLS1_MESONH_DEF_driver.nc')

if lverbose:
    case.info()

################################################
# 5. Ploting, if asked
################################################

if lplot:
    case.plot(rep_images='./images/driver_DEF/',timeunits='hours')
