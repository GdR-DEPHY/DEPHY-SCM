#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on 5 Sept 2024

@author: Frederique Cheruy 

Modifications:
    19/05/2025, R. Roehrig: Update with the new way to add wind advection + some cleaning
"""

## DICE 3 nights ....
## 
## 


import netCDF4 as nc
import numpy as np

from dephycf.Case import Case

################################################
# 0. General configuration of the present script
################################################

lplot = True  # plot all the variables
lverbose = False  # print information about variables and case

################################################
# 1. General information about the case
################################################

case_name = 'DICE'
subcase_name = 'REF'

# case = Case('{0}/{1}'.format(case_name,subcase_name),
case = Case('{0}/{1}'.format(case_name, subcase_name),
            lat=37.65,
            lon=263.26,
            startDate="19991023185959",
            endDate="19991026185959",
            surfaceType='land',
            zorog=436.)

case.set_title("Forcing and initial conditions for DICE case")
case.set_reference("Svenson et al. (2011, BLM)")
case.set_author("A. Lock, M. Best. Conversion to DEPHY standard by F. Cheruy")
case.set_script("DEPHY-SCM/{0}/{1}/driver_DEF.py".format(case_name, subcase_name))

################################################
# 2. Initial netcdf file
################################################

fin = nc.Dataset('dice_driver.nc', 'r')

# Surface pressure
ps = 97509.  # dans dice_driver
case.add_init_ps(ps)

nlev, = fin['height'].shape

# Height
height = fin['height'][:]  # reverse altitude order
case.add_init_height(height, lev=height, levtype='altitude')

# PREssure
pressure = fin['pf'][:]
case.add_init_pressure(pressure, lev=height, levtype='altitude')

# Zonal and meridional wind (same value for the first two levels)
u = fin['u'][:]
v = fin['v'][:]

case.add_init_wind(u=u, v=v, lev=height, levtype='altitude')

# Potential temperature
theta = fin['theta'][:]
case.add_init_theta(theta, lev=height, levtype='altitude')

# Vapor water content
qv = fin['qv'][:]
case.add_init_qv(qv, lev=height, levtype='altitude')

# Water vapor mixing ratio
# rv = rv*0.001 # conversion g/kg en kg/kg

# case.add_init_rv(rv, lev=z, levtype='altitude')

################################################
# 3. Forcing
################################################

timeForc = fin['time'][:]
nt, = timeForc.shape
# ici c est la dimension qui ne va pas
# On veut que la premiere dimension soit la vertical
# /!\ on rajoute le niveau de surface
ug_lev = fin['Ug'][:]
vg_lev = fin['Vg'][:]

# Convert to 2D and swap axes so that time is first
ug = np.swapaxes(np.tile(ug_lev, (nlev, 1)), 0, 1)
vg = np.swapaxes(np.tile(vg_lev, (nlev, 1)), 0, 1)
case.add_geostrophic_wind(ug=ug, vg=vg, time=timeForc, uglev=height, vglev=height, levtype='altitude')

# large-scale vertical velocity
w = fin['w'][:, :]
case.add_vertical_velocity(w=w, lev=height, levtype='altitude', time=timeForc)

omega = fin['omega'][:, :]
case.add_vertical_velocity(omega=omega, lev=height, levtype='altitude', time=timeForc)

# temperature advection and radiation
hT = fin['hadvT'][:, :]
hT = hT / 86400.  # conversion K/day en K/s
case.add_temp_advection(hT, include_rad=False, time=timeForc, lev=height, levtype='altitude')

# Moisture advection
hq = fin['hadvq'][:, :]
hq = hq / 86400.  # conversion kg/kg/day en kg/kg/s
case.add_qv_advection(hq, time=timeForc, lev=height, levtype='altitude')

# wind advection
hu = fin['hadvu'][:, :]
hu = hu / 86400.  # conversion m/s/day en m/s2
hv = fin['hadvv'][:, :]
hv = hv / 86400.  # conversion m/s/day en m/s2
case.add_wind_advection(ua_adv=hu, va_adv=hv, lev=height, time=timeForc, levtype='altitude')

# No radiation
#case.deactivate_radiation()

# Surface fluxes
sens = fin['shf'][:]
flat = fin['lhf'][:]
# voir pour ustar
ustar = fin['ustar'][:]
case.add_surface_fluxes(sens, flat, time=timeForc, forc_wind='ustar', time_ustar=timeForc, ustar=ustar)

# Add surface temperature
ts = fin['Tg'][:]
case.add_forcing_ts(ts, time=timeForc)


################################################
# 4. Writing file
################################################

case.write('DICE_REF_DEF_driver.nc')

if lverbose:
    case.info()

################################################
# 5. Ploting, if asked
################################################

if lplot:
    case.plot(rep_images='./images/driver_DEF/', timeunits='hours')
