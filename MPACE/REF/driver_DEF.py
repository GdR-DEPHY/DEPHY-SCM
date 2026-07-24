#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 27 November 2019

@author: Romain Roehrig

Modification
  2020/11/11, R. Roehrig: update for improved case definition interface.
  2025/01/26, N. Villefranque: merge with MESONH and clean for publication.
"""
#MPACE original case definition
#From Klein et al. 2009; QJRMS, DOI: 10.1002/qj.416
#In  the  baseline  simulation, longitude is 209.0, latitude is 71.75. The lower initial boundary condition is specified as an ocean surface with temperature 274.01 K Models were asked to simulate the 12h starting from 1700 UTC 9 October 2004. Initial profiles of ice-liquid water temperature and total water are prescribed and correspond to a cloudy (purely liquid) convective boundary layer topped by an inversion. Surfaces heat and latent fluxes are assumed constant throughout the simulation and vertical velocity and horizontal advections of heat and water vapor are prescribed.


import numpy as np

from datetime import datetime, timedelta
from dephycf.Case import Case
from dephycf import constants

################################################
# 0. General configuration of the present script
################################################

lplot    = True  # plot all the variables
lverbose = False # print information about variables and case

################################################
# 1. General information about the case
################################################
# This is an idealized case so date is arbritrary
# lat/lon are fixed to a lat representative of Arctic conditions
# 8h duration with a constant surface cooling

case = Case('MPACE/REF',
        lat=71.75,
        lon=209.00,
        startDate="20041009170000",
        endDate="20041010050000",
        surfaceType="ocean",
        zorog=0.)

case.set_title("Forcing and initial conditions for MPACE case - Original definition")
case.set_reference("Klein et al. (2009, QJRMS)")
case.set_author("E. Vignon")
case.set_script("driver_DEF.py")
case.set_modifications("In the original case, their is no struct constraint on how to nudge the wind towards the initial profile. Here we apply a nudging time scale of 1h.")
case.set_comment("For model with an explicit aerosol-cloud coupling, it is recommended to use the bimodal lognormal \n size distribution for dry aerosols given in Klein et al. 2009. Aerosol composition was assumed to be ammonium bisulphate\n with an insoluble fraction of about 30%. A concentration value of 0.16 L-1 is also recommended for INPs in the\ndeposition, condensation-freezing, and immersion-freezing modes")


################################################
# 2. Initial state
################################################
pinv=85000. # inversion pressure in Pa
ps=101000.  # surface pressure
pf = np.logspace(np.log10(ps), np.log10(30000), 101)


# Surface pressure
ps = 101000.
case.add_init_ps(ps)
# surface temperature
ts=274.01
case.add_init_ts(ts)

# Zonal and meridional wind

u=pf*0.-13.
v=pf*0.-3
case.add_init_wind(u=u,ulev=pf,v=v,vlev=pf,levtype='pressure')

# Liquid potential temperature and total humidity
thetal = np.where(pf > pinv, 269.2, 275.33+0.0791*(815-pf/100.))
qt=np.where(pf>pinv, 1.95/1000., (0.291+0.00204*(pf/100-590))/1000.)
qt=np.maximum(0.,qt)
case.add_init_thetal(thetal,lev=pf,levtype='pressure')

case.add_init_qt(qt,lev=pf,levtype='pressure') 



################################################
# 3. Forcing
################################################

# Constant surface pressure
case.add_surface_pressure_forcing(ps,timeid='time')

# Constant SST [K]
case.add_surface_temp(ts,timeid='time')

# Surface forcing in fluxes
hs=136.5
hl=107.7
case.add_surface_fluxes(sens=hs,lat=hl,timeid='time',forc_wind='z0',z0=0.01)
# vertical velocity
Div=5.8e-6 # large scale divergence in s-1
omega=np.minimum(Div*(ps-pf),Div*(ps-pinv))

case.add_vertical_velocity(omega=omega,timeid='time',lev=pf,levtype='pressure',levid='lev')

# temperature and vapor advection

tadv=np.minimum(-4., -15.*(1.-(ps-pf)/21818.))/86400.
qadv=np.minimum(-0.164,-3.*(1.-(ps-pf)/15171.))/86400./1000

case.add_temp_advection(tadv,timeid='time',lev=pf,levtype='pressure',levid='lev')
case.add_qv_advection(qadv,timeid='time',lev=pf,levtype='pressure',levid='lev')


# nudging of wind towards constant values
u=np.zeros((len(pf)))-13.
v=np.zeros((len(pf)))-3.
case.add_wind_nudging(unudg=u,vnudg=v,timescale=3600.,p_nudging=110000.,timeid='time',lev=pf,levtype='pressure',levid='lev')



################################################
# 4. Writing file
################################################

case.write('MPACE_REF_DEF_driver.nc')

if lverbose:
    case.info()

################################################
# 5. Ploting, if asked
################################################

if lplot:
    case.plot(rep_images='./images/driver_DEF/',timeunits='hours')
