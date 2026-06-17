#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 10 Juin 2026

@author: Fleur Couvreux

Modification
"""

## LBA/REF original case definition

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


case = Case('GOAMAZON/SINGLE',
        lat=-8.,
        lon=-63.,
        startDate="20141005120000",
        endDate="20141006000000",
        surfaceType='land',
        zorog=0.)


case.set_title("Forcing and initial conditions for the single pulse case - Original definition")
case.set_reference("Y. Tian, Y. Zhang 2025. Factors Controlling Precipitation Onset and Maintenance: Inferences from Large Eddy Simulations of Two Afternoon Deep Convective Regimes over Amazon Geophysical Research Letters. DOI: 10.1029/2024GL113920")
case.set_author("F. Couvreux")
case.set_script("DEPHY-SCM/GOAMAZON/SINGLE/driver_DEF.py")

################################################
# 2. Initial state
################################################

# Surface pressure
ps = 100580.
case.add_init_ps(ps)

# Initial profile in init.txt
# Altitude above the ground
z = np.genfromtxt('sing_sounding_278_5.txt',skip_header=2,dtype=np.float32,usecols=0)

# Zonal and meridional wind
u = np.genfromtxt('sing_sounding_278_5.txt',skip_header=2,dtype=np.float32,usecols=4)
v = np.genfromtxt('sing_sounding_278_5.txt',skip_header=2,dtype=np.float32,usecols=5)

case.add_init_wind(u=u, v=v, lev=z, levtype='altitude')

# Potential temperature
theta = np.genfromtxt('sing_sounding_278_5.txt',skip_header=2,dtype=np.float32,usecols=2)

case.add_init_theta(theta, lev=z, levtype='altitude')

# Water vapor mixing ratio
q = np.genfromtxt('sing_sounding_278_5.txt',skip_header=2,dtype=np.float32,usecols=3)
q=q*0.001# to transfer the specific humidity in kg/kg

case.add_init_qv(q, lev=z, levtype='altitude')

################################################
# 3. Forcing
################################################
forc = nc.Dataset('sing_forcing.nc','r')
# Forcing time axis
timeForc = forc['time'][24:]
timeForc=(timeForc-278.5)*86400. # conversion in seconds since 05/10/2004 at 12TU

print('timeForc',timeForc[0])
nt = len(timeForc)

#height_forcb = forc['height'][:,:]#2D field is this a problem?
height_forc = forc['height'][0,:]#2D field is this a problem?
#print('h_forc.shape',height_forc.shape,height_forcb[0,1],min(height_forcb[:,1]),max(height_forcb[:,1]),height_forcb[0,39])

tend_T=forc['tls'][24:,:]
tend_q=forc['qls'][24:,:]
U_forc=forc['U'][24:,:]
V_forc=forc['V'][24:,:]
W_forc=forc['W'][24:,:]

# Potential temperature advection, which includes radiation
case.add_theta_advection(tend_T,include_rad=False,lev=height_forc,levtype='altitude',time=timeForc)

# Water vapour mixing ratio advection
case.add_rv_advection(tend_q,lev=height_forc,levtype='altitude',time=timeForc)

#Vertical velocity
case.add_vertical_velocity(w=W_forc,lev=height_forc,levtype='altitude',time=timeForc)

# Wind nudging
case.add_wind_nudging(unudg=U_forc,vnudg=V_forc,timescale=3600.,time=timeForc,lev=height_forc,levtype='altitude')

# Surface fluxes
timeflux = np.genfromtxt('sing_sfc_flux_278_5.txt',skip_header=1,dtype=np.float32,usecols=0)
timeflux=(timeflux-278.5)*86400. # conversion in s and startitng from the beginning of the simulation =278.5
fsens = np.genfromtxt('sing_sfc_flux_278_5.txt',skip_header=1,dtype=np.float32,usecols=2)
flat = np.genfromtxt('sing_sfc_flux_278_5.txt',skip_header=1,dtype=np.float32,usecols=3)
ftau = np.genfromtxt('sing_sfc_flux_278_5.txt',skip_header=1,dtype=np.float32,usecols=4)
#ATTENTION BESOIN DE REGLER PRESCRIPTION De Quantite de mouvement
case.add_surface_fluxes(fsens,flat,time=timeflux,forc_wind='z0',z0=0.035)

################################################
# 4. Writing file
################################################

case.write('GOAMAZON_SINGLE_DEF_driver.nc')

if lverbose:
    case.info()

################################################
# 6. Ploting, if asked
################################################

if lplot:
    case.plot(rep_images='./images/driver_DEF/',timeunits='hours')
