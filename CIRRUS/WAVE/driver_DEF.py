#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 09 June 2026

@author: DEPHY team

Modification
"""

## Cirrus case (idealized) - Borella 2025

import netCDF4 as nc
import numpy as np

from datetime import datetime, timedelta
from dephycf.Case import Case
from dephycf import constants


###################################################
def MurphKoop(T):
# saturation vapor pressure with respect to liquid/ice according to Murphy and Koop
# esl valid for 123<T<332K, esi for 110<T<273.16
    esl = np.exp(54.842763 - 6763.22/T - 4.210*np.log(T)+ 0.000367*T+np.tanh(0.0415*(T-218.8))*(53.878-1331.22/T-9.44523*np.log(T)+0.014025*T))
    esi = np.exp(9.550426 - 5723.265/T + 3.53068*np.log(T)-0.00728332*T)
    return(esl,esi)

##################################################
def rhi2qv(rhi,temp,pres):
   Md=28.0*0.80+0.20*32.0
   Mw=18.0 
   esl,esi=MurphKoop(temp)
   ei=esi*rhi
   x=ei/pres
   rv=x*Mw/Md
   qv=rv/(1.+rv)
   return qv


###################################################
def read_file(file):
#file='07110.2022122700.HR_complet.csv'
#file='07145.2022122700.HR_complet.csv'

   data=np.loadtxt(file,usecols=(0,1,2,3,4,5,6,7,8),skiprows=3,delimiter=',',max_rows=3160)
   zz=data[:,1]
   pres=data[:,0]
   temp=data[:,2]
   tdew=data[:,3]
   dd=data[:,4]
   ff=data[:,5]

   esl,esi=MurphKoop(tdew)
   x=esl/pres
   Md=28.0*0.80+0.20*32.0
   Mw=18.0
   rv=x*Mw/Md
   qv=rv/(1.+rv)
   Rd=287.0
   cp=1004.0
   p0=101600.0
   theta=temp*(pres/p0)**(-Rd/cp)

   return zz,pres,temp,theta,qv,ff,dd

################################################
# 0. General configuration of the present script
################################################

lplot    = True  # plot all the variables
lverbose = False # print information about variables and case

################################################
# 1. General information about the case
################################################

duration=12
tmin = datetime(2022, 12, 27, 00, 00)
tmax = tmin + timedelta(hours=duration)

case = Case('CIRRUS/WAVE',
        lat=48.4, # BREST 
        lon=-4.5,
        startDate=tmin,
        endDate=tmax,
        surfaceType='land',
        zorog=0)

case.set_title("Forcing and initial conditions for CIRRUS")
case.set_reference("Borella et al. 2025")
case.set_author("DEPHY Team")
case.set_script("DEPHY-SCM/CIRRUS/WAVE/driver_DEF.py")

################################################
# 2. Initial state
################################################


file='07145.2022122700.HR_complet.csv'
zz,pres,temp,theta,qv,ff,dd=read_file(file)

# Surface pressure
ps = 100650.
case.add_init_ps(ps)

htop = 12000

# Zonal and meridional wind
zu = zv = [0, 600, 2700, 7500, htop]
u  = [0, 0, 0, 0, 0]
v  = [0]*len(u)

case.add_init_wind(u=np.array(u),ulev=np.array(zu),v=np.array(v),vlev=np.array(zv),levtype='altitude')

#  temperature
ztemp = np.array(zz)
ttemp  = np.array(temp)
case.add_init_temp(ttemp,lev=ztemp,levtype='altitude')

# qv
qvprof=np.array(qv)
rhi_issr=120.0
rhi_below=50.0
rhi_above=100.0
h1=9500.
h2=10000.

for k in range(len(qvprof)):
    if (zz[k]>h1 and zz[k]<h2):
        tmp=rhi2qv(rhi_issr/100.,temp[k],pres[k])
        qvprof[k]=tmp
    elif (zz[k]<h1):
        tmp=rhi2qv(rhi_below/100.,temp[k],pres[k])
        qvprof[k]=tmp
    else:
        tmp=rhi2qv(rhi_above/100.,temp[k],pres[k])
        qvprof[k]=tmp


zqv=np.array(zz)
case.add_init_qv(np.array(qvprof),lev=zqv,levtype='altitude') 

################################################
# 3. Forcing
################################################

# vertical velocities

wasc = 0.5  # m/s
h1 = 8000
h2 = 11000

zw = np.array([0., h1-100,h1, h2,h2+100, htop])

t1 = 7200
t2 = 10800
tw = np.array([0, t1-10,t1, t2,t2+10, 43200])

ww = np.zeros((len(tw), len(zw)))

mask_z = (zw >= h1) & (zw <= h2)
mask_t = (tw >= t1) & (tw <= t2)

ww[np.ix_(mask_t, mask_z)] = wasc

case.add_vertical_velocity(w=ww,time=tw,lev=zw,levtype='altitude')

# add zero tendencies 
case.add_theta_advection([0]*len(zw), lev=zw, levtype="altitude", include_rad=False)

# Surface Forcing
case.add_surface_fluxes(sens=0,lat=0,forc_wind='z0',z0=0.01)

tskin = 280
case.add_surface_skin_temp(tskin)

################################################
# 4. Writing file
################################################

case.write('CIRRUS_WAVE_DEF_driver.nc')

if lverbose:
    case.info()

################################################
# 5. Ploting, if asked
################################################

if lplot:
    case.plot(rep_images='./images/driver_DEF/',timeunits='hours')
