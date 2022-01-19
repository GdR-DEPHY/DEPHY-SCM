#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on 10 January 2022

@author: Maike Ahlgrimm

Modification

"""

## MAGIC Leg15 case definition, based on ship-following
## SCM forcing. Original ship-following LES forcing
## derived by Jeremy McGibbon and Chris Bretherton
## (JAMES 2017), adapted for SCM use by M Ahlgrimm

import os
import sys
sys.path = ['../../utils/',] + sys.path

import netCDF4 as nc
from netCDF4 import Dataset
import numpy as np

from datetime import datetime, timedelta

#import constants
from Case import Case
import xarray as xr

from scipy import interpolate
import matplotlib.pyplot as plt

def ___convert_to_datetime(d):
#    ss=len(d)
  return datetime.strptime(np.datetime_as_string(d,unit='s'), '%Y-%m-%dT%H:%M:%S')

################################################
# 0. General configuration of the present script
################################################

lplot = False # plot all the variables
lverbose = False # print information about variables and case
lblend = False

################################################
# 1. General information about the case
################################################

case = Case('MAGIC/LEG15A',
        startDate="20130720173000",
        endDate="20130725040000",
        surfaceType='ocean',
        zorog=0.)

case.set_title("Forcing and initial conditions for ship-following MAGIC LEG15A case")
case.set_reference("J. McGibbon, C. Bretherton (JAMES 2017)")
case.set_author("M. Ahlgrimm")
case.set_script("DEPHY-SCM/MAGIC/LEG15A/driver_DEF.py")

# radiosonde profiles
snd=xr.open_dataset('../aux/15A/15A_snd.nc')
# large-scale forcing derived from IFS model data
lsf=xr.open_dataset('../aux/15A/15A_lsf.nc')
# SST
sfc=xr.open_dataset('../aux/15A/15A_sfc.nc')

################################################
# 2. Initial state
################################################
basetime=datetime(2013,7,20,17,30)

sndtimeorig=snd.time #minutes since 2013-07-20 17:29:00
lsftimeorig=lsf.time #hours   since 2013-07-20 13:00:00
sfctimeorig=sfc.time #minutes since 2013-07-20 13:15:00

# Convert times into seconds since 2013-07-20T173000
delta=___convert_to_datetime(sndtimeorig[0])-basetime
sndtime=delta.total_seconds()
for i in np.arange(1,len(sndtimeorig),1):
  delta=___convert_to_datetime(sndtimeorig[i])-basetime
  sndtime= np.append(sndtime,delta.total_seconds())

delta=___convert_to_datetime(lsftimeorig[0])-basetime
lsftime=delta.total_seconds()
for i in np.arange(1,len(lsftimeorig),1):
  delta=___convert_to_datetime(lsftimeorig[i])-basetime
  lsftime= np.append(lsftime,delta.total_seconds())

delta=___convert_to_datetime(sfctimeorig[0])-basetime
sfctime=delta.total_seconds()
for i in np.arange(1,len(sfctimeorig),1):
  delta=___convert_to_datetime(sfctimeorig[i])-basetime
  sfctime= np.append(sfctime,delta.total_seconds())  
#
# Model-derived forcing is available from 13UTC, but
# the initial profile from sonde only at 17:30UTC.
# The case should start at 17:30

#Profiles from sonde
zsnd_temp=snd.z.values         #m
qsnd_temp=snd.q.values/1000.   #kg/kg
psnd_temp=snd.p.values*100.    #Pa
usnd_temp=snd.u.values         #m/s
vsnd_temp=snd.v.values         #m/s
tsnd_temp=snd.T.values         #K
thetasnd_temp=snd.theta.values #K
rhsnd_temp=snd.RH.values       #1

# Initial state from surface file
sst0 = sfc.sst[255].values
case.add_init_ts(sst0)

# Since sonde profile starts at 7.5m above surface,
# add values at surface
nsnd,nlevsnd=tsnd_temp.shape
nlevsnd=nlevsnd+1
zsnd=np.zeros((nlevsnd),dtype=np.float64)
qsnd=np.zeros((nsnd,nlevsnd),dtype=np.float64)
psnd=np.zeros((nsnd,nlevsnd),dtype=np.float64)
usnd=np.zeros((nsnd,nlevsnd),dtype=np.float64)
vsnd=np.zeros((nsnd,nlevsnd),dtype=np.float64)
tsnd=np.zeros((nsnd,nlevsnd),dtype=np.float64)
thetasnd=np.zeros((nsnd,nlevsnd),dtype=np.float64)
rhsnd=np.zeros((nsnd,nlevsnd),dtype=np.float64)

zsnd[1:] = zsnd_temp
zsnd[0]  = 0.        # add lowest height level 0m

for istep in range(0,len(sndtime)):
  qsnd[istep,1:]=qsnd_temp[istep,:]
  psnd[istep,1:]=psnd_temp[istep,:]
  usnd[istep,1:]=usnd_temp[istep,:]
  vsnd[istep,1:]=vsnd_temp[istep,:]
  tsnd[istep,1:]=tsnd_temp[istep,:]
  thetasnd[istep,1:]=thetasnd_temp[istep,:]
  rhsnd[istep,1:]=rhsnd_temp[istep,:]
  
  qsnd[istep,0]     = qsnd_temp[istep,0]      # assume same value as at 7.5m
  psnd[istep,0]     = psnd_temp[istep,0]+10.  # assume same value as at 7.5m plus 10Pa
  usnd[istep,0]     = 0.                      # assume no-slip surface
  vsnd[istep,0]     = 0.                      # assume no-slip surface
  tsnd[istep,0]     = tsnd_temp[istep,0]      # assume same value as at 7.5m
  thetasnd[istep,0] = thetasnd_temp[istep,0]  # assume same value as at 7.5m
  rhsnd[istep,0]    = rhsnd_temp[istep,0]     # assume same value as at 7.5m


# Large-scale state - model derived
zlsf=lsf.z.values[4:,:] #m
ulsf=lsf.u.values[4:,:] #m/s
vlsf=lsf.v.values[4:,:] #m/s
qlsf=lsf.q.values[4:,:] #kg/kg
tlsf=lsf.T.values[4:,:] #k

# In order to merge sonde data with model data, first
# interpolate hourly model data in time to sonde time

# dimensions of model data
ntim,nlevin=zlsf.shape
# time dimension of sonde data
nsnd,nlevsnd=tsnd.shape

ulsf_sondetime=np.zeros((nsnd,nlevin),dtype=np.float64)
vlsf_sondetime=np.zeros((nsnd,nlevin),dtype=np.float64)
qlsf_sondetime=np.zeros((nsnd,nlevin),dtype=np.float64)
tlsf_sondetime=np.zeros((nsnd,nlevin),dtype=np.float64)
zlsf_sondetime=np.zeros((nsnd,nlevin),dtype=np.float64)

for ilev in range(0,nlevin):
  ff = interpolate.interp1d(lsftime[4:], tlsf[:,ilev],
                            bounds_error=False, fill_value=tlsf[-1,ilev]) # Pad after end date with the last value, if necessary
  tlsf_sondetime[:,ilev] = ff(sndtime)
  vv = interpolate.interp1d(lsftime[4:], vlsf[:,ilev],
                            bounds_error=False, fill_value=vlsf[-1,ilev]) # Pad after end date with the last value, if necessary
  vlsf_sondetime[:,ilev] = vv(sndtime)
  uu = interpolate.interp1d(lsftime[4:], ulsf[:,ilev],
                            bounds_error=False, fill_value=ulsf[-1,ilev]) # Pad after end date with the last value, if necessary
  ulsf_sondetime[:,ilev] = uu(sndtime)
  qq = interpolate.interp1d(lsftime[4:], qlsf[:,ilev],
                            bounds_error=False, fill_value=qlsf[-1,ilev]) # Pad after end date with the last value, if necessary
  qlsf_sondetime[:,ilev] = qq(sndtime)
  zz = interpolate.interp1d(lsftime[4:], zlsf[:,ilev],
                            bounds_error=False, fill_value=zlsf[-1,ilev]) # Pad after end date with the last value, if necessary
  zlsf_sondetime[:,ilev] = zz(sndtime)

  
    
# Now interpolate model data to sonde vertical grid
ulsf_sonde=np.zeros((nsnd,nlevsnd),dtype=np.float64)
vlsf_sonde=np.zeros((nsnd,nlevsnd),dtype=np.float64)
qlsf_sonde=np.zeros((nsnd,nlevsnd),dtype=np.float64)
tlsf_sonde=np.zeros((nsnd,nlevsnd),dtype=np.float64)
for istep in range(0,nsnd):
  zlsf_sondetime[istep,:]=zlsf_sondetime[10,:]
  ff = interpolate.interp1d(zlsf_sondetime[istep,:], tlsf_sondetime[istep,:],
                            bounds_error=False, fill_value=tlsf_sondetime[istep,-1]) # Pad lowest level, if necessary
  tlsf_sonde[istep,:] = ff(zsnd)
  uu = interpolate.interp1d(zlsf_sondetime[istep,:], ulsf_sondetime[istep,:],
                            bounds_error=False, fill_value=ulsf_sondetime[istep,-1]) # Pad lowest level, if necessary
  ulsf_sonde[istep,:] = uu(zsnd)
  vv = interpolate.interp1d(zlsf_sondetime[istep,:], vlsf_sondetime[istep,:],
                            bounds_error=False, fill_value=vlsf_sondetime[istep,-1]) # Pad lowest level, if necessary
  vlsf_sonde[istep,:] = vv(zsnd)
  qq = interpolate.interp1d(zlsf_sondetime[istep,:], qlsf_sondetime[istep,:],
                            bounds_error=False, fill_value=qlsf_sondetime[istep,-1]) # Pad lowest level, if necessary
  qlsf_sonde[istep,:] = qq(zsnd)

# Define "background" fields and blend sonde data with model data
ubg=np.zeros((nsnd,nlevsnd),dtype=np.float64)
vbg=np.zeros((nsnd,nlevsnd),dtype=np.float64)
tbg=np.zeros((nsnd,nlevsnd),dtype=np.float64)
qbg=np.zeros((nsnd,nlevsnd),dtype=np.float64)

for istep in range(0,nsnd):
  #Blend model data in with sonde data in lowest 5 levels
  ubg[istep,7:]=usnd[istep,7:]
  vbg[istep,7:]=vsnd[istep,7:]
  qbg[istep,:]=qsnd[istep,:]
  tbg[istep,:]=tsnd[istep,:]  
  weight=[0.,0.,0.2,.4,.6,.8,1.]
  for ilev in range(0,7):
    ubg[istep,ilev]=(1.-weight[ilev])*ulsf_sonde[istep,ilev]+weight[ilev]*usnd[istep,ilev]
    vbg[istep,ilev]=(1.-weight[ilev])*vlsf_sonde[istep,ilev]+weight[ilev]*vsnd[istep,ilev]
#    qbg[istep,ilev]=(1.-weight[ilev])*qlsf_sonde[istep,ilev]+weight[ilev]*qsnd[istep,ilev]
#    tbg[istep,ilev]=(1.-weight[ilev])*tlsf_sonde[istep,ilev]+weight[ilev]*tsnd[istep,ilev]


if (lblend):
  for tidx in range(0,nsnd): 
    plt.figure(figsize=(14,6))
    plt.suptitle("sonde nr "+str(tidx))
    c=4
    r=1
    w=1
    #tidx=10
    ytop=500
    plt.subplot(r,c,w)
    plt.plot(tlsf_sondetime[tidx,:],zlsf_sondetime[tidx,:],label='lsf')
    plt.plot(tsnd[tidx,:],zsnd,color='red',label='sonde')
    plt.plot(tbg[tidx,:],zsnd,color='green',label='blend')
    plt.legend(loc='best')
    plt.ylim([0,ytop])
    plt.xlim([285,295])
    w=w+1

    plt.subplot(r,c,w)
    plt.plot(qlsf_sondetime[tidx,:],zlsf_sondetime[tidx,:],label='lsf')
    plt.plot(qsnd[tidx,:],zsnd,color='red',label='sonde')
    plt.plot(qbg[tidx,:],zsnd,color='green',label='blend')
    plt.legend(loc='best')
    plt.ylim([0,ytop])
    plt.xlim([0.005,.015])
    w=w+1

    plt.subplot(r,c,w)
    plt.plot(ulsf_sondetime[tidx,:],zlsf_sondetime[tidx,:],label='lsf')
    plt.plot(usnd[tidx,:],zsnd,color='red',label='sonde')
    plt.plot(ubg[tidx,:],zsnd,color='green',label='blend')

    plt.legend(loc='best')
    plt.ylim([0,ytop])
    plt.xlim([-10.,10])
    w=w+1

    plt.subplot(r,c,w)
    plt.plot(vlsf_sondetime[tidx,:],zlsf_sondetime[tidx,:],label='lsf')
    plt.plot(vsnd[tidx,:],zsnd,color='red',label='sonde')
    plt.plot(vbg[tidx,:],zsnd,color='green',label='blend')
    plt.legend(loc='best')
    plt.ylim([0,ytop])
    plt.xlim([-10.,0.])
    w=w+1
    plt.savefig('/hpc/uhome/mahlgrim/DEPHY-SCM/MAGIC/LEG15A/images/driver_DEF/merged_background_surface_'+str(tidx).zfill(2)+'.png')
    plt.close()
  #  plt.show()


    plt.figure(figsize=(14,6))
    plt.suptitle("sonde nr "+str(tidx))      

    c=4
    r=1
    w=1
  #  tidx=10
    ytop=20000
    plt.subplot(r,c,w)
    plt.plot(tlsf_sondetime[tidx,:],zlsf_sondetime[tidx,:],label='lsf')
    plt.plot(tsnd[tidx,:],zsnd,color='red',label='sonde')
    plt.plot(tbg[tidx,:],zsnd,color='green',label='blend')
    plt.legend(loc='best')
    plt.ylim([0,ytop])
    w=w+1

    plt.subplot(r,c,w)
    plt.plot(qlsf_sondetime[tidx,:],zlsf_sondetime[tidx,:],label='lsf')
    plt.plot(qsnd[tidx,:],zsnd,color='red',label='sonde')
    plt.plot(qbg[tidx,:],zsnd,color='green',label='blend')
    plt.legend(loc='best')
    plt.ylim([0,ytop])
    w=w+1

    plt.subplot(r,c,w)
    plt.plot(ulsf_sondetime[tidx,:],zlsf_sondetime[tidx,:],label='lsf')
    plt.plot(usnd[tidx,:],zsnd,color='red',label='sonde')
    plt.plot(ubg[tidx,:],zsnd,color='green',label='blend')

    plt.legend(loc='best')
    plt.ylim([0,ytop])
    w=w+1

    plt.subplot(r,c,w)
    plt.plot(vlsf_sondetime[tidx,:],zlsf_sondetime[tidx,:],label='lsf')
    plt.plot(vsnd[tidx,:],zsnd,color='red',label='sonde')
    plt.plot(vbg[tidx,:],zsnd,color='green',label='blend')
    plt.legend(loc='best')
    plt.ylim([0,ytop])
    w=w+1
    plt.savefig('/hpc/uhome/mahlgrim/DEPHY-SCM/MAGIC/LEG15A/images/driver_DEF/merged_background_all_'+str(tidx).zfill(2)+'.png')
    plt.close()
  #  plt.show()

  
# Winds
case.add_init_wind(u=ubg[0,:], v=vbg[0,:], lev=zsnd, levtype='altitude')

# Temperature
case.add_init_temp(tsnd[0,:], lev=zsnd, levtype='altitude')
#print(zsnd[0:10])
#print(tsnd[0:10])
#print(qsnd[0:10])
#print(usnd[0:10])
#print(vsnd[0:10])
# Potential temperature
case.add_init_theta(thetasnd[0,:], lev=zsnd, levtype='altitude')

# Water vapor mixing ratio
case.add_init_rv(qsnd[0,:], lev=zsnd, levtype='altitude')

# Pressure
case.add_init_pressure(psnd[0,:],lev=zsnd, levtype='altitude')


# Initial state from large scale forcing

# Surface pressure
ps0 = psnd[0,0]
case.add_init_ps(ps0)

################################################
# 3. Forcing
################################################

# Model-field-derived large-scale forcing
# Based on 0.5deg lat/lon gridded ECMWF data
# A smoothing (Gaussian kernel with standard
# deviation of 2deg) was applied to the model
# fields

# height on model levels (hybrid pressure coordinate,
# varying in time)

lat = lsf.lat.values[4:]
lon = lsf.lon.values[4:]

case.add_latitude(lat,time=lsftime[4:])
case.add_longitude(lon,time=lsftime[4:])

#  Geostrophic wind
ug=lsf.u_geo.values[4:,:]
vg=lsf.v_geo.values[4:,:]
case.add_geostrophic_wind(ug=ug,vg=vg,time=lsftime[4:],lev=zlsf_sondetime[0,:],height=zlsf_sondetime[0,:],levtype='altitude')

# large-scale horizontal advective tendencies
t_adv  = lsf.tls[4:,:] #K/s
rv_adv = lsf.qls[4:,:] #kg/kg/s
#case.add_rv_advection(rv_adv,time=lsftime[4:],lev=zlsf[-1,:],height=zlsf,levtype='altitude')
#case.add_temp_advection(t_adv,time=lsftime[4:],lev=zlsf[-1,:],height=zlsf,levtype='altitude',include_rad=False)

case.add_rv_advection(rv_adv,time=lsftime[4:],lev=zlsf_sondetime[0,:],height=zlsf_sondetime[0,:],levtype='altitude')
case.add_temp_advection(t_adv,time=lsftime[4:],lev=zlsf_sondetime[0,:],height=zlsf_sondetime[0,:],levtype='altitude',include_rad=False)

# Pressure levels 
#plsf=lsf.p[4:,:]
#case.add_pressure_forcing(plsf,time=lsftime[4:],lev=zlsf[-1,:],height=zlsf,levtype='altitude')

# Surface pressure
#ps=lsf.p_surf[4:]*100. #Pa
#case.add_surface_pressure_forcing(ps,time=lsftime[4:])
case.add_surface_pressure_forcing(psnd[:,0],time=sndtime)


#plt.figure(figsize=(14,6))
#plt.plot(sndtime,psnd[:,0]/100.)
#plt.plot(lsftime[4:],ps/100.,color='red')
#plt.show()


# Large-scale vertical velocity
w=lsf.wls[4:,:] #m/s
case.add_vertical_velocity(w=w,lev=zlsf_sondetime[0,:],height=zlsf_sondetime[0,:],time=lsftime[4:],levtype='altitude')


#case.add_wind_nudging(unudg=unudg,vnudg=vnudg,timescale=3600.*3.,p_nudging=110000.,time=timeForc,timeid='time',lev=levForc,levtype='pressure',levid='lev')
case.add_temp_nudging(tbg,timescale=3600.*3.,z_nudging=0.,time=sndtime,timeid='time',lev=zsnd,levtype='altitude')
case.add_qv_nudging(qbg,timescale=3600.*3.,z_nudging=0.,time=sndtime,timeid='time',lev=zsnd,levtype='altitude')

# Surface Forcing
sst     = sfc.sst[255:] #K
case.add_forcing_ts(sst, time=sfctime[255:])

################################################
# 4. Writing file
################################################

case.write('MAGIC_LEG15A_DEF_driver.nc')

if lverbose:
    case.info()

################################################
# 5. Ploting, if asked
################################################

if lplot:
    case.plot(rep_images='./images/driver_DEF/',timeunits='hours')
