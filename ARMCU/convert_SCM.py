#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on 27 November 2019

@author: Romain Roehrig
"""

import os
import sys
sys.path = ['../utils/',] + sys.path

import time

import numpy as np
import netCDF4 as nc

import SCM_utils as utils
import thermo

levout = np.array(range(0,6001,10),dtype=np.float64) # New vertical grid, 10-m resolution from surface to 6000 m (above the surface)
levout = utils.Axis('lev',levout,name='Altitude',units='m')

tunits = 'seconds since 1997-06-21 11:30:0.0'
t0 = 0 # 11:30 UTC, 21 June 1997
t0 = utils.Axis('t0',[t0,],name='Initial time',units=tunits)

timeout = np.array(range(0,86400+2*3600+1-41400,1800),dtype=np.float64) # up to 02:00 UTC, 22 June 1997, 30-min timestep
timeout = utils.Axis('time',timeout,name='time',units=tunits)

datain = {}
dataout = {}

f = nc.Dataset('ARMCU_REF_orig.nc','r')

for var in f.variables:
    if not(var in f.dimensions):
        print var
        datain[var] = utils.read(var,f)
        #datain[var].info()
        
        if not(var in ['ps',]):
            #print 'interpo'
            dataout[var] = utils.interpol(datain[var],levout=levout,timeout=timeout)
            #dataout[var].info()
            #dataout[var].plot(rep_images=rep_images)
        else:
            dataout[var] = datain[var]


# get latitude and longitude axes
lat = datain['ps'].lat
lon = datain['ps'].lon

# get a few variables required for computing new variables
ps=datain['ps'].data[0,0,0]
z = dataout['height'].data[0,:,0,0]
theta = dataout['theta'].data[0,:,0,0]

# Compute pressure as a function of z and theta
pressure = thermo.z2p(theta=theta,z=z,ps=ps)
pressure = np.reshape(pressure,(1,levout.length,1,1))
dataout['pressure'] = utils.Variable('pressure',name='Pressure',units='Pa',data=pressure,level=levout,time=t0,lat=lat,lon=lon)

# Compute temperature from theta
temp = thermo.theta2t(p=pressure[0,:,0,0],theta=theta)
temp = np.reshape(temp,(1,levout.length,1,1))
dataout['temp'] = utils.Variable('temp',name='Temperature',units='K',data=temp,level=levout,time=t0,lat=lat,lon=lon)

# Convert mixing ratio to specific humidity
qt = thermo.rt2qt(dataout['rt'].data[0,:,0,0],units='kg kg-1')
qt = np.reshape(qt,(1,levout.length,1,1))
dataout['qv'] = utils.Variable('qv',name='Specific Humidity',units='kg kg-1',data=qt,level=levout,time=t0,lat=lat,lon=lon)

# Add initial profiles for ql and qi
ql = dataout['qv'].data*0.
dataout['ql'] = utils.Variable('ql',name='Liquid Water Content',units='kg kg-1',data=ql,level=levout,time=t0,lat=lat,lon=lon)
dataout['qi'] = utils.Variable('qi',name='Ice Water Content',   units='kg kg-1',data=ql,level=levout,time=t0,lat=lat,lon=lon)

# Add Forcing surface pressure, constant in time
nt,nlev,nlat,nlon = dataout['thadv'].data.shape
ps_forc = np.zeros((nt,nlat,nlon),dtype=np.float32)
for it in range(0,nt):
    ps_forc[it,0,0] = ps

dataout['ps_forc']   = utils.Variable('ps_forc',name='Forcing surface pressure',units='Pa',data=ps_forc,time=timeout,lat=lat,lon=lon)

# Add Forcing surface temperature, constant in time to be revisited...
ts = np.zeros((nt,nlat,nlon),dtype=np.float32) + 320.

dataout['ts']   = utils.Variable('ts',name='Surface temperature',units='K',data=ts,time=timeout,lat=lat,lon=lon)

# Add Forcing pressure and height, constant in time
height_forc = np.zeros((nt,nlev,nlat,nlon),dtype=np.float32)
pressure_forc = np.zeros((nt,nlev,nlat,nlon),dtype=np.float32)
for it in range(0,nt):
    height_forc[it,:,0,0] = dataout['height'].data[0,:,0,0]
    pressure_forc[it,:,0,0] = dataout['pressure'].data[0,:,0,0]

dataout['height_forc']   = utils.Variable('height_forc',  name='Forcing height',  units='m', data=height_forc,  level=levout,time=timeout,lat=lat,lon=lon)
dataout['pressure_forc'] = utils.Variable('pressure_forc',name='Forcing pressure',units='Pa',data=pressure_forc,level=levout,time=timeout,lat=lat,lon=lon)

# Compute temperature advection from theta advection
tadv = thermo.theta2t(p=dataout['pressure_forc'].data,theta=dataout['thadv'].data)
dataout['tadv'] = utils.Variable('tadv',name='Advection of temperature',units='K s-1',data=tadv,level=levout,time=timeout,lat=lat,lon=lon)

# Suppose that rv advection equals rt advection...
dataout['rvadv'] = utils.Variable('rvadv',name='Advection of water vapor mixing ratio',units='kg kg s-1',data=dataout['rtadv'].data,level=levout,time=timeout,lat=lat,lon=lon)


# Compute qt advection from rt advection and suppose qv advection equals qt advection.
qtadv =  thermo.advrt2advqt(rt=dataout['rt'].data,advrt=dataout['rtadv'].data)
dataout['qtadv'] = utils.Variable('qtadv',name='Advection of total specific humidity',units='kg kg s-1',data=qtadv,level=levout,time=timeout,lat=lat,lon=lon)
dataout['qvadv'] = utils.Variable('qvadv',name='Advection of specific humidity',units='kg kg s-1',data=qtadv,level=levout,time=timeout,lat=lat,lon=lon)



g = nc.Dataset('ARMCU_REF_1D.nc','w',format='NETCDF3_CLASSIC')

for var in ['ps','height','pressure','u','v','temp','theta','qv','rt','ql','qi','tke','ps_forc','height_forc','pressure_forc','ug','vg','tadv','thadv','qvadv','qtadv','rvadv','rtadv','ts','sfc_sens_flx','sfc_lat_flx']:
    print var
    #data[var].info()
    dataout[var].write(g)

g.Conventions = "CF-1.0"
g.comment = "Forcing and initial conditions for ARMCu case - Original definition - SCM-enabled version"
g.reference = f.reference
g.author = "R. Roehrig"
g.modifications = ""
g.script = 'DEPHY-SCM/ARMCU/convert_SCM.py'
g.history = "Created " + time.ctime(time.time())
g.case = f.case
g.startDate = f.startDate ;
g.endDate = f.endDate ;
g.tadv = 1
g.qvadv = 1
g.qtadv = 1
g.thadv = 1
g.rvadv = 1
g.rtadv = 1
g.trad = f.thrad 
g.thrad = f.thrad 
g.forc_omega = f.forc_omega 
g.forc_w = f.forc_w 
g.forc_geo = f.forc_geo 
g.nudging_u = f.nudging_u 
g.nudging_v = f.nudging_v 
g.nudging_t = f.nudging_th 
g.nudging_th = f.nudging_th 
g.nudging_qv = f.nudging_rt
g.nudging_qt = f.nudging_rt
g.nudging_rv = f.nudging_rt 
g.nudging_rt = f.nudging_rt 
g.zorog = f.zorog 
g.z0 = f.z0 
g.surfaceType = f.surfaceType 
g.surfaceForcing = f.surfaceForcing 
g.surfaceForcingWind = f.surfaceForcingWind 

g.close()

f.close()
