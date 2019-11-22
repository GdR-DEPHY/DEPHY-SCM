import os
import sys
sys.path.append('../utils/')

import numpy as np
import netCDF4 as nc
import SCM_utils as utils

import time

levout = np.array(range(0,6001,10),dtype=np.float64)
levout = utils.Axis('lev',levout,name='Altitude',units='m')

t0 = 41400
tunits = 'seconds since 1997-06-21 0:0:0.0'
t0 = utils.Axis('t0',[t0,],name='Initial time',units=tunits)

timeout = np.array(range(41400,93600+1,1800),dtype=np.float64)
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


lat = datain['ps'].lat
lon = datain['ps'].lon

ps=datain['ps'].data[0,0,0]
z = dataout['height'].data[0,:,0,0]
theta = dataout['theta'].data[0,:,0,0]

# Compute pressure as a function of z and theta
pressure = utils.z2p(theta=theta,z=z,ps=ps)
pressure = np.reshape(pressure,(1,levout.length,1,1))
dataout['pressure'] = utils.Variable('pressure',name='Pressure',units='Pa',data=pressure,level=levout,time=t0,lat=lat,lon=lon)

# Compute temperature from theta
temp = utils.theta2t(p=pressure[0,:,0,0],theta=theta)
temp = np.reshape(temp,(1,levout.length,1,1))
dataout['temp'] = utils.Variable('temp',name='Temperature',units='K',data=temp,level=levout,time=t0,lat=lat,lon=lon)

# Convert mixing ratio to specific humidity
qt = utils.rt2qt(dataout['rt'].data[0,:,0,0],units='kg kg-1')
qt = np.reshape(qt,(1,levout.length,1,1))
dataout['qv'] = utils.Variable('qv',name='Specific Humidity',units='kg kg-1',data=qt,level=levout,time=t0,lat=lat,lon=lon)

# Add initial profiles for ql and qi
ql = dataout['qv'].data*0.
dataout['ql'] = utils.Variable('ql',name='Liquid Water Content',units='kg kg-1',data=ql,level=levout,time=t0,lat=lat,lon=lon)
dataout['qi'] = utils.Variable('qi',name='Ice Water Content',   units='kg kg-1',data=ql,level=levout,time=t0,lat=lat,lon=lon)

# Add Forcing pressure and height, constant in time
nt,nlev,nlat,nlon = dataout['advth'].data.shape
height_forc = np.zeros((nt,nlev,nlat,nlon),dtype=np.float32)
pressure_forc = np.zeros((nt,nlev,nlat,nlon),dtype=np.float32)
for it in range(0,nt):
    height_forc[it,:,0,0] = dataout['height'].data[0,:,0,0]
    pressure_forc[it,:,0,0] = dataout['pressure'].data[0,:,0,0]

dataout['height_forc']   = utils.Variable('height_forc',  name='Forcing height',  units='m', data=height_forc,  level=levout,time=timeout,lat=lat,lon=lon)
dataout['pressure_forc'] = utils.Variable('pressure_forc',name='Forcing pressure',units='Pa',data=pressure_forc,level=levout,time=timeout,lat=lat,lon=lon)

# Compute temperature advection from theta advection
advt = utils.theta2t(p=dataout['pressure_forc'].data,theta=dataout['advth'].data)
dataout['advt'] = utils.Variable('advt',name='Advection of Temperature',units='K s-1',data=advt,level=levout,time=timeout,lat=lat,lon=lon)

# Suppose that qv advection equals rt advection...
dataout['advqv'] = utils.Variable('advqv',name='Advection of Specific humidity',units='kg kg s-1',data=dataout['advrt'].data,level=levout,time=timeout,lat=lat,lon=lon)



g = nc.Dataset('ARMCU_REF_1D.nc','w',format='NETCDF3_CLASSIC')

for var in ['ps','height','pressure','u','v','temp','theta','qv','rt','ql','qi','tke','height_forc','pressure_forc','ug','vg','advt','advqv','shf','lhf']:
    print var
    #data[var].info()
    dataout[var].write(g)

g.Conventions = "CF-1.0"
g.comment = "Forcing and initial conditions for ARMCu case - Original definition - SCM-enabled version"
g.reference = f.reference
g.author = "R. Roehrig"
g.modifications = ""
g.history = "Created " + time.ctime(time.time())
g.case = f.case
g.startDate = f.startDate ;
g.endDate = f.endDate ;
g.tadv = 1
g.qvadv = 1
g.tadvh = 0
g.qvadvh = 0
g.tadvv = 0
g.qvadvv = 0
g.trad = f.trad 
g.forc_omega = f.forc_omega 
g.forc_w = f.forc_w 
g.forc_geo = f.forc_geo 
g.nudging_u = f.nudging_u 
g.nudging_v = f.nudging_v 
g.nudging_t = f.nudging_t 
g.nudging_q = f.nudging_q 
g.zorog = f.zorog 
g.z0 = f.z0 
g.surfaceType = f.surfaceType 
g.surfaceForcing = f.surfaceForcing 
g.surfaceForcingWind = f.surfaceForcingWind 

g.close()

f.close()
