## ARM-Cumulus original case definition
## From http://projects.knmi.nl/eurocs/ARM/case_ARM_html/

import os
import sys
sys.path.append('../utils/')
import time

import netCDF4 as nc
import numpy as np

import SCM_utils as utils

data = {}

###############################
# 1. General information about the case
###############################


lat = 36 #35.762
lon = -97.48

startDate = "19970621113000"
endDate =   "19970622020000"

tunits = 'seconds since 1997-06-21 11:30:0.0'
t0 = 0 # 11:30 UTC, 21 June 1997
t1 = 86400 + 2*3600 - 41400 # 02:00 UTC, 22 June 1997

ps = 97000. # Surface pressure (Pa)
zorog = 0 # Altitude above geoide (m)

z0 = 0.035 # Roughness length (m)

# Shared axes
t0Axis = utils.Axis('t0',[t0,],name='Initial time',units=tunits,calendar='gregorian')
latAxis = utils.Axis('lat',[lat,],name='Latitude',units='degrees_north')
lonAxis = utils.Axis('lon',[lon,],name='Longitude',units='degrees_east')


###############################
# 2. Initial state
###############################

data['ps'] = utils.Variable('ps',name='Surface Pressure',units='Pa',data=[ps,],time=t0Axis,lat=latAxis,lon=lonAxis)

#         z (m) theta (K) rt (g kg-1) u (m s-1) v (m s-1)
init = [   0.0,   299.00,   15.20,      10.0,     0.0,\
         50.0,   301.50,   15.17,      10.0,     0.0,\
        350.0,   302.50,   14.98,      10.0,     0.0,\
        650.0,   303.53,   14.80,      10.0,     0.0,\
        700.0,   303.70,   14.70,      10.0,     0.0,\
       1300.0,   307.13,   13.50,      10.0,     0.0,\
       2500.0,   314.00,    3.00,      10.0,     0.0,\
       5500.0,   343.20,    3.00,      10.0,     0.0]

init = np.array(init,dtype=np.float64)

z = init[0::5]
levAxes = {}
for var in ['height','theta','rt','u','v']:
    levAxes[var] = utils.Axis('lev_{0}'.format(var),z,name='Altitude for variable {0}'.format(var),units='m')

data['height'] = utils.Variable('height',name='Height',                     units='m',      data=init[0::5],      level=levAxes['height'],time=t0Axis,lat=latAxis,lon=lonAxis)
data['theta']  = utils.Variable('theta', name='Potential Temperature',      units='m',      data=init[1::5],      level=levAxes['theta'], time=t0Axis,lat=latAxis,lon=lonAxis)
data['rt']     = utils.Variable('rt',    name='Total Water Mixing Ratio',   units='kg kg-1',data=init[2::5]/1000.,level=levAxes['rt'],    time=t0Axis,lat=latAxis,lon=lonAxis)
data['u']      = utils.Variable('u',     name='Zonal Wind',                 units='m s-1',  data=init[3::5],      level=levAxes['u'],     time=t0Axis,lat=latAxis,lon=lonAxis)
data['v']      = utils.Variable('v',     name='Meridional Wind',            units='m s-1',  data=init[4::5],      level=levAxes['v'],     time=t0Axis,lat=latAxis,lon=lonAxis)

# Convert mixing ratio to specific humidity
#qt = utils.rt2qt(data['rt'].data,units='g kg-1')
#data['qt'] = utils.Variable('qt',name='Total Water Specific Mass',units='g kg-1',data=qt,level=z,levunits='m',time=[t0,],timeunits=tunits,timename='t0')

# Compute pressure as a function of z and theta
#pa = utils.z2p(theta=data['theta'].data,z=z,ps=ps)
#data['pa'] = utils.Variable('pa',name='Pressure',units='Pa',data=pa,level=z,levunits='m',time=[t0,],timeunits=tunits,timename='t0')

# Compute temperature from theta
#temp = utils.theta2t(p=data['pa'].data,theta=data['theta'].data)
#data['temp'] = utils.Variable('temp',name='Temperature',units='K',data=temp,level=z,levunits='m',time=[t0,],timeunits=tunits,timename='t0')

# Turbulent Kinetic Energy
ztke = range(0,6000+1,10)
nztke = len(ztke)
tke = np.zeros(nztke,dtype=np.float64)

for iz in range(0,nztke):
    if ztke[iz] < 150:
      tke[iz] = 0.15*(1.-ztke[iz]/150.)
    else:
      tke[iz] = 0.

levAxes['tke'] = utils.Axis('lev_tke',ztke,name='Altitude for variable tke',units='m')
data['tke'] = utils.Variable('tke',name='Turbulent Kinetic Energy',units='m2 s-2',data=tke,level=levAxes['tke'],time=t0Axis,lat=latAxis,lon=lonAxis)


###############################
# 3. Forcing
###############################

# Constant Geostrophic wind across the simulation
ug = np.zeros((2,8),dtype=np.float64)
ug[0,:] = init[3::5]
ug[1,:] = init[3::5]

vg = np.zeros((2,8),dtype=np.float64)
vg[0] = init[4::5]
vg[1] = init[4::5]

timeAxes = {}
for var in ['ug','vg']:
    timeAxes[var] = utils.Axis('time_{0}'.format(var),[t0,t1],name='Forcing time for {0}'.format(var),units=tunits)
    levAxes[var] = utils.Axis('lev_{0}'.format(var),z,name='Altitude for variable {0}'.format(var),units='m')

data['ug']     = utils.Variable('ug',    name='Geostrophic Zonal Wind',     units='m s-1', data=ug,level=levAxes['ug'],time=timeAxes['ug'],lat=latAxis,lon=lonAxis)
data['vg']     = utils.Variable('vg',    name='Geostrophic Meridional Wind',units='m s-1', data=vg,level=levAxes['vg'],time=timeAxes['vg'],lat=latAxis,lon=lonAxis)

# Surface Forcing
#            t (s) H (W m-2) LE (W m-2)
sfcForc= [  41400,  -30,       5,\
            55800,   90,     250,\
            64800,  140,     450,\
            68400,  140,     500,\
            77400,  100,     420,\
            86400,  -10,     180,\
            93600,  -10,       0]

sfcForc = np.array(sfcForc,dtype=np.float64)

timeSfc = sfcForc[0::3] - 41400
for var in ['sfc_sens_flx','sfc_lat_flx']:
    timeAxes[var] = utils.Axis('time_{0}'.format(var),timeSfc,name='Forcing time for {0}'.format(var),units=tunits)

data['sfc_sens_flx'] = utils.Variable('sfc_sens_flx',name='Surface Sensible Heat Flux',units='W m-2',data=sfcForc[1::3],time=timeAxes['sfc_sens_flx'],lat=latAxis,lon=lonAxis)
data['sfc_lat_flx'] = utils.Variable('sfc_lat_flx',name='Surface Latent Heat Flux',  units='W m-2',data=sfcForc[2::3],time=timeAxes['sfc_lat_flx'],lat=latAxis,lon=lonAxis)

# Advection forcing (+ radiative tendency)
#       t (s), A_theta (K hour-1) R_theta (K hour-1) A_rt (g kg-1 hour-1)
advF = [ 41400,      0.000,            -0.125,           0.080,\
         52200,      0.000,             0.000,           0.020,\
         63000,      0.000,             0.000,          -0.040,\
         73800,     -0.080,             0.000,          -0.100,\
         84600,     -0.160,             0.000,          -0.160,\
         93600,     -0.160,            -0.100,          -0.300]

advF = np.array(advF,dtype=np.float64)

timeF = advF[0::4] - 41400
ntf = len(timeF)
A_theta = advF[1::4]
R_theta = advF[2::4]
A_rt = advF[3::4]

zforc = range(0,6000+1,10)
nzf = len(zforc)
forc_theta = np.zeros((ntf,nzf),dtype=np.float64)
forc_rt = np.zeros((ntf,nzf),dtype=np.float64)

for var in ['thadv','rtadv']:
    timeAxes[var] = utils.Axis('time_{0}'.format(var),timeF,name='Forcing time for {0}'.format(var),units=tunits)
    levAxes[var] = utils.Axis('lev_{0}'.format(var),zforc,name='Altitude for {0}'.format(var),units='m')

# 2000 (Brown et al. 2002, QJRMS) ou 3000 m (http://projects.knmi.nl/eurocs/ARM/case_ARM_html/ and used in Meso-NH)
# We take 3000 m here.

for it in range(0,ntf):
    for iz in range(0,nzf):
        if zforc[iz] < 1000.:
          forc_theta[it,iz] = A_theta[it]+R_theta[it]
          forc_rt[it,iz] = A_rt[it]
        elif zforc[iz] <= 3000. :
          forc_theta[it,iz] = (A_theta[it]+R_theta[it])*(1.-(zforc[iz]-1000.)/2000.)
          forc_rt[it,iz] = A_rt[it]*(1.-(zforc[iz]-1000.)/2000.)
        else:
          forc_theta[it,iz] = 0.
          forc_rt[it,iz] = 0.

data['thadv'] = utils.Variable('thadv',name='Advection of Potential Temperature',   units='K s-1',      data=forc_theta/3600.,      level=levAxes['thadv'],time=timeAxes['thadv'],lat=latAxis,lon=lonAxis)
data['rtadv'] = utils.Variable('rtadv',name='Advection of Total Water Mixing Ratio',units='kg kg-1 s-1',data=forc_rt/1000./3600.,   level=levAxes['rtadv'],time=timeAxes['rtadv'],lat=latAxis,lon=lonAxis)


###############################
# 4. Writing file
###############################


g = nc.Dataset('ARMCU_REF_orig.nc','w',format='NETCDF3_CLASSIC')

for var in ['ps','height','u','v','theta','rt','tke','ug','vg','thadv','rtadv','sfc_sens_flx','sfc_lat_flx']: #data.keys():
    print var
    #data[var].info()
    #data[var].plot(rep_images=rep_images)
    data[var].write(g)

g.Conventions = "CF-1.0" 
g.comment = "Forcing and initial conditions for ARM-Cumulus case - Original definition" 
g.reference = "http://projects.knmi.nl/eurocs/ARM/case_ARM_html ; Brown et al. (2002, QJRMS)" 
g.author = "R. Roehrig" 
g.modifications = ""
g.script = 'DEPHY-SCM/ARMCU/setup_orig.py'
g.history = "Created " + time.ctime(time.time())
g.case = "ARMCU/REF" 
g.startDate = startDate
g.endDate = endDate
g.rtadv = 1 
g.thadv = 1 
g.thrad = "adv"
g.forc_omega = 0
g.forc_w = 0
g.forc_geo = 1
g.nudging_u = 0
g.nudging_v = 0
g.nudging_th = 0
g.nudging_rt = 0
g.zorog = zorog
g.z0 = z0
g.surfaceType = "land"
g.surfaceForcing = "surfaceFlux"
g.surfaceForcingWind = "z0"

g.close()


