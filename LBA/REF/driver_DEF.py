<<<<<<< HEAD
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 11 December 2019

@author: Romain Roehrig

Modification
  2023/05/24, F. Couvreux: LBA case
"""

## AMMA/REF original case definition

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


case = Case('LBA/REF',
        lat=-8.,
        lon=-63.,
        startDate="19990223073000",
        endDate="19990223150000",
        surfaceType='land',
        zorog=0.)


case.set_title("Forcing and initial conditions for LBA case - Original definition")
case.set_reference(" https://rmets.onlinelibrary.wiley.com/doi/full/10.1256/qj.04.147; Grabowski et al. (2006, QJRMS)")
case.set_author("F. Couvreux")
case.set_script("DEPHY-SCM/LBA/REF/driver_DEF.py")

# time units are expected to be seconds since startDate
#t0 = 0 # 18:00 UTC, 15 July 2006
#t1 = 72*3600 # 72-hour long simulation

################################################
# 2. Input netCDF file
################################################


################################################
# 2. Initial state
################################################

# Surface pressure
ps = 99130.
case.add_init_ps(ps)

#initial profile in init.txt
#Altitude above the ground
z = np.genfromtxt('init.txt',dtype=None,skip_header=1,usecols=0)

# Zonal and meridional wind
u = np.genfromtxt('init.txt',dtype=None,skip_header=1,usecols=3)
v = np.genfromtxt('init.txt',dtype=None,skip_header=1,usecols=4)


case.add_init_wind(u=u, v=v, lev=z, levtype='altitude')

# Potential temperature
theta = np.genfromtxt('init.txt',dtype=None,skip_header=1,usecols=1)

case.add_init_theta(theta, lev=z, levtype='altitude')

# Water vapor mixing ratio
rv = np.genfromtxt('init.txt',dtype=None,skip_header=1,usecols=2)

case.add_init_rv(rv, lev=z, levtype='altitude')

################################################
# 3. Forcing
################################################

# Forcing time axis
timeForc = [0.,3600.,7200.,10800.,14400.,18000.,21600.]
height=np.genfromtxt('LBA_formatcommun_ZFR_1',dtype=None,usecols=0)
nk,=height.shape
u1=np.genfromtxt('LBA_formatcommun_ZFR_1',dtype=None,usecols=1)
v1=np.genfromtxt('LBA_formatcommun_ZFR_1',dtype=None,usecols=2)
dthdt1=np.genfromtxt('LBA_formatcommun_ZFR_1',dtype=None,usecols=6)
u2=np.genfromtxt('LBA_formatcommun_ZFR_2',dtype=None,usecols=1)
v2=np.genfromtxt('LBA_formatcommun_ZFR_2',dtype=None,usecols=2)
dthdt2=np.genfromtxt('LBA_formatcommun_ZFR_2',dtype=None,usecols=6)
u3=np.genfromtxt('LBA_formatcommun_ZFR_3',dtype=None,usecols=1)
v3=np.genfromtxt('LBA_formatcommun_ZFR_3',dtype=None,usecols=2)
dthdt3=np.genfromtxt('LBA_formatcommun_ZFR_3',dtype=None,usecols=6)
u4=np.genfromtxt('LBA_formatcommun_ZFR_4',dtype=None,usecols=1)
v4=np.genfromtxt('LBA_formatcommun_ZFR_4',dtype=None,usecols=2)
dthdt4=np.genfromtxt('LBA_formatcommun_ZFR_4',dtype=None,usecols=6)
u5=np.genfromtxt('LBA_formatcommun_ZFR_5',dtype=None,usecols=1)
v5=np.genfromtxt('LBA_formatcommun_ZFR_5',dtype=None,usecols=2)
dthdt5=np.genfromtxt('LBA_formatcommun_ZFR_5',dtype=None,usecols=6)
u6=np.genfromtxt('LBA_formatcommun_ZFR_6',dtype=None,usecols=1)
v6=np.genfromtxt('LBA_formatcommun_ZFR_6',dtype=None,usecols=2)
dthdt6=np.genfromtxt('LBA_formatcommun_ZFR_6',dtype=None,usecols=6)
u7=np.genfromtxt('LBA_formatcommun_ZFR_7',dtype=None,usecols=1)
v7=np.genfromtxt('LBA_formatcommun_ZFR_7',dtype=None,usecols=2)
dthdt7=np.genfromtxt('LBA_formatcommun_ZFR_7',dtype=None,usecols=6)
u=np.zeros((7,nk))
v=np.zeros((7,nk))
dthdt=np.zeros((7,nk))
u[0,:]=u1
u[1,:]=u2
u[2,:]=u3
u[3,:]=u4
u[4,:]=u5
u[5,:]=u6
u[6,:]=u7
v[0,:]=v1
v[1,:]=v2
v[2,:]=v3
v[3,:]=v4
v[4,:]=v5
v[5,:]=v6
v[6,:]=v7
dthdt[0,:]=dthdt1
dthdt[1,:]=dthdt2
dthdt[2,:]=dthdt3
dthdt[3,:]=dthdt4
dthdt[4,:]=dthdt5
dthdt[5,:]=dthdt6
dthdt[6,:]=dthdt7


#temperature advection and radiation
case.add_theta_advection(dthdt,include_rad=True,lev=height,levtype='altitude',time=timeForc)

#wind nudging
case.add_wind_nudging(unudg=u,vnudg=v,timescale=3600.,time=timeForc,timeid='time',lev=height,levtype='altitude',levid='lev')

# Surface fluxes
timeflux=np.genfromtxt('flux_LBA',dtype=None,usecols=0)
sens=np.genfromtxt('flux_LBA',dtype=None,usecols=1)
flat=np.genfromtxt('flux_LBA',dtype=None,usecols=2)
case.add_surface_fluxes(sens,flat,time=timeflux,forc_wind='z0',z0=0.035)

################################################
# 4. Writing file
################################################

case.write('LBA_REF_DEF_driver.nc')

if lverbose:
    case.info()

################################################
# 6. Ploting, if asked
################################################

if lplot:
    case.plot(rep_images='./images/driver_DEF/',timeunits='hours')
=======
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 24 May 2023

@author: Fleur Couvreux

Modification
  2023/06/05, R. Roehrig, some fixes and cleaning
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


case = Case('LBA/REF',
        lat=-8.,
        lon=-63.,
        startDate="19990223073000",
        endDate="19990223150000",
        surfaceType='land',
        zorog=0.)


case.set_title("Forcing and initial conditions for LBA case - Original definition")
case.set_reference(" https://rmets.onlinelibrary.wiley.com/doi/full/10.1256/qj.04.147; Grabowski et al. (2006, QJRMS)")
case.set_author("F. Couvreux")
case.set_script("DEPHY-SCM/LBA/REF/driver_DEF.py")

################################################
# 2. Initial state
################################################

# Surface pressure
ps = 99130.
case.add_init_ps(ps)

# Initial profile in init.txt
# Altitude above the ground
z = np.genfromtxt('init.txt',dtype=np.float32,usecols=0)

# Zonal and meridional wind
u = np.genfromtxt('init.txt',dtype=np.float32,usecols=3)
v = np.genfromtxt('init.txt',dtype=np.float32,usecols=4)

case.add_init_wind(u=u, v=v, lev=z, levtype='altitude')

# Potential temperature
theta = np.genfromtxt('init.txt',dtype=np.float32,usecols=1)

case.add_init_theta(theta, lev=z, levtype='altitude')

# Water vapor mixing ratio
rv = np.genfromtxt('init.txt',dtype=np.float32,usecols=2)

case.add_init_rv(rv, lev=z, levtype='altitude')

################################################
# 3. Forcing
################################################

# Forcing time axis
timeForc = [0.,3600.,7200.,10800.,14400.,18000.,21600.]
nt = len(timeForc)

height_forc = np.genfromtxt('LBA_formatcommun_ZFR_1',dtype=np.float32,usecols=0)
nlev, = height_forc.shape

data = {}
for var in ['u','v','dthdt']:
    data[var] = np.zeros((nt, nlev), dtype=np.float32)

for it in range(0,nt):
    fin = f'LBA_formatcommun_ZFR_{it+1}'
    data['u'][it,:] = np.genfromtxt(fin,dtype=np.float32,usecols=1)
    data['v'][it,:] = np.genfromtxt(fin,dtype=np.float32,usecols=2)
    data['dthdt'][it,:] = np.genfromtxt(fin,dtype=np.float32,usecols=6)

# Potential temperature advection, which includes radiation
case.add_theta_advection(data['dthdt'],include_rad=True,lev=height_forc,levtype='altitude',time=timeForc)

# Wind nudging
case.add_wind_nudging(unudg=data['u'],vnudg=data['v'],timescale=3600.,time=timeForc,lev=height_forc,levtype='altitude')

# Surface fluxes
timeflux = np.genfromtxt('flux_LBA',dtype=np.float32,usecols=0)
sens = np.genfromtxt('flux_LBA',dtype=np.float32,usecols=1)
flat = np.genfromtxt('flux_LBA',dtype=np.float32,usecols=2)
case.add_surface_fluxes(sens,flat,time=timeflux,forc_wind='z0',z0=0.035)

################################################
# 4. Writing file
################################################

case.write('LBA_REF_DEF_driver.nc')

if lverbose:
    case.info()

################################################
# 6. Ploting, if asked
################################################

if lplot:
    case.plot(rep_images='./images/driver_DEF/',timeunits='hours')
>>>>>>> 390ecda956c8fb0ac9d37e29161a0d2e7c3c3634
