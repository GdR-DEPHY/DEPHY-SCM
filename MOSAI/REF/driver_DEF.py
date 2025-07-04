#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 29 Mai 2024

@author: DEPHY team "atelier cas 1D"

Modifications
27/05/2025 Alice Maison : various options added, including geostrophic wind, wind advection and forcing smoothing
"""

import os, sys

import netCDF4 as nc
import numpy as np
import scipy as sp
import scipy.ndimage

import argparse
parser=argparse.ArgumentParser()

from dephycf.Case import Case
from dephycf.other_thermo import rm_from_hu, td_from_rm, hu_from_td
from dephycf.utils import kernel_mean_gauss

# function to pass True|False
def str2bool(v):
  if isinstance(v, bool):
    return v
  if v.lower() in ('yes', 'true', 't', 'y', '1'):
    return True
  elif v.lower() in ('no', 'false', 'f', 'n', '0'):
    return False
  else:
    raise argparse.ArgumentTypeError('Boolean value expected.')

################################################
# 0. General configuration of the present script
################################################

parser.add_argument("cover", type=str, metavar="cover", choices=['MAIZE','DECIDUOUS'], help="cover type: MAIZE|DECIDUOUS")
parser.add_argument("initpf", type=str, metavar="initpf", choices=['rs_smooth','rs_idea'], help="initial profile: rs_smooth|rs_idea")
parser.add_argument("--advTq", type=str2bool, nargs='?', const=True, default=False, help="flag to activate advection tendencies of temperature and humidity")
parser.add_argument("--advuv", type=str2bool, nargs='?', const=True, default=False, help="flag to activate advection tendencies of horizontal wind")
parser.add_argument("--vertvel", type=str2bool, nargs='?', const=True, default=False, help="flag to add vertical velocity")
parser.add_argument("fadv", type=str, metavar="adv_from", choices=['ARPEGEoper','ERA5'], help="advection from: ARPEGEoper|ERA5")
parser.add_argument("sadv", type=float, default=0., metavar="smooth_adv", help="smooth advection tendencies")
parser.add_argument("zadv", type=float, default=50000., metavar="max_alt_adv", help="maximum altitude of advection tendencies")
parser.add_argument("--geo", type=str2bool, nargs='?', const=True, default=False, help="flag to activate geostrophic wind")
parser.add_argument("sgeo", type=float, default=0., metavar="smooth_geos_wind", help="smooth geostrophic wind")
parser.add_argument("zgeo", type=float, default=16000., metavar="max_alt_geos_wind", help="maximum altitude of geostrophic wind")
parser.add_argument("--rad", type=str2bool, nargs='?', const=True, default=False, help="flag to activate radiation")
parser.add_argument("--ffx", type=str2bool, nargs='?', const=True, default=False, help="flag to activate surface flux forcing")
parser.add_argument("--fts", type=str2bool, nargs='?', const=True, default=False, help="flag to activate surface temperature forcing")
parser.add_argument("-s", type=float, default=-9999, metavar="thresh_sensib", help="sensible heat flux threshold")
args=parser.parse_args()

lplot = True # plot all the variables
lverbose = False # print information about variables and case

################################################
# 1. General information about the case
################################################

start_date = "20230819050000"
end_date  =  "20230820180000"
Zorog = 594.

cover = args.cover
initial_profile = args.initpf
advTq = args.advTq
advuv = args.advuv
vertvel = args.vertvel
adv_from = args.fadv
smooth_adv = args.sadv
zmax_adv = args.zadv
geoswd = args.geo
smooth_geo = args.sgeo
zmax_geo = args.zgeo
radia = args.rad
forc_flux = args.ffx
forc_ts = args.fts
sensib_thresh = args.s

# name of the case
add_nam = ''
if ((advTq) | (advuv) | (vertvel)):
  if (advTq): add_nam = add_nam+'_advTq'
  if (advuv): add_nam = add_nam+'_advuv'
  if (vertvel): add_nam = add_nam+'_w'
  add_nam = add_nam+'_'+adv_from[0:3]
if (geoswd): add_nam = add_nam+'_geo'
scase = cover+'_'+args.initpf+add_nam

print('cover:', cover)
print('initial profiles:', initial_profile)
print('case:', scase)
print('advection of temperature and humidity:', advTq)
print('advection of horizontal wind:', advuv)
print('vertical velocity:', vertvel)
print('from:', adv_from)
print('smooth adv:', smooth_adv)
print('Zmax adv:', zmax_adv)
print('geostrophic wind:', geoswd)
print('smooth geo wind:', smooth_geo)
print('Zmax geo wind:', zmax_geo)
print('radiation:', radia)
print('flux forcing:', forc_flux)
print('surface temp. forcing:', forc_ts)
print('sensible heat flux threshold:', sensib_thresh)

from datetime import datetime, timedelta
timeref = datetime.strptime(start_date, "%Y%m%d%H%M%S")
timeend = datetime.strptime(end_date, "%Y%m%d%H%M%S")

case = Case('MOSAI/%s'%scase,
        lat=43.1,
        lon=0.36,
        startDate=start_date,
        endDate=end_date,
        surfaceType='land',
        zorog=0.)
        #zorog=Zorog)

case.set_title("Forcing and initial conditions for MOSAI case - %s case"%scase)
case.set_reference("MOSAI campaign")
case.set_author("DEPHY team")
case.set_script("DEPHY-SCM/MOSAI/%s/driver_DEF.py"%scase)

################################################
# 2. Initial state
################################################
# Surface pressure
ps = 95260. # Approximate value from observations
case.add_init_ps(ps)

if ((initial_profile == 'rs_idea') | (initial_profile == 'rs_idea_ql')):
  ### initial profiles computed from observed data from radiosoundings, 60m-mast and 6m-mast at 5 am UTC
  z = np.genfromtxt('input.csv', delimiter=';', skip_header=1, usecols=0)
  temp = np.genfromtxt('input.csv', delimiter=';', skip_header=1, usecols=1)
  pre = np.genfromtxt('input.csv', delimiter=';', skip_header=1, usecols=2)
  hur = np.genfromtxt('input.csv', delimiter=';', skip_header=1, usecols=3)
  u = np.genfromtxt('input.csv', delimiter=';', skip_header=1, usecols=4)
  v = np.genfromtxt('input.csv', delimiter=';', skip_header=1, usecols=5)
  # in rm_from_hu, p is in hPa, T in K and h in %
  # returns g/kg
  rt = [0.001*rm_from_hu(p, t, h) for (p,t,h) in zip(pre/100.,temp,hur)]

elif (initial_profile == 'rs_smooth'):
  ds = nc.Dataset("input.nc", "r")
  # utility functions to read netCDF variables, smooth the soundings and keep
  # only the ascending part of the profiles
  def getvar(var):
    if (var == "height"):
      return ds.variables[var][:] - Zorog
    else:
      return ds.variables[var][:]
  def smoothvar(var, npts=40):
    return kernel_mean_gauss(getvar(var)[imin:imax],z,npts)

  z = getvar("height")    # m
  imin = 0
  zdiff = 0
  while zdiff < 2:
    imin+=1
    zdiff=z[imin]-z[imin-1]
  imax = np.argmax(z)     # to select only the upward sounding
  z = z[imin:imax]

  # get all variables and smooth them
  u = smoothvar("ua")     # m/s
  v = smoothvar("va")     # m/s
  temp = smoothvar("ta")  # °C
  hur = smoothvar("hur")  # %
  pre = smoothvar("pa")   # hPa
  z = smoothvar("height") # m
  temp = [t + 273.15 for t in temp] # => to K
  # in rm_from_hu, p is in hPa, T in K and h in %
  # returns g/kg
  rt = [0.001*rm_from_hu(p, t, h) for (p,t,h) in zip(pre,temp,hur)]
  pre = [p*100 for p in pre] # in common format, P is in Pa

else:
  print('Please choose between rs_idea and rs_smooth for initpf')
  exit()

# Wind initial profiles
case.add_init_wind(u=u, v=v, lev=z, levtype='altitude')

# Pressure
case.add_init_pressure(pre, lev=z, levtype='altitude')

# Temperature
case.add_init_temp(temp, lev=z, levtype='altitude')

# Relative humidity
hur = hur/100. # humidity between 0 and 1
case.add_init_hur(hur, lev=z, levtype='altitude')

# Total water mixing ratio
case.add_init_rt(rt, lev=z, levtype='altitude')

################################################
# 3. Forcing
################################################

# Geostrophic wind
if geoswd:
  # radar VHF between 1.5 and 16 km
  dst = nc.Dataset("geostrophic_wind.nc", "r")
  time_counter_vhf = np.array(dst.variables['time'][:])
  lev_vhf = np.array(dst.variables['level'][:]) # above mean sea level
  u_vhf = np.array(dst.variables['UWE'][:])
  v_vhf = np.array(dst.variables['VSN'][:])
  z_vhf = lev_vhf - 587.

  time_start_vhf = datetime(2023,8,1,0,0,0)
  Nt_vhf = time_counter_vhf.shape[0]
  time_vhf_ = np.array([time_start_vhf + timedelta(seconds = int(time_counter_vhf[t])) for t in range(Nt_vhf)])
  mask = np.array(np.where((time_vhf_ >= timeref) & (time_vhf_ <= timeend))).squeeze()
  sec = 800
  tmin = np.array(np.where((time_vhf_ >= timeref - timedelta(seconds = sec)) & (time_vhf_ <= timeref + timedelta(seconds = sec)))).min()
  tmax = np.array(np.where((time_vhf_ >= timeend - timedelta(seconds = sec)) & (time_vhf_ <= timeend + timedelta(seconds = sec)))).max()

  time_vhf = time_counter_vhf[tmin:tmax]
  time_vhf = time_vhf - time_vhf.min()
  u_vhf = u_vhf[tmin:tmax,:]
  v_vhf = v_vhf[tmin:tmax,:]

  # smooth
  sigma_x = smooth_geo
  sigma_y = smooth_geo
  sigma = [sigma_y, sigma_x]
  u_vhf = sp.ndimage.filters.gaussian_filter(u_vhf, sigma, mode='constant')
  v_vhf = sp.ndimage.filters.gaussian_filter(v_vhf, sigma, mode='constant')

  # assume constant wind over this altitude
  zmax_vhf = np.max(z_vhf[z_vhf<=zmax_geo])
  u_vhf[:,z_vhf>zmax_vhf] = u_vhf[:,z_vhf==zmax_vhf]
  v_vhf[:,z_vhf>zmax_vhf] = v_vhf[:,z_vhf==zmax_vhf]

  # the wind value at 1500m is extended to the surface
  u_vhf_zmin = np.reshape(u_vhf[:,0], (u_vhf.shape[0],1))
  v_vhf_zmin = np.reshape(v_vhf[:,0], (v_vhf.shape[0],1))
  z_vhf = np.insert(z_vhf, 0, 0., axis=0)
  u_vhf = np.concatenate((u_vhf_zmin, u_vhf), axis=1)
  v_vhf = np.concatenate((v_vhf_zmin, v_vhf), axis=1)

  # VHF radar
  case.add_geostrophic_wind(ug=u_vhf,vg=v_vhf,time=time_vhf,lev=z_vhf,levtype='altitude')

# Advection
if advTq | advuv | vertvel :

  # Advection tendencies from ARPEGE oper files
  # reverse altitude dimension to get height-increasing variables
  if (adv_from == 'ARPEGEoper'):
    fin = nc.Dataset("tendencies.nc", "r")
    timeForc= fin['time'][:]-fin['time'][0]
    levForc = fin['height_f'][0,::-1]
    levForc_h = fin['height_h'][0,::-1]
    #Compute and subtract orography altitude
    Zorog_adv = levForc_h[0] - (levForc[0] - levForc_h[0])
    levForc = levForc - Zorog_adv
    if advTq:
      tadv = fin['tadv'][:,::-1]
      qadv = fin['qadv'][:,::-1]
    if advuv:
      uadv = fin['uadv'][:,::-1]
      vadv = fin['vadv'][:,::-1]
    if vertvel:
      vtvl = fin['omega'][:,::-1]

  # Advection tendencies from ERA5 reanalysis
  elif (adv_from == 'ERA5'):
    fin = nc.Dataset("tendencies.nc", "r")
    timeForc = fin['time'][:]
    levForc = fin['z'][:,:]
    levForc = levForc.mean(axis=0)
    levForc = levForc - levForc[0]
    if advTq:
      tadv = fin['tadv'][:,:]
      qadv = fin['qadv'][:,:]
    if advuv:
      uadv = fin['uadv'][:,:]
      vadv = fin['vadv'][:,:]
    if vertvel:
      vtvl = fin['w'][:,:]

  else:
    print('Please choose between ARPEGEoper and ERA5 for fadv')
    exit()

  # Smooth advection
  sigma_x = smooth_adv
  sigma_y = smooth_adv
  sigma = [sigma_y, sigma_x]
  # Remove tendencies above this altitude
  levForc2 = levForc[levForc<=zmax_adv]
 
  # Temperature and humidity advection
  if advTq:
    tadv = sp.ndimage.filters.gaussian_filter(tadv, sigma, mode='constant')
    qadv = sp.ndimage.filters.gaussian_filter(qadv, sigma, mode='constant')
    tadv2 = tadv[:,levForc<=zmax_adv]
    qadv2 = qadv[:,levForc<=zmax_adv]
    case.add_temp_advection(tadv2,time=timeForc,timeid='time',lev=levForc2,levtype='altitude')
    case.add_qv_advection(qadv2,time=timeForc,timeid='time',lev=levForc2,levtype='altitude')

  if advuv:
    uadv = sp.ndimage.filters.gaussian_filter(uadv, sigma, mode='constant')
    vadv = sp.ndimage.filters.gaussian_filter(vadv, sigma, mode='constant')
    uadv2 = uadv[:,levForc<=zmax_adv]
    vadv2 = vadv[:,levForc<=zmax_adv]
    case.add_wind_advection(ua_adv=uadv2, va_adv=vadv2,time=timeForc,timeid='time',lev=levForc2,levtype='altitude')

  if vertvel:
    vtvl = sp.ndimage.filters.gaussian_filter(vtvl, sigma, mode='constant')
    vtvl2 = vtvl[:,levForc<=zmax_adv]
    if (adv_from == 'ARPEGEoper'):
      case.add_vertical_velocity(omega=vtvl2,time=timeForc,lev=levForc2,levtype='altitude')
    elif (adv_from == 'ERA5'):
      case.add_vertical_velocity(w=vtvl2,time=timeForc,lev=levForc2,levtype='altitude')

else:
  # No advection forcing => set to 0
  z_tend_adv = [0,15000]
  tadv = [0,0]
  qadv = [0,0]
  case.add_temp_advection(tadv, lev=z_tend_adv, levtype="altitude")
  case.add_qv_advection(qadv, lev=z_tend_adv, levtype="altitude")

# Radiation
if not(radia):
  case.deactivate_radiation()

if (forc_flux | forc_ts):
  # Surface Forcings : day time latent sensible ustar tskin
  flux_file = "flux.txt"
  dateSfc = np.genfromtxt(flux_file,dtype=str,skip_header=1,usecols=0)
  hourSfc = np.genfromtxt(flux_file,dtype=str,skip_header=1,usecols=1)
  lhf     = np.genfromtxt(flux_file,dtype=float,skip_header=1,usecols=2)
  shf     = np.genfromtxt(flux_file,dtype=float,skip_header=1,usecols=3)
  us      = np.genfromtxt(flux_file,dtype=float,skip_header=1,usecols=4)
  ts      = np.genfromtxt(flux_file,dtype=float,skip_header=1,usecols=5)

  time_from_start_sec = []
  latent_flux = []
  sensib_flux = []
  ustar = []
  surface_temp = []
  for it in range(0, dateSfc.shape[0]):
    datetime_str = dateSfc[it]+" "+hourSfc[it]
    datetime_object = datetime.strptime(datetime_str, '"%Y-%m-%d %H:%M:%S"')
    delta_t_sec = (datetime_object-timeref).total_seconds()
    if delta_t_sec >= 0:
      time_from_start_sec += [delta_t_sec]
      latent_flux += [lhf[it]]
      sensib_flux += [max(sensib_thresh,shf[it])]
      ustar += [us[it]]
      surface_temp += [ts[it]]

  if forc_flux:
    case.add_surface_fluxes(sens=sensib_flux,lat=latent_flux,time=time_from_start_sec,
            forc_wind='ustar',ustar=ustar,time_ustar=time_from_start_sec)
  if forc_ts:
    case.add_surface_skin_temp(surface_temp, time=time_from_start_sec)

################################################
# 4. Writing file
################################################

case.write('MOSAI_%s_DEF_driver.nc'%scase)

if lverbose:
    case.info()

################################################
# 5. Ploting, if asked
################################################

if lplot:
    case.plot(rep_images='./images/driver_DEF/',timeunits='hours')
