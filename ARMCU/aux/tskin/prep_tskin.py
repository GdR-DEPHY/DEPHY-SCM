#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on 05 July 2021

@author: Romain Roehrig

Modification
"""

import netCDF4 as nc
import numpy as np


missing = 1.e20

# Reading original data downloaded 
# from https://adc.arm.gov/discovery/#/results/instrument_code::irt10m
with nc.Dataset('sgpirt10mC1.b1.19970621.000200.cdf') as f1:
    data1 = f1['sfc_ir_temp'][:]
    time1 = f1['time'][:]
    tunits1 = f1['time'].units
    date1 = nc.num2date(f1['time'][:], units=f1['time'].units)

with nc.Dataset('sgpirt10mC1.b1.19970622.000200.cdf') as f2:
    data2 = f2['sfc_ir_temp'][:]
    time2 = f2['time'][:]
    tunits2 = f2['time'].units
    date2 = nc.num2date(f2['time'][:], units=f2['time'].units)


nt1, = data1.shape
nt2, = data2.shape

# Original data have 1-minute time resolution
time = np.arange(0,60*60*24*2,60, dtype=np.float64)
nt, = time.shape
tunits = 'seconds since 1997-06-21 00:00:00'

time1_new = nc.date2num(date1, units=tunits)
time2_new = nc.date2num(date2, units=tunits)

# For the sake of simplicity, data are supposed to be valid over the minute 
# indicated by the time axis in original data.
timeBnds = np.zeros((nt,2), dtype=np.float64)

tskin = np.zeros((nt,), dtype=np.float32) + missing

it1 = 0
it2 = 0
for it in range(0,nt):
    timeBnds[it,0] = it*60
    timeBnds[it,1] = (it+1)*60
    if it1 < nt1 and time1_new[it1] == time[it]:
        tskin[it] = data1[it1]
        it1 += 1
    if it1 >= nt1 and it2 < nt2 and time2_new[it2] == time[it]:
        tskin[it] = data2[it2]
        it2 += 1

with nc.Dataset('tskin_SGP_C1_irt10m_19970621000030-19970622235930.nc','w') as g:
    g.createDimension('time', None)
    g.createDimension('bnds', 2)

    timeVar = g.createVariable('time', np.float64, ('time',))
    # Time stamps are supposed to be at the center of each time interval
    timeVar[:] = time[:] + 30. 
    timeVar.standard_name = 'time'
    timeVar.units = tunits
    timeVar.calendar = 'gregorian'
    timeVar.bounds = 'time_bnds'

    timeBndsVar = g.createVariable('time_bnds', np.float64, ('time','bnds'))
    timeBndsVar[:,:] = timeBnds 

    tskinVar = g.createVariable('tskin', np.float32, ('time',), fill_value=missing)
    tskinVar[:] = tskin
    tskinVar.standard_name = 'skin_temperature'
    tskinVar.units = 'K'
    

    g.comment = 'Surface Infrared temperature at SGP C1 (Central Facility) site. IRT10m measurements'
    g.reference = 'Atmospheric Radiation Measurement (ARM) user facility. 1996, updated hourly. Infrared Thermometer (IRT10M). 1996-04-16 to 2021-05-27, Southern Great Plains (SGP) Central Facility, Lamont, OK (C1). Compiled by V. Morris and J. Howie. ARM Data Center. Data set accessed 2021-05-29 at http://dx.doi.org/10.5439/1025203.'

# To compute hourly mean
# cdo hourmean tskin_SGP_C1_irt10m_19970621000030-19970622235930.nc tskin_SGP_C1_irt10m_19970621003000-19970622233000.nc


