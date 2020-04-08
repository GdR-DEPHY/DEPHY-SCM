#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on 7 April 2020

@author: Romain Roehrig
"""

import os
dir_path = os.path.dirname(os.path.realpath(__file__))

import netCDF4 as nc

import numpy as np
from scipy import interpolate

levmin = 49.

fERAI = nc.Dataset('{0}/dynamo_nsa3_ERAI_6hourly.nc'.format(dir_path))

levERAI = fERAI['level'][:]*100.

indz = levERAI < 5000.  # CSU data does not provide data above 50 hPa

nlevERAI, = levERAI[indz].shape

timeERAI = nc.num2date(fERAI['time'][:],fERAI['time'].units,calendar=fERAI['time'].calendar)
timeERAI = nc.date2num(timeERAI,units='seconds since 2011-10-01 0:0:0.0',calendar='gregorian')

fERAI.close()

fStd = nc.Dataset('{0}/Standard_Atmosphere.nc'.format(dir_path))

paStd = fStd['pa'][:]
ind = paStd < 100. # ERA-Interim does not provide data above 1 hPa
paStd = fStd['pa'][ind]
zgStd = fStd['zg'][ind]
taStd = fStd['ta'][ind]

nlevStd, = paStd.shape
fStd.close()

def extend(var,data,init=True,time=None):

    fERAI = nc.Dataset('{0}/dynamo_nsa3_ERAI_6hourly.nc'.format(dir_path))

    if init and (var == 'p'):

        nlev, = data.shape
        tmp = np.zeros(nlev+nlevERAI+nlevStd)
        tmp[0:nlev] = data
        tmp[nlev:nlev+nlevERAI] = levERAI[indz][::-1]
        tmp[nlev+nlevERAI:] = paStd

    elif init and (var in ['ta','hus','ua','va','wap','zg']):
        
        nlev, = data.shape
        tmp = np.zeros(nlev+nlevERAI+nlevStd)
        tmp[0:nlev] = data
        tmp[nlev:nlev+nlevERAI] = fERAI[var][0,:][indz][::-1]
        if var == 'zg':
            tmp[nlev:nlev+nlevERAI] = tmp[nlev:nlev+nlevERAI]/9.80665
            tmp[nlev+nlevERAI:] = zgStd[:]

        if var == 'ta':
            tmp[nlev+nlevERAI:] = taStd[:]

    elif var == 'p':

        nt,nlev = data.shape
        tmp = np.zeros((nt,nlev+nlevERAI+nlevStd))
        tmp[:,0:nlev] = data  
        for it in range(0,nt):
            tmp[it,nlev:nlev+nlevERAI] = levERAI[indz][::-1]
            tmp[it,nlev+nlevERAI:] = paStd


    elif var in ['ta','hus','ua','va','wap','zg']:

        nt,nlev = data.shape
        tmp = np.zeros((nt,nlev+nlevERAI+nlevStd))
        tmp[:,0:nlev] = data        

        dataERAI = fERAI[var][:,indz][:,::-1]

        for i in range(0,nlevERAI):
            ff = interpolate.interp1d(timeERAI,dataERAI[:,i])
            tmp[:,nlev+i] = ff(time)

        if var == 'zg':
            tmp[:,nlev:nlev+nlevERAI] = tmp[:,nlev:nlev+nlevERAI]/9.80665
            for it in range(0,nt):
                tmp[it,nlev+nlevERAI:] = zgStd[:]

        if var == 'ta':
            for it in range(0,nt):
                tmp[it,nlev+nlevERAI:] = taStd[:]


    else:
        nt,nlev = data.shape
        tmp = np.zeros((nt,nlev+nlevERAI+nlevStd))
        tmp[:,0:nlev] = data  

    fERAI.close()

    return tmp




