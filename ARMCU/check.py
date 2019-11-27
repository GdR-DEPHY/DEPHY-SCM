#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on 27 November 2019

@author: Romain Roehrig
"""

import os
import sys
sys.path = ['../utils/',] + sys.path

import numpy as np
import netCDF4 as nc

import SCM_utils as utils

rep_images = './images/check_1D/'

if not(os.path.exists(rep_images)):
    os.makedirs(rep_images)

data0 = {}

f = nc.Dataset('ARMCU_REF_orig.nc','r')

for var in f.variables:
    if not(var in f.dimensions) and (f[var].ndim <= 3 or f[var][:].shape[0] == 1):
        print var
        data0[var] = utils.read(var,f)

f.close()

data = {}

f = nc.Dataset('ARMCU_REF_1D.nc','r')

for var in f.variables:
    if not(var in f.dimensions) and (f[var].ndim <= 3  or f[var][:].shape[0] == 1):
        print var
        data[var] = utils.read(var,f)
        #data[var].info()
        if data0.has_key(var):
            data0[var].plot(rep_images=rep_images,var2=data[var],label="orig",label2="interpolated")

f.close()
