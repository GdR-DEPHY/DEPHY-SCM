#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on 27 November 2019

@author: Romain Roehrig
"""

import os
import sys
sys.path = ['../../../utils/',] + sys.path

import numpy as np
import netCDF4 as nc

from Variable import read as readvar

rep_images = './images/setup_old/'

if not(os.path.exists(rep_images)):
    os.makedirs(rep_images)

data = {}

f = nc.Dataset('ARMCu_driver_RR_new3.nc','r')

for var in f.variables:
    if not(var in f.dimensions) and not(var in ['bounds_lat','bounds_lon']):
        print var
        data[var] = readvar(var,f)
        #data[var].info()
        data[var].plot(rep_images=rep_images)

f.close()
