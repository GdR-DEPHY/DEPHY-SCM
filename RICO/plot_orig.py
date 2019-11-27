#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on 27 November 2019

@author: Romain Roehrig
"""

import os
import sys
sys.path = ['../utils/',] + sys.path

import netCDF4 as nc
import SCM_utils as utils

rep_images = './images/setup_orig/'

if not(os.path.exists(rep_images)):
    os.makedirs(rep_images)

data = {}

f = nc.Dataset('RICO_REF_orig.nc','r')

for var in f.variables:
    if not(var in f.dimensions):
        print var
        data[var] = utils.read(var,f)
        #data[var].info()
        data[var].plot(rep_images=rep_images)

f.close()
