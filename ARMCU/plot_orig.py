import os
import sys
sys.path.append('../utils/')

import numpy as np
import netCDF4 as nc
import SCM_utils as utils

import time

rep_images = './images/setup_orig/'

if not(os.path.exists(rep_images)):
    os.makedirs(rep_images)

data = {}

f = nc.Dataset('ARMCU_REF_orig.nc','r')

for var in f.variables:
    if not(var in f.dimensions):
        print var
        data[var] = utils.read(var,f)
        #data[var].info()
        data[var].plot(rep_images=rep_images)

f.close()
