#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on 5 Decemberr 2019

@author: Romain Roehrig
"""

## Check old and new forcing for RICO

import sys
sys.path = ['../../../utils/',] + sys.path

import numpy as np

from Case import Case

################################################
# 0. General configuration of the present script
################################################

lplot = True    # plot the new version of the case
lcompare = True # plot comparisons between original and new versions

################################################
# 1. Get the new version of the case
################################################

# initialize the case structure for the original version
newcase = Case('SANDU_REF')

# read case information in file
newcase.read('../SANDU_REF_SCM_driver.nc')

# display some information about the case
newcase.info()

################################################
# 2. Get the old version of the case
################################################

# initialize the case structure for the old version
oldcase = Case('SANDU/REF_OLD')

# read case information in file
oldcase.read('Composite_REF_driver_RR_new3_converted.nc')

# display some information about the case
oldcase.info()

################################################
# 4. Plots if asked
################################################

if lplot:
    newcase.plot(rep_images='./images/new_SCM/',timeunits='hours')
    oldcase.plot(rep_images='./images/old_SCM/',timeunits='hours')

if lcompare:
    newcase.plot_compare(oldcase,rep_images='./images/compare/',label1="New",label2="Old",timeunits='hours')

