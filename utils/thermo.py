#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on 27 November 2019

@author: Romain Roehrig
"""

import sys

import math

import numpy as np

import constants as cc

def rt2qt(rt,units='kg kg-1'):

    if units == 'kg kg-1':
        return rt/(1+rt)
    elif units == 'g kg-1':
        return rt/(1.+rt/1000.)
    else:
        print 'units unknown:', units
        sys.exit()


def advrt2advqt(rt=None,advrt=None,rt_units='kg kg-1'):

    if rt is None:
        print "rt is missing"
        sys.exit()
    if advrt is None:
        print "advrt is missing"
        sys.exit()

    if rt_units == 'kg kg-1':
        return advrt/((1+rt)*(1+rt))
    elif rt_units == 'g kg-1':
        return advrt/((1.+rt/1000.)*(1.+rt/1000.))
    else:
        print 'units unknown for rt:', rt_units
        sys.exit()


def theta2t(p=None,theta=None,p0=cc.p0,kappa=cc.kappa):

    if theta is None:
        print "theta is missing"
        sys.exit()
    if p is None:
        print "p is missing"
        sys.exit()

    tmp = theta*(p/p0)**kappa
    return tmp

def t2theta(p,t,p0=cc.p0,kappa=cc.kappa):

    if t is None:
        print "t is missing"
        sys.exit()
    if p is None:
        print "p is missing"
        sys.exit()

    tmp = t*(p0/p)**kappa
    return tmp


def z2p(theta=None,z=None,ps=None,g=cc.g,Rd=cc.Rd,p0=cc.p0,kappa=cc.kappa):

    if theta is None:
        print "theta is missing"
        sys.exit()
    if z is None:
        print "z is missing"
        sys.exit()
    if ps is None:
        print "ps is missing"
        sys.exit()

    nlev, = theta.shape

    integ = 0.

    p = np.zeros((nlev,),dtype=np.float64)
    p[0] = ps

    for ilev in range(1,nlev):
      dz = z[ilev]- z[ilev-1]    
#      print 'dz =', dz
      integ = integ + (g/(Rd*theta[ilev-1])+g/(Rd*theta[ilev]))/2*dz
#      print 'integ =', integ
      tmp = ps**kappa-p0**kappa*kappa*integ
      p[ilev] = math.exp(math.log(tmp)/kappa)

    return p

