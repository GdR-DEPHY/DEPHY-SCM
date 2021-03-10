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

def qt2rt(qt,units='kg kg-1'):

    if units == 'kg kg-1':
        return qt/(1-qt)
    elif units == 'g kg-1':
        return qt/(1.-qt/1000.)
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

def advqt2advrt(qt=None,advqt=None,qt_units='kg kg-1'):

    if qt is None:
        print "qt is missing"
        sys.exit()
    if advqt is None:
        print "advqt is missing"
        sys.exit()

    if qt_units == 'kg kg-1':
        return advqt/((1-qt)*(1-qt))
    elif qt_units == 'g kg-1':
        return advqt/((1.-qt/1000.)*(1.-qt/1000.))
    else:
        print 'units unknown for qt:', qt_units
        sys.exit()

def theta2t(p=None,theta=None,p0=cc.p0,kappa=cc.kappa):

    if theta is None:
        print "theta is missing"
        sys.exit()
    if p is None:
        print "p is missing"
        sys.exit()

    temp = theta*(p/p0)**kappa
    return temp

def t2theta(p=None,temp=None,p0=cc.p0,kappa=cc.kappa):

    if temp is None:
        print "t is missing"
        sys.exit()
    if p is None:
        print "p is missing"
        sys.exit()

    theta = temp*(p0/p)**kappa
    return theta


def z2p(thetal=None, theta=None, ta=None,
        z=None, ps=None, qv=None,
        g=cc.g, Rd=cc.Rd, Rv=cc.Rv, p0=cc.p0, kappa=cc.kappa):

    if (thetal is None) and (theta is None) and (ta is None):
        print "thetal, theta and ta are missing. At least one of them should be given"
        sys.exit()
    if z is None:
        print "z is missing"
        sys.exit()
    if ps is None:
        print "ps is missing"
        sys.exit()

    if thetal is not None and (theta, ta) is (None, None):
        theta = thetal

    if theta is not None:
        nlev, = theta.shape

        if qv is None:
            R = theta*0.+Rd
        else:
            R = Rd+qv*(Rv-Rd)

        integ = 0.

        p = np.zeros((nlev,),dtype=np.float64)
        p[0] = ps

        for ilev in range(1,nlev):
            dz = z[ilev]- z[ilev-1]    
#            print 'dz =', dz
            integ = integ + (g/(R[ilev-1]*theta[ilev-1])+g/(R[ilev]*theta[ilev]))/2*dz
#            print 'integ =', integ
            tmp = ps**kappa-p0**kappa*kappa*integ
            p[ilev] = math.exp(math.log(tmp)/kappa)
    else: # Use ta instead
        nlev, = ta.shape

        if qv is None:
            R = ta*0.+Rd
        else:
            R = Rd+qv*(Rv-Rd)

        integ = 0.

        p = np.zeros((nlev,),dtype=np.float64)
        p[0] = ps

        for ilev in range(1,nlev):
            dz = z[ilev]- z[ilev-1]    
#            print 'dz =', dz
            integ = integ + (g/(R[ilev-1]*ta[ilev-1])+g/(R[ilev]*ta[ilev]))/2*dz
#            print 'integ =', integ
            tmp = math.log(ps)-integ
            p[ilev] = math.exp(tmp)


    return p

def p2z(thetal=None, theta=None, ta=None,
        p=None, qv=None,
        g=cc.g, Rd=cc.Rd, Rv=cc.Rv, p0=cc.p0, kappa=cc.kappa):

    if (thetal is None) and (theta is None) and (ta is None):
        print "thetal, theta and ta are missing. At least one of them should be given"
        sys.exit()
    if p is None:
        print "p is missing"
        sys.exit()

    if theta is not None and ta is None:
        ta = theta2t(theta=theta,p=p)

    nlev, = ta.shape

    if qv is None:
        R = ta*0.+Rd
    else:
        R = Rd+qv*(Rv-Rd)

    z = np.zeros((nlev,),dtype=np.float64)
    z[0] = 0. 

    for ilev in range(1,nlev):
        dz = (R[ilev-1]*ta[ilev-1])+(R[ilev]*ta[ilev])/(2*g)*(math.log(p[ilev-1])-math.log(p[ilev]))
        z[ilev] = z[ilev-1] + dz

    return z

def zlev2plev(zlev,z,p):
    nlev, = z.shape
    plev = 110000
    if z[0] < z[1]:
        for ilev in range(0,nlev-1):
            if z[ilev] <= zlev and z[ilev+1] > zlev:
                plev = p[ilev] + (p[ilev+1]-p[ilev])/(z[ilev+1]-z[ilev])*(zlev-z[ilev])
    else:
        for ilev in range(0,nlev-1):
            if z[ilev+1] <= zlev and z[ilev] > zlev:
                plev = p[ilev] + (p[ilev+1]-p[ilev])/(z[ilev+1]-z[ilev])*(zlev-z[ilev])

    return plev



def plev2zlev(plev,z,p):
    nlev, = z.shape
    zlev = 0
    if p[0] > p[1]:
        for ilev in range(0,nlev-1):
            if p[ilev] >= plev and p[ilev+1] < plev:
                zlev = z[ilev] + (z[ilev+1]-z[ilev])/(p[ilev+1]-p[ilev])*(plev-p[ilev])
    else:
        for ilev in range(0,nlev-1):
            if p[ilev+1] >= plev and p[ilev] < plev:
                zlev = z[ilev] + (z[ilev+1]-z[ilev])/(p[ilev+1]-p[ilev])*(plev-p[ilev])

    return zlev

def rh2qv(rh,temp,pres):
    """
    Compute qv from relative humidity (rh in %) and temperature (in K) and pressure (in Pa)
    Based on Goff-Gratch equation
    From http://climatologie.u-bourgogne.fr/data/matlab/goff_gratch.m
    """

    # Saturation water vapor pressure against ice
    eilog = -9.09718 * ((273.16/temp) -1.)
    eilog2 = -3.5654 * np.log10(273.16/temp)
    eilog3 = 0.876793 * (1. - (temp/273.16))
    es1=6.1071*np.exp((eilog+eilog2+eilog3)*math.log(10.))

    # against liquid water
    eilog=-7.90298*((373.16/temp) - 1.)
    eilog2=5.02808*np.log10((373.16/temp))
    eilog3=-1.3816e-7*(np.exp((11.344*(1.-(temp/373.16)))*math.log(10.)) -1.)
    eilog4=8.1328e-3*(np.exp((-3.49149*((373.16/temp) - 1.) )*math.log(10)) -1.)
    es2=1013.246*np.exp((eilog+eilog2+eilog3+eilog4)*math.log(10.))

    es = np.where(temp < 273.15, es1,es2)

    ws = 0.62197* es/(pres/100. - 0.378*es) # qsat

    # specific humidity in kg kg-1
    return (rh/100.0)*ws
