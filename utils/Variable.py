#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on 27 November 2019

@author: Romain Roehrig
"""

import os
import sys
import time

import netCDF4 as nc

import numpy as np
from scipy import interpolate

from Axis import Axis

import plotbasics


class Variable:

    def __init__(self,varid,data=None,units=None,name=None,level=None,time=None,lat=None,lon=None,axlist=None,axes=None,plotcoef=1.,plotunits=None):

        self.id = varid
        self.units = units
        self.name = name

        self.data = None
        if not(data is None): 
            self.data = np.array(data,dtype=np.float32)
            self.sh = len(self.data.shape)

        if axes is None:
            self.axes = []
            self.axlist = []

            self.time = time
            self.level = level
            self.lat = lat
            self.lon = lon

            for ax in ['time','level','lat','lon']:
                if not(self.__dict__[ax] is None):
                    self.axes.append(self.__dict__[ax])
                    self.axlist.append(self.__dict__[ax].id)
        else:
            self.axes = axes
            self.axlist = axlist

            self.time = None
            self.level = None
            self.lat = None
            self.lon = None

            for ax in axes:
                if ax.id == 't0' or ax.id[0:4] == 'time':
                    self.time = ax
                elif ax.id[0:3] == 'lev' or ax.id == 'nlev':
                    self.level = ax
                elif ax.id == 'lat':
                    self.lat = ax
                elif ax.id == 'lon':
                    self.lon = ax
                else:
                    print 'Axis unexpected:', ax.id
                    ax.info()
                    sys.exit()

        self.coord = " ".join(self.axlist + ['lat','lon'])

        self.plotcoef = plotcoef
        if plotunits is None:
            self.plotunits = self.units
        else:
            self.plotunits = plotunits

    def info(self):
        print '-'*5, 'Variable:', self.id
        print '-'*10, 'Name:', self.name
        print '-'*10, 'Units:', self.units
        print '-'*10, 'Axes:', self.axlist
        print '-'*10, 'Coordinates:', self.coord
        print '-'*10, 'mean: {0}; min: {1}; max: {2}'.format(np.average(self.data),np.amin(self.data),np.amax(self.data))
        

    def set_coordinates(self,*coord):
        self.axlist = tuple(coord)
        self.coord = " ".join(coord)
 
    def write(self,filein):

        for ax in self.axes:
            ax.write(filein)

        if not(self.data is None):
          tmp = filein.createVariable(self.id,"f8",self.axlist)
          tmp[:] = self.data
          tmp.standard_name = self.name
          tmp.units = self.units
          tmp.coordinates = self.coord

    def plot(self,rep_images=None,var2=None,label="",label2="",timeunits=None,levunits=None):

        coef = self.plotcoef

        if not(self.level is None):
            if levunits is None:
                levs = self.level.data
                levunits = self.level.units
                zlabel = 'Altitude [{0}]'.format(levunits)
            elif levunits == 'hPa' and self.level.units == 'Pa':
                levs = self.level.data/100.
                zlabel = 'Pressure [{0}]'.format(levunits)
            elif (levunits == 'hPa' and self.level.units == 'hPa') or (levunits == 'Pa' and self.level.units == 'Pa'):
                levs = self.level.data
                zlabel = 'Pressure [{0}]'.format(levunits)
            elif levunits == 'km' and self.level.units == 'm':
                levs = self.level.data/1000.
                zlabel = 'Altitude [{0}]'.format(levunits)
            elif (levunits == 'km' and self.level.units == 'km') or (levunits == 'm' and self.level.units == 'm'):
                levs = self.level.data
                zlabel = 'Altitude [{0}]'.format(levunits)
            else:
                print "ERROR: unexpected case for levunits:", levunits, self.level.units
            if not(var2 is None):
                if levunits is None:
                    levs2 = var2.level.data
                elif levunits == 'hPa' and var2.level.units == 'Pa':
                    levs2 = var2.level.data/100.
                elif (levunits == 'hPa' and var2.level.units == 'hPa') or (levunits == 'Pa' and var2.level.units == 'Pa'):
                    levs2 = var2.level.data
                elif levunits == 'km' and var2.level.units == 'm':
                    levs2 = var2.level.data/1000.
                elif (levunits == 'km' and var2.level.units == 'km') or (levunits == 'm' and var2.level.units == 'm'):
                    levs2 = var2.level.data
                else:
                    print "ERROR: unexpected case for levunits (var2):", levunits, var2.level.units


        if not(self.time is None) and not(self.level is None):

            if self.time.length == 1:
                if var2 is None:
                    plotbasics.plot(self.data[0,:]*coef,levs,
                           xlabel='{0} [{1}]'.format(self.id,self.plotunits),
                           ylabel=zlabel,
                           title='{0} ({1})'.format(self.name,self.time.name),
                           rep_images=rep_images,name='{0}.png'.format(self.id),
                           yunits=levunits)
                else:
                    plotbasics.plot(self.data[0,:]*coef,levs,
                            x2=var2.data[0,:]*coef,
                            y2=levs2,xlabel='{0} [{1}]'.format(self.id,self.plotunits),
                            ylabel=zlabel,
                            title='{0} ({1})'.format(self.name,self.time.name),
                            rep_images=rep_images,name='{0}.png'.format(self.id),
                            label=label,label2=label2,
                            yunits=levunits)
            else:
                if timeunits is None:
                    time = self.time.data
                    tunits = self.time.units
                else:
                    if timeunits == 'hours':
                        time = self.time.data/3600.
                        tunits = self.time.units.replace("seconds","hours")
                    elif timeunits == 'days':
                        time = self.time.data/86400.
                        tunits = self.time.units.replace("seconds","days")
                    else:
                        print "ERROR: timeunits unexpected for plotting:", timeunits
                        sys.exit()

                plotbasics.plot2D(time,levs,self.data[:,:]*coef,
                        xlabel=tunits,
                        ylabel=zlabel,
                        title='{0} [{1}]'.format(self.id,self.plotunits),
                        rep_images=rep_images,name='{0}.png'.format(self.id),
                        yunits=levunits)

        elif not(self.time is None):

            if self.time.length > 1:

                if timeunits is None:
                    time = self.time.data
                    if not(var2 is None):
                        time2 = var2.time.data
                    tunits = self.time.units
                else:
                    if timeunits == 'hours':
                        time = self.time.data/3600.
                        if not(var2 is None):
                            time2 = var2.time.data/3600.
                        tunits = self.time.units.replace("seconds","hours")
                    elif timeunits == 'days':
                        time = self.time.data/86400.
                        if not(var2 is None):
                            time2 = var2.time.data/86400.
                        tunits = self.time.units.replace("seconds","days")
                    else:
                        print "ERROR: timeunits unexpected for plotting:", timeunits
                        sys.exit()

                if var2 is None:
                    plotbasics.plot(time,self.data[:]*coef,
                            ylabel='{0} [{1}]'.format(self.id,self.plotunits),
                            xlabel=tunits,
                            title=self.name,
                            rep_images=rep_images,name='{0}.png'.format(self.id)) 
                else:
                    plotbasics.plot(time,self.data[:]*coef,
                            x2=time2,
                            y2=var2.data[:]*coef,
                            ylabel='{0} [{1}]'.format(self.id,self.plotunits),
                            xlabel=tunits,
                            title=self.name,
                            rep_images=rep_images,name='{0}.png'.format(self.id),
                            label=label,label2=label2)
            else:
                print 'no plot for variable', self.id

        elif not(self.level is None):

            if var2 is None:
                plotbasics.plot(self.data[:]*coef,levs,
                        xlabel='{0} [{1}]'.format(self.id,self.plotunits),
                        ylabel=zlabel,
                        title=self.name,
                        rep_images=rep_images,name='{0}.png'.format(self.id),
                        yunits=levunits)
            else:
                plotbasics.plot(self.data[:]*coef,levs,
                        x2=var2.data[:]*coef,
                        y2=levs2,
                        xlabel='{0} [{1}]'.format(self.id,self.plotunits),
                        ylabel=zlabel,
                        title='{0} ({1})'.format(self.name,self.time.name),
                        rep_images=rep_images,name='{0}.png'.format(self.id),
                        label=label,label2=label2,
                        yunits=levunits)

        else:
            print 'no plot for variable', self.id



def read(name,filein):

    tmp = filein[name]

    axlist = []
    axes = []
    for ax in tmp.dimensions:
        try:
            axes.append(Axis(ax,filein[ax][:],name=filein[ax].long_name,units=filein[ax].units))
        except AttributeError:
            axes.append(Axis(ax,filein[ax][:],units=filein[ax].units))
        except:
            raise

        axlist.append(ax)

    try:
        varout = Variable(name,data=tmp[:],name=tmp.long_name,units=tmp.units,axes=axes,axlist=axlist)
    except AttributeError:
        varout = Variable(name,data=tmp[:],units=tmp.units,axes=axes,axlist=axlist)
    except: 
        raise

    return varout


def interpol(var,levout=None,timeout=None,log=False):

    # TODO: Extrapolation to be revisited

    if not(var.time is None) and not(var.level is None):

        ntin, nlevin = var.time.length, var.level.length    
        if not(levout is None):
            nlevout = levout.length
            tmp = np.zeros((ntin,nlevout,1,1),dtype=np.float64)
            for it in range(0,ntin):
                #ff = interpolate.interp1d(var.level.data,var.data[it,:,0,0],bounds_error=False,fill_value="extrapolate")
                if log:
                    ff = interpolate.interp1d(var.level.data,np.log(var.data[it,:,0,0]),bounds_error=False,fill_value=np.log(var.data[it,-1,0,0]))
                    tmp[it,:,0,0] = np.exp(ff(levout.data))
                else:
                    ff = interpolate.interp1d(var.level.data,var.data[it,:,0,0],bounds_error=False,fill_value=var.data[it,-1,0,0])
                    tmp[it,:,0,0] = ff(levout.data)
        else:
            tmp = var.data
            nlevout = nlevin
            levout = var.data.level

        if not(timeout is None) and ntin > 1:
            ntout = timeout.length
            tmp2 = np.zeros((ntout,nlevout,1,1),dtype=np.float64)
            for ilev in range(0,nlevout):
                #ff = interpolate.interp1d(var.time.data,tmp[:,ilev,0,0],bounds_error=False,fill_value="extrapolate")
                ff = interpolate.interp1d(var.time.data,tmp[:,ilev,0,0],bounds_error=False,fill_value=tmp[-1,ilev,0,0])
                tmp2[:,ilev,0,0] = ff(timeout.data)
        else:
            tmp2 = tmp
            ntout = ntin
            timeout = var.time

        varout = Variable(var.id,name=var.name,units=var.units,data=tmp2,time=timeout,level=levout,lat=var.lat,lon=var.lon)

    elif not(var.time is None):

        ntin = var.time.length
        if not(timeout is None) and ntin > 1:
            ntout = timeout.length
            tmp = np.zeros((ntout,1,1),dtype=np.float64)
            #ff = interpolate.interp1d(var.time.data,var.data[:,0,0],bounds_error=False,fill_value="extrapolate")
            ff = interpolate.interp1d(var.time.data,var.data[:,0,0],bounds_error=False,fill_value=var.data[-1,0,0])
            tmp[:,0,0] = ff(timeout.data)
        else:
            tmp = np.log(var.data)
            tmp = var.data
            ntout = ntin
            timeout = var.time

        varout = Variable(var.id,name=var.name,units=var.units,data=tmp,time=timeout,lat=var.lat,lon=var.lon)

    elif not(var.level is None):

        nlevin = var.level.length    
        if not(levout is None):
            nlevout = levout.length
            tmp = np.zeros((nlevout,1,1),dtype=np.float64)
            #ff = interpolate.interp1d(var.level.data,var.data[:,0,0],bounds_error=False,fill_value="extrapolate")
            if log:
                ff = interpolate.interp1d(var.level.data,np.log(var.data[:,0,0]),bounds_error=False,fill_value=var.data[-1,0,0])
                tmp[:,0,0] = np.exp(ff(levout.data))
            else:
                ff = interpolate.interp1d(var.level.data,var.data[:,0,0],bounds_error=False,fill_value=var.data[-1,0,0])
                tmp[:,0,0] = ff(levout.data)                
        else:
            tmp = var.data
            nlevout = nlevin
            levout = data.level

        varout = Variable(var.id,name=var.name,units=var.units,data=tmp,level=levout,lat=var.lat,lon=var.lon)

    else:
        print 'ERROR: Weird, time and level are None for var=', var
        sys.exit()

    return varout
 

