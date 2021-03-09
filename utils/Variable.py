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

    def __init__(self, varid, data=None, name=None, units=None,
            height=None, height_id=None, height_units=None,
            pressure=None, pressure_id=None, pressure_units=None,
            level=None, time=None, lat=None, lon=None,
            axlist=None, axes=None,
            plotcoef=1., plotunits=None):

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
            #self.lat = lat
            #self.lon = lon

            for ax in ['time','level']:#,'lat','lon']:
                if not(self.__dict__[ax] is None):
                    self.axes.append(self.__dict__[ax])
                    self.axlist.append(self.__dict__[ax].id)
        else:
            self.axes = axes
            self.axlist = axlist

            self.time = None
            self.level = None
            #self.lat = None
            #self.lon = None

            for ax in axes:
                if ax.id == 't0' or ax.id[0:4] == 'time':
                    self.time = ax
                elif ax.id[0:3] == 'lev' or ax.id == 'nlev':
                    self.level = ax
                #elif ax.id == 'lat':
                #    self.lat = ax
                #elif ax.id == 'lon':
                #    self.lon = ax
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

        self.height = None
        self.pressure = None

        print self.id, height is None

        if height is not None: # height is privileged over pressure
            if isinstance(height,Axis) or isinstance(height,Variable):
                self.height = Variable(height.id, data=height.data, units=height.units, name=height.name,
                        level=self.level, time=self.time)#, lat=self.lat, lon=self.lon)
                self.coord = " ".join([self.time.id,height.id,'lat','lon'])
                print self.id, self.coord
            else:
                if height_id is None:
                    height_id = 'zh_{0}'.format(self.id)
                height_name = 'height_for_{0}'.format(self.id)
                if height_units is None:
                    height_units = 'm'
                self.height = Variable(height_id, data=height, units=height_units, name=height_name,
                        level=self.level, time=self.time)#, lat=self.lat, lon=self.lon)
                self.coord = " ".join([self.time.id,height_id,'lat','lon'])

            self.height.set_coordinates(self.coord)

        elif pressure is not None:
            
            if isinstance(pressure,Axis) or isinstance(pressure,Variable):
                self.pressure = Variable(pressure.id, data=pressure.data, units=pressure.units, name=pressure.name,
                        level=self.level, time=self.time)#, lat=self.lat, lon=self.lon)
                self.coord = " ".join([self.time.id,pressure.id,'lat','lon'])
            else:
                if pressure_id is None:
                    pressure_id = 'pa_{0}'.format(self.id)
                pressure_name = 'air_pressure_for_{0}'.format(self.id)
                if pressure_units is None:
                    pressure_units = 'Pa'
                self.pressure = Variable(pressure_id, data=pressure, units=pressure_units, name=pressure_name,
                        level=self.level, time=self.time)#, lat=self.lat, lon=self.lon)
                self.coord = " ".join([self.time.id,pressure_id,'lat','lon'])

            self.pressure.set_coordinates(self.coord)

    def info(self):
        print '-'*5, 'Variable:', self.id
        print '-'*10, 'Name:', self.name
        print '-'*10, 'Units:', self.units
        print '-'*10, 'Axes:', self.axlist
        print '-'*10, 'Coordinates:', self.coord
        print '-'*10, 'mean: {0}; min: {1}; max: {2}'.format(np.average(self.data),np.amin(self.data),np.amax(self.data))
        

    def set_coordinates(self, *coord):
        #self.axlist = tuple(coord)
        self.coord = " ".join(coord)

    def set_level(self,lev=None):

        if lev is None:
            print 'WARNING: level is None. Nothing to do'
        else:
            self.level = lev
            newaxlist = []
            newaxes = []
            for ax in self.axes:
                if ax.id[0:3] == 'lev':
                    newaxlist.append(lev.id)
                    newaxes.append(lev)
                else:
                    newaxlist.append(ax.id)
                    newaxes.append(ax)
            self.axes = newaxes
            self.axlist = newaxlist
 
    def write(self, filein,
            write_time_axes=True, write_level_axes=True,
            write_data=True, write_vertical=True):

        if write_time_axes:
            for ax in self.axes:
                if ax.id[0:4] == 'time' or ax.id == 't0':
                    ax.write(filein)

        if write_level_axes:
            for ax in self.axes:
                if ax.id[0:3] == 'lev':
                    ax.write(filein)

        if write_data:
            if self.data is not None:
                if self.id in filein.variables:
                    print 'WARNING: {0} already if netCDF file. Not overwritten'.format(self.id)
                else:
                    tmp = filein.createVariable(self.id, "f8", self.axlist)
                    tmp[:] = self.data
                    tmp.standard_name = self.name
                    tmp.units = self.units
                    tmp.coordinates = self.coord

        if write_vertical:
            if self.height is not None:
                if self.height.id in filein.variables:
                    print 'WARNING: {0} already if netCDF file. Not overwritten'.format(self.height.id)
                else:
                    self.height.write(filein)

            if self.pressure is not None:
                if self.pressure.id in filein.variables:
                    print 'WARNING: {0} already if netCDF file. Not overwritten'.format(self.pressure.id)
                else:
                    self.pressure.write(filein)

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

    def interpol_time(self,time=None):

        if self.time is None:
            print 'ERROR: time interpolation requested for variable {0} which does not have a time axis'.format(self.id)
            raise ValueError
 
        if time is None:
            print 'WARNING: time is None. Thus no time interpolation'
            return self

        ntout, = time.data.shape
        linit = self.data.shape[0] == 1
        l2D = (len(self.data.shape) == 2) and (self.data.shape[0] > 1)

        height = None
        height_id = None
        height_units = None

        pressure = None
        pressure_id = None
        pressure_units = None

        if linit:

            print 'WARNING: Variable "{0}" is an initial state variable. No need for time interpolation.'.format(self.id)
            return self

        elif l2D: # time,level variable

            ntin, nlevin = self.data.shape
            data = np.zeros((ntout,nlevin),dtype=np.float64)
            for ilev in range(0,nlevin):
                ff = interpolate.interp1d(self.time.data, self.data[:,ilev],
                        bounds_error=False, fill_value=self.data[-1,ilev]) # Pad after end date with the last value, if necessary
                data[:,ilev] = ff(time.data)

        
            if self.height is not None:
                height_id = self.height.id
                height_units = self.height.units
                height = np.zeros((ntout,nlevin),dtype=np.float64)
                for ilev in range(0,nlevin):
                    ff = interpolate.interp1d(self.height.time.data, self.height.data[:,ilev],
                            bounds_error=False, fill_value=self.height.data[-1,ilev]) # Pad after end date with the last value, if necessary
                    height[:,ilev] = ff(time.data)

            if self.pressure is not None:
                pressure_id = self.pressure.id
                pressure_units = self.pressure.units
                pressure = np.zeros((ntout,nlevin),dtype=np.float64)
                for ilev in range(0,nlevin):
                    ff = interpolate.interp1d(self.pressure.time.data, self.pressure.data[:,ilev],
                            bounds_error=False, fill_value=self.pressure.data[-1,ilev]) # Pad after end date with the last value, if necessary
                    pressure[:,ilev] = ff(time.data)

        else: # time only variable
 
            ff = interpolate.interp1d(self.time.data, self.data[:],
                        bounds_error=False, fill_value=self.data[-1]) # Pad after end date with the last value, if necessary
            data = ff(time.data)

        return Variable(self.id, data=data, name=self.name, units=self.units,
                level=self.level, time=time,
                height=height, height_id=height_id, height_units=height_units,
                pressure=pressure, pressure_id=pressure_id, pressure_units=pressure_units)

    def interpol_vert(self,height=None,pressure=None,log=False):

        if self.level is None:
            print 'WARNING: vertical interpolation requested for variable {0}, which does not have a level axis'.format(self.id)
            print 'WARNING: simply return original variable'
            return self

        if height is None and pressure is None:
            print "WARNING: height and pressure are None. Thus no vertical interpolaton"
            return self

        ntin, nlevin = self.data.shape

        _height = None
        _height_id = None
        _height_units = None

        _pressure = None
        _pressure_id = None
        _pressure_units = None

        if height is not None and self.height is not None:

            if len(height.shape) == 1:
                _height = np.tile(height,(ntin,1))
            else:
                _height = height

            _, nlevout = _height.shape
            _height_id = self.height.id
            _height_units = self.height.units
            _level = Axis('lev_{0}'.format(self.id), _height[0,:],
                    name='height_for_{0}'.format(self.id), units='m')

            data = np.zeros((ntin,nlevout),dtype=np.float64)
            for it in range(0,ntin):
                ff = interpolate.interp1d(self.height.data[it,:], self.data[it,:],
                        bounds_error=False, fill_value=self.data[it,-1]) # Pad after end date with the last value, if necessary
                data[it,:] = ff(_height[it,:])

        elif pressure is not None and self.pressure is not None:

            if len(pressure.shape) == 1:
                _pressure = np.tile(pressure,(ntin,1))
            else:
                _pressure = pressure

            _, nlevout = _pressure.shape
            _pressure_id = self.pressure.id
            _pressure_units = self.pressure.units
            _level = Axis('lev_{0}'.format(self.id), _pressure[0,:], 
                    name='air_pressure_for_{0}'.format(self.id), units='Pa')

            data = np.zeros((ntin,nlevout),dtype=np.float64)
            for it in range(0,ntin):
                ff = interpolate.interp1d(self.pressure.data[it,:], self.data[it,:],
                        bounds_error=False, fill_value=self.data[it,-1]) # Pad after end date with the last value, if necessary
                data[it,:] = ff(_pressure[it,:])

        else:

            print 'ERROR: case unexpected for vertical interpolation of variable', self.id
            raise ValueError

        return Variable(self.id, data=data, name=self.name, units=self.units,
                level=_level, time=self.time,
                height=_height, height_id=_height_id, height_units=_height_units,
                pressure=_pressure, pressure_id=_pressure_id, pressure_units=_pressure_units)

def read(name,filein):

    tmp = filein[name]

    axlist = []
    axes = []
    for ax in tmp.dimensions:
        kwargs = {}
        for att in filein[ax].ncattrs():
            if att == 'standard_name':
                kwargs['name'] = filein[ax].getncattr(att)
            else:
                kwargs[att] = filein[ax].getncattr(att)

        tmp_ax = Axis(ax,filein[ax][:],**kwargs)
        axes.append(tmp_ax)

        axlist.append(ax)

    coord = tmp.coordinates
    time_id = coord.split()[0]
    var_vert_id = coord.split()[1]

    height = None
    height_id = None
    height_units = None

    pressure = None
    pressure_id = None
    pressure_units = None

    if var_vert_id[0:2] == 'zh':
        height = filein[var_vert_id][:]
        height_id = var_vert_id
        height_units = filein[var_vert_id].units
    elif var_vert_id[0:2] == 'pa':
        pressure = filein[var_vert_id][:]
        pressure_id = var_vert_id
        pressure_units = filein[var_vert_id].units
    #else:
    #    print 'ERROR: vertical variable for {0} is unexpected'.format(name)
    #    raise ValueError

    try:
        varout = Variable(name, data=tmp[:], name=tmp.standard_name, units=tmp.units,
                height=height, height_id=height_id, height_units=height_units,
                pressure=pressure, pressure_id=pressure_id, pressure_units=pressure_units,
                axes=axes, axlist=axlist)
    except AttributeError:
        varout = Variable(name, data=tmp[:], units=tmp.units,
                height=height, height_id=height_id, height_units=height_units,
                pressure=pressure, pressure_id=pressure_id, pressure_units=pressure_units,
                axes=axes, axlist=axlist)
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
 







