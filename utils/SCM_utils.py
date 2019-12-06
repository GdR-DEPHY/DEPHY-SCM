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

from variables_attributes import attributes as var_attributes
from attributes import known_attributes, required_attributes

import thermo

import plotbasics

class Axis:

    def __init__(self,axisid,data,name=None,units=None,**kwargs):

        self.id = axisid
        self.data = np.array(data,dtype=np.float64)
        self.length, = self.data.shape
        self.units = units
        self.name = name
        for x in kwargs.keys():
            self.__dict__[x] = kwargs[x]

    def info(self,data=False):
        print '-'*10, 'axis id:', self.id
        print '-'*10, 'name:', self.name
        print '-'*10, 'units:', self.units
        print '-'*10, 'length:', self.length
        for x in self.__dict__.keys():
            if not(x in ['id','name','units','length','data']):
                print '-'*10, '{0}: {1}'.format(x,self.__dict__[x])
        if data:
            print '-'*10, 'data:', self.data


    def write(self,filein):

        if not(self.id in filein.dimensions):
            dim = filein.createDimension(self.id, self.length)
            ax = filein.createVariable(self.id,"f8",(self.id,))
            ax[:] = self.data
            ax.long_name = self.name
            ax.units = self.units
            for x in self.__dict__.keys():
                if not(x in ['id','units','name','data','length']):
                    ax.setncattr(x,self.__dict__[x])


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

        self.coord = " ".join(self.axlist)

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
          tmp = filein.createVariable(self.id,"f4",self.axlist)
          tmp[:] = self.data
          tmp.long_name = self.name
          tmp.units = self.units
          #tmp.coordinates = self.coord

    def plot(self,rep_images=None,var2=None,label="",label2="",timeunits=None):

        coef = self.plotcoef

        if not(self.time is None) and not(self.level is None):

            if self.time.length == 1:
                if var2 is None:
                    plotbasics.plot(self.data[0,:,0,0]*coef,self.level.data,
                            xlabel='{0} [{1}]'.format(self.id,self.plotunits),
                            ylabel='Altitude [{0}]'.format(self.level.units),
                            title='{0} ({1})'.format(self.name,self.time.name),
                            rep_images=rep_images,name='{0}.png'.format(self.id))
                else:
                    plotbasics.plot(self.data[0,:,0,0]*coef,self.level.data,
                            x2=var2.data[0,:,0,0]*coef,
                            y2=var2.level.data,xlabel='{0} [{1}]'.format(self.id,self.plotunits),
                            ylabel='Altitude ({0})'.format(self.level.units),
                            title='{0} ({1})'.format(self.name,self.time.name),
                            rep_images=rep_images,name='{0}.png'.format(self.id),
                            label=label,label2=label2)
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
                        print "timeunits unexpected for plotting:", timeunits
                        sys.exit()

                plotbasics.plot2D(time,self.level.data,self.data[:,:,0,0]*coef,
                        xlabel=tunits,
                        ylabel='Altitude [{0}]'.format(self.level.units),
                        title='{0} [{1}]'.format(self.id,self.plotunits),
                        rep_images=rep_images,name='{0}.png'.format(self.id))

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
                        print "timeunits unexpected for plotting:", timeunits
                        sys.exit()

                if var2 is None:
                    plotbasics.plot(time,self.data[:,0,0]*coef,
                            ylabel='{0} [{1}]'.format(self.id,self.plotunits),
                            xlabel=tunits,
                            title=self.name,
                            rep_images=rep_images,name='{0}.png'.format(self.id)) 
                else:
                    plotbasics.plot(time,self.data[:,0,0]*coef,
                            x2=time2,
                            y2=var2.data[:,0,0]*coef,
                            ylabel='{0} [{1}]'.format(self.id,self.plotunits),
                            xlabel=tunits,
                            title=self.name,
                            rep_images=rep_images,name='{0}.png'.format(self.id),
                            label=label,label2=label2)
            else:
                print 'no plot for variable', self.id

        elif not(self.level is None):

            if var2 is None:
                plotbasics.plot(self.data[:,0,0]*coef,self.level.data,
                        xlabel='{0} [{1}]'.format(self.id,self.plotunits),
                        ylabel='Altitude ({0})'.format(self.level.units),
                        title=self.name,
                        rep_images=rep_images,name='{0}.png'.format(self.id))
            else:
                plotbasics.plot(self.data[:,0,0]*coef,self.level.data,
                        x2=var2.data[:,0,0]*coef,
                        y2=var2.level.data,
                        xlabel='{0} [{1}]'.format(self.id,self.plotunits),
                        ylabel='Altitude [{0}]'.format(self.level.units),
                        title='{0} ({1})'.format(self.name,self.time.name),
                        rep_images=rep_images,name='{0}.png'.format(self.id),
                        label=label,label2=label2)

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


def interpol(var,levout=None,timeout=None):

    # TODO: Extrapolation to be revisited

    if not(var.time is None) and not(var.level is None):

        ntin, nlevin = var.time.length, var.level.length    
        if not(levout is None):
            nlevout = levout.length
            tmp = np.zeros((ntin,nlevout,1,1),dtype=np.float64)
            for it in range(0,ntin):
                #ff = interpolate.interp1d(var.level.data,var.data[it,:,0,0],bounds_error=False,fill_value="extrapolate")
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
            tmp = var.data
            ntout = ntin
            timeout = var.time

        varout = Variable(var.id,name=var.name,units=var.units,data=tmp,time=timeout,lat=var.lat,lon=var.lon)

    elif not(var.level is None):

        nlevin = var.level.length    
        if not(levout is None):
            nlevout = levout.length
            tmp = np.zeros((nlevout,1,1),dtype=np.float64)
            ff = interpolate.interp1d(var.level.data,var.data[:,0,0],bounds_error=False,fill_value="extrapolate")
            ff = interpolate.interp1d(var.level.data,var.data[:,0,0],bounds_error=False,fill_value=var.data[-1,0,0])
            tmp[:,0,0] = ff(levout.data)
        else:
            tmp = var.data
            nlevout = nlevin
            levout = data.level

        varout = Variable(var.id,name=var.name,units=var.units,data=tmp,level=levout,lat=var.lat,lon=var.lon)

    else:
        print 'Weird, time and level are None for var=', var
        sys.exit()

    return varout
 
class Case:

    def __init__(self,caseid,lat=None,lon=None,startDate=None,endDate=None,zorog=0.,z0=None):

        self.id = caseid

        # Latitude (degrees_noth) and Longitude (degrees_east)
        self.lat = lat
        self.lon = lon

        if not(lat is None):
            self.latAxis = Axis('lat',[lat,],name='Latitude',units='degrees_north')
        if not(lon is None):
            self.lonAxis = Axis('lon',[lon,],name='Longitude',units='degrees_east')

        # Start en end dates of the simulation
        self.startDate = startDate
        self.endDate = endDate

        # Altitude above sea level (m)
        self.zorog = zorog

        # Roughness length
        self.z0 = z0

        # Initial time axis
        self.t0 = 0
        if not(startDate is None):
            self.tunits = 'seconds since {0}-{1:0>2}-{2:0>2} {3:0>2}:{4:0>2}:{5:3.1f}'.format(startDate[0:4],startDate[4:6],startDate[6:8],startDate[8:10],startDate[10:12],float(startDate[12:14]))

            self.t0Axis = Axis('t0',[self.t0,],name='Initial time',units=self.tunits,calendar='gregorian')

        # Variables
        self.varlist = []
        self.variables = {}

        # Attributes
        self.attlist = ['case','title','reference','author','version','modifications','script','comment',
                'startDate','endDate',
                'zorog'
                ]
        self.attributes = {
                'case': self.id,
                'title': "",
                'reference': "",
                'author': "",
                'version': "Created on " + time.ctime(time.time()),
                'modifications': "",
                'script': "",
                'comment': "",
                'startDate': self.startDate,
                'endDate': self.endDate,
                'zorog': self.zorog
                }
        for att in set(required_attributes).difference(set(self.attlist)):
            self.attlist.append(att)
            if att in ['surfaceType','surfaceForcing','surfaceForcingWind']:
                self.attributes[att] = ""
            else:
                self.attributes[att] = 0

        if not(self.z0 is None):
            self.attlist.append('z0')
            self.attributes['z0'] = self.z0

    def set_dates(self,startDate,endDate):

        self.startDate = startDate
        self.endDate = endDate

        self.t0 = 0
        self.tunits = 'seconds since {0}-{1:0>2}-{2:0>2} {3:0>2}:{4:0>2}:{5:3.1f}'.format(startDate[0:4],startDate[4:6],startDate[6:8],startDate[8:10],startDate[10:12],float(startDate[12:14]))
        self.t0Axis = Axis('t0',[self.t0,],name='Initial time',units=self.tunits,calendar='gregorian')

    def set_latlon(self,lat,lon):

        self.lat = lat
        self.lon = lon

        self.latAxis = Axis('lat',[lat,],name='Latitude',units='degrees_north')
        self.lonAxis = Axis('lon',[lon,],name='Longitude',units='degrees_east')

    def set_attribute(self,attid,attvalue):

        if not(attid in known_attributes):
            print "Warning, attribute {0} is not known. It might not be written in the case output file.".format(attid)

        if not(attid in self.attlist):
            self.attlist.append(attid)

        self.attributes[attid] = attvalue

    def set_title(self,title):

        self.attributes['title'] = title

    def set_comment(self,comment):

        self.attributes['comment'] = comment

    def set_reference(self,ref):

        self.attributes['reference'] = ref

    def set_author(self,author):

        self.attributes['author'] = author

    def set_modifications(self,modifications):

        self.attributes['modifications'] = modifications

    def set_script(self,script):

        self.attributes['script'] = script

    def add_variable(self,varid,vardata,lev=None,levtype=None,levid=None,time=None,timeid=None,name=None,units=None):

        if varid in self.varlist:
            print 'Variable {0} is already defined'.format(varid)
            sys.exit()

        self.varlist.append(varid)

        ######################
        # Level axis, if needed
        if lev is None:
            levAxis=None
            nlev = None
        else:
            try:
                nlev, = lev.shape
            except AttributeError:
                nlev = len(lev)
            except:
                raise
            if levtype == 'altitude':
                if levid is None:
                    levAxis = Axis('lev_{0}'.format(varid),lev,name='{0} for variable {1}'.format(levtype,varid),units='m')
                else:
                    levAxis = Axis(levid,lev,name='{0}'.format(levtype),units='m')
            elif levtype == 'pressure':
                print 'Pressure level type is not coded yet'
                sys.exit()
                if levid is None:
                    levAxis = Axis('lev_{0}'.format(varid),lev,name='{0} for variable {1}'.format(levtype,varid),units='Pa')
                else:
                    levAxis = Axis(levid,lev,name='{0}'.format(levtype),units='Pa')
            else:
                print 'levtype unexpected:', levtype
                print 'levtype should be defined and in altitude, pressure:'
                sys.exit()

        ######################
        # Time axis
        if time is None:
            timeAxis = None
            nt = None
        else:
            # time is supposed to be given in seconds since beginning
            if timeid is None:
                timeAxis = Axis('time_{0}'.format(varid),time,
                        name='Forcing time for variable {0}'.format(varid),
                        units=self.tunits)
            else:
                timeAxis = Axis(timeid,time,name='Forcing time',units=self.tunits)
            try:
                nt, = time.shape
            except AttributeError:
                nt = len(time)
            except:
                raise

        # if variable is an initial variable, impose t0Axis as time axis
        if varid in ['ps','height','pressure','u','v','temp','theta','thetal','qv','qt','rv','rt','ql','qi','tke']:
            timeAxis = self.t0Axis
            nt = 1

        ######################
        # Variable attributes
        if name is None:
            varname = var_attributes[varid]['name']
        else:
            varname = name

        if units is None:
            varunits = var_attributes[varid]['units']
        else:
            print 'Warning: the framework expects SI units'
            varunits = units

        if var_attributes[varid].has_key('plotcoef'):
            plotcoef = var_attributes[varid]['plotcoef']
            plotunits = var_attributes[varid]['plotunits']
        else:
            plotcoef = 1.
            plotunits = None

        ######################
        # Create variable and add to case variables dictionnary
        if nlev is None:
            if nt is None:
                print 'nt=None and nlev=None unexpected'
                sys.exit()
            else:
                tmp = np.reshape(vardata,(nt,1,1))
        else:
            if nt is None:
                tmp = np.reshape(vardata,(nlev,1,1))
            else:
                tmp = np.reshape(vardata,(nt,nlev,1,1))

        self.variables[varid] = Variable(varid,name=varname,units=varunits,data=tmp,level=levAxis,time=timeAxis,lat=self.latAxis,lon=self.lonAxis,plotcoef=plotcoef,plotunits=plotunits)


    def info(self):

        for att in known_attributes:
            if att in self.attlist:
                print '{0}: {1}'.format(att,self.attributes[att])
     
        print "######################"
        print "# Variable information"
        for var in self.varlist:
            self.variables[var].info()

    def write(self,fileout,verbose=False):

        g = nc.Dataset(fileout,'w',format='NETCDF3_CLASSIC')

        for var in var_attributes.keys():
            if var in self.varlist:
                if verbose:
                    self.variables[var].info()
                self.variables[var].write(g)

        for att in known_attributes:
            if att in self.attlist:
                g.setncattr(att,self.attributes[att])

        g.close()

    def read(self,filein,verbose=False):

        f = nc.Dataset(filein,'r')

        self.set_dates(f.startDate,f.endDate)
        self.set_latlon(f['lat'][0],f['lon'][0])

        for var in f.variables:
            if not(var in f.dimensions) and not(var[0:6] == 'bounds') :
                #print var
                tmp = read(var,f)
                if tmp.level is None:
                    self.add_variable(var,tmp.data,time=tmp.time.data)
                else:
                    self.add_variable(var,tmp.data,time=tmp.time.data,lev=tmp.level.data,levtype='altitude')
                if verbose:
                    data[var].info()

        #print self.attlist
        for att in known_attributes:
            #print att
            try:
                self.set_attribute(att,f.getncattr(att))
            except AttributeError:
                pass
            except:
                raise

        f.close()

    def plot(self,rep_images='./images/',timeunits=None):

        if not(os.path.exists(rep_images)):
            os.makedirs(rep_images)

        for var in self.varlist:
            self.variables[var].plot(rep_images=rep_images,timeunits=timeunits)

    def plot_compare(self,cc,rep_images='./images/',label1=None,label2=None,timeunits=None):

        if not(os.path.exists(rep_images)):
            os.makedirs(rep_images)

        for var in self.varlist:
            if var in cc.varlist: 
                if len(self.variables[var].data.shape) <= 3  or self.variables[var].data.shape[0] == 1:
                    self.variables[var].plot(rep_images=rep_images,
                            var2=cc.variables[var],
                            label=label1,label2=label2,
                            timeunits=timeunits)


    def convert2SCM(self,time=None,lev=None,levtype=None):

        if lev is None:
            print 'Output level should be given'
            sys.exit()

        if levtype == 'altitude':
            levout = Axis('lev',lev,name='Altitude',units='m')
        elif levtype == 'pressure':
            levout = Axis('lev',lev,name='pressure',units='Pa')
            print 'Pressure level type is not coded yet'
            sys.exit()
        else:
            print 'levtype unexpected:', levtype
            sys.exit()

        if time is None:
            print 'Output time should be given'
            sys.exit()

        # time should be given in same units as original
        timeout = Axis('time',time,name='time',units=self.tunits)
        ntout = timeout.length

        ###########################
        # Init new case structure
        ###########################

        caseSCM = Case(self.id,lat=self.lat,lon=self.lon,startDate=self.startDate,endDate=self.endDate,zorog=self.zorog,z0=self.z0)
        for att in self.attlist:
            caseSCM.set_attribute(att,self.attributes[att])

        ###########################
        # Interpolation
        ###########################

        dataout = {}
        for var in self.varlist:
            if not(var in ['ps',]):
                dataout[var] = interpol(self.variables[var],levout=levout,timeout=timeout)
            else:
                dataout[var] = self.variables[var]

        ###########################
        # Initial state
        ###########################
 
        caseSCM.add_variable('ps',dataout['ps'].data)

        nt0,nlev,nlat,nlon = dataout['u'].data.shape

        for var in ['height','pressure','u','v','temp','theta','qv','qt','rv','rt','ql','qi','tke']:
            if var in dataout.keys():
                caseSCM.add_variable(var,dataout[var].data,lev=lev,levtype=levtype,levid='lev')
            elif var == 'height':
                caseSCM.add_variable(var,dataout['u'].level.data,lev=lev,levtype=levtype,levid='lev')
            elif var == 'pressure':
                ps = dataout['ps'].data[0,0,0]
                if 'theta' in dataout.keys():
                    z = dataout['theta'].level.data
                    theta = dataout['theta'].data[0,:,0,0]
                    print 'compute pressure from altitude, potential temperature and surface pressure'
                    pressure = thermo.z2p(theta=theta,z=z,ps=ps)

                elif 'temp' in dataout.keys():
                    z = dataout['temp'].level.data
                    temp = dataout['temp'].data[0,:,0,0]
                    print 'compute pressure from altitude, temperature and surface pressure'
                    pressure = thermo.z2p(temp=temp,z=z,ps=ps)                    
                else:
                    print 'theta or temp should be defined'
                    sys.exit()
                pressure = np.reshape(pressure,(1,nlev,1,1))
                caseSCM.add_variable(var,pressure,lev=lev,levtype=levtype,levid='lev')
            elif var == 'theta':
                pressure = caseSCM.variables['pressure'].data[0,:,0,0]
                temp = dataout['temp'].data[0,:,0,0]
                print 'compute potential temperature from pressure and temperature'
                theta = thermo.t2theta(p=pressure,temp=temp)
                caseSCM.add_variable(var,theta,lev=lev,levtype=levtype,levid='lev')                
            elif var == 'temp':
                pressure = caseSCM.variables['pressure'].data[0,:,0,0]
                theta = dataout['theta'].data[0,:,0,0]
                print 'compute temperature from pressure and potential temperature'
                temp = thermo.theta2t(p=pressure,theta=theta)
                caseSCM.add_variable(var,temp,lev=lev,levtype=levtype,levid='lev')
            elif var == 'qv':
                if 'qt' in dataout.keys():
                    print 'Suppose qv=qt'
                    caseSCM.add_variable(var,dataout['qt'].data,lev=lev,levtype=levtype,levid='lev')
                elif 'rv'in dataout.keys():
                    rv = dataout['rv'].data[0,:,0,0]
                    print 'compute qv from rv'
                    qv = thermo.rt2qt(rv)
                    caseSCM.add_variable(var,qv,lev=lev,levtype=levtype,levid='lev')
                elif 'rt' in dataout.keys():
                    rt = dataout['rt'].data[0,:,0,0]
                    print 'compute qt from rt and suppose qv=qt'
                    qv = thermo.rt2qt(rt)
                    caseSCM.add_variable(var,qv,lev=lev,levtype=levtype,levid='lev')
                else:
                    print 'Either qt, rv or rt should be defined'
                    sys.exit()
            elif var == 'qt':
                if 'qv' in dataout.keys():
                    print 'suppose qt=qv'
                    caseSCM.add_variable(var,dataout['qv'].data,lev=lev,levtype=levtype,levid='lev')
                elif 'rv'in dataout.keys():
                    rv = dataout['rv'].data[0,:,0,0]
                    print 'compute qv from rv and suppose qt=qv'
                    qt = thermo.rt2qt(rv)
                    caseSCM.add_variable(var,qt,lev=lev,levtype=levtype,levid='lev')
                elif 'rt' in dataout.keys():
                    rt = dataout['rt'].data[0,:,0,0]
                    print 'compute qt from rt'
                    qt = thermo.rt2qt(rt)
                    caseSCM.add_variable(var,qt,lev=lev,levtype=levtype,levid='lev')
                else:
                    print 'Either qv, rv or rt should be defined'
                    sys.exit()
            elif var == 'rv':
                if 'qv' in dataout.keys():
                    qv = dataout['qv'].data[0,:,0,0]
                    print 'compute rv from qv'
                    rv = thermo.qt2rt(qv)
                    caseSCM.add_variable(var,rv,lev=lev,levtype=levtype,levid='lev')                    
                elif 'qt'in dataout.keys():
                    qt = dataout['qt'].data[0,:,0,0]
                    print 'compute rt from qt and suppose rv=rt'
                    rt = thermo.qt2rt(qt)
                    caseSCM.add_variable(var,rt,lev=lev,levtype=levtype,levid='lev')                     
                elif 'rt' in dataout.keys():
                    print 'suppose rv=rt'
                    caseSCM.add_variable(var,dataout['rt'].data,lev=lev,levtype=levtype,levid='lev')
                else:
                    print 'Either qv, qt or rt should be defined'
                    sys.exit()
            elif var == 'rt':
                if 'qv' in dataout.keys():
                    qv = dataout['qv'].data[0,:,0,0]
                    print 'compute rv from qv and suppose rt=rv'
                    rv = thermo.qt2rt(qv)
                    caseSCM.add_variable(var,rv,lev=lev,levtype=levtype,levid='lev')                    
                elif 'qt'in dataout.keys():
                    qt = dataout['qt'].data[0,:,0,0]
                    print 'compute rt from qt'
                    rt = thermo.qt2rt(qt)
                    caseSCM.add_variable(var,rt,lev=lev,levtype=levtype,levid='lev')                     
                elif 'rv' in dataout.keys():
                    print 'suppose rt=rv'
                    caseSCM.add_variable(var,dataout['rv'].data,lev=lev,levtype=levtype,levid='lev')
                else:
                    print 'either qv, qt or rv should be defined'
                    sys.exit()                
            elif var in ['ql','qi','tke']:
                caseSCM.add_variable(var,dataout['u'].data*0,lev=lev,levtype=levtype,levid='lev')
            else:
                print 'Case unexpected: variable {0} have to be defined'.format(var)
                sys.exit()

        ###########################
        # Forcing
        ###########################

        #---- Surface pressure forcing
        var = 'ps_forc'
        if var in self.varlist:
            caseSCM.add_variable(var,dataout[var].data,time=time,timeid='time')
        else:
            tmp = np.zeros((ntout,1,1),dtype=np.float32) + dataout['ps'].data[0,0,0]
            caseSCM.add_variable(var,tmp,time=time,timeid='time')

        #---- Forcing height
        var = 'height_forc'
        if var in self.varlist:
            caseSCM.add_variable(var,dataout[var].data,lev=lev,levtype=levtype,levid='lev',time=time,timeid='time')
        else:
            tmp = np.zeros((ntout,nlev,1,1),dtype=np.float32)
            for it in range(0,ntout):
                tmp[it,:,0,0] = lev[:]
            caseSCM.add_variable(var,tmp,lev=lev,levtype=levtype,levid='lev',time=time,timeid='time')

        #---- Forcing presssure
        var = 'pressure_forc'
        if var in self.varlist:
            caseSCM.add_variable(var,dataout[var].data,lev=lev,levtype=levtype,levid='lev',time=time,timeid='time')
        else:
            tmp = np.zeros((ntout,nlev,1,1),dtype=np.float32)
            for it in range(0,ntout):
                tmp[it,:,0,0] = caseSCM.variables['pressure'].data[0,:,0,0]
            caseSCM.add_variable(var,tmp,lev=lev,levtype=levtype,levid='lev',time=time,timeid='time')

        #---- Large-scale temperature advection
        var = 'tadv'
        att = 'tadv'
        if att in self.attlist and self.attributes[att] == 1:
            caseSCM.add_variable(var,dataout[var].data,lev=lev,levtype=levtype,levid='lev',time=time,timeid='time')
            print 'compute thadv from tadv'
            pressure = caseSCM.variables['pressure_forc'].data
            tadv = caseSCM.variables['tadv'].data
            thadv = thermo.t2theta(p=pressure,temp=tadv)
            caseSCM.add_variable('thadv',thadv,lev=lev,levtype=levtype,levid='lev',time=time,timeid='time')
            caseSCM.set_attribute('thadv',1)

        var = 'thadv'
        att = 'thadv'
        if att in self.attlist and self.attributes[att] == 1:
            caseSCM.add_variable(var,dataout[var].data,lev=lev,levtype=levtype,levid='lev',time=time,timeid='time')
            print 'compute tadv from thadv'
            pressure = caseSCM.variables['pressure_forc'].data
            thadv = caseSCM.variables['thadv'].data
            tadv = thermo.theta2t(p=pressure,theta=thadv)
            caseSCM.add_variable('tadv',tadv,lev=lev,levtype=levtype,levid='lev',time=time,timeid='time')
            caseSCM.set_attribute('tadv',1)

        #---- Temperature radiative tendency
        var = 'trad'
        att = 'trad'
        if att in self.attlist and self.attributes[att] == 1:
            caseSCM.add_variable(var,dataout[var].data,lev=lev,levtype=levtype,levid='lev',time=time,timeid='time')
            print 'compute thrad from trad'
            pressure = caseSCM.variables['pressure_forc'].data
            trad = caseSCM.variables['trad'].data
            thrad = thermo.t2theta(p=pressure,temp=trad)
            caseSCM.add_variable('thrad',thrad,lev=lev,levtype=levtype,levid='lev',time=time,timeid='time')
            caseSCM.set_attribute('thrad',1)

        if att in self.attlist and self.attributes[att] == "adv":
            caseSCM.set_attribute('thrad',"adv")

        var = 'thrad'
        att = 'thrad'
        if att in self.attlist and self.attributes[att] == 1:
            caseSCM.add_variable(var,dataout[var].data,lev=lev,levtype=levtype,levid='lev',time=time,timeid='time')
            print 'compute trad from thrad'
            pressure = caseSCM.variables['pressure_forc'].data
            thrad = caseSCM.variables['thrad'].data
            trad = thermo.theta2t(p=pressure,theta=thrad)
            caseSCM.add_variable('trad',trad,lev=lev,levtype=levtype,levid='lev',time=time,timeid='time')
            caseSCM.set_attribute('trad',1)

        if att in self.attlist and self.attributes[att] == "adv":
            caseSCM.set_attribute('trad',"adv")

        #---- Large-scale humidity advection
        var = 'qvadv'
        att = 'qvadv'
        if att in self.attlist and self.attributes[att] == 1:
            caseSCM.add_variable(var,dataout[var].data,lev=lev,levtype=levtype,levid='lev',time=time,timeid='time')
            print 'suppose qtadv=qvadv'
            caseSCM.add_variable('qtadv',dataout[var].data,lev=lev,levtype=levtype,levid='lev',time=time,timeid='time')
            caseSCM.set_attribute('qtadv',1)
            print 'Warning: Should compute rvadv and rtadv from qvadv - Supposed equal yet'
            caseSCM.add_variable('rvadv',dataout[var].data,lev=lev,levtype=levtype,levid='lev',time=time,timeid='time')
            caseSCM.set_attribute('rvadv',1)
            caseSCM.add_variable('rtadv',dataout[var].data,lev=lev,levtype=levtype,levid='lev',time=time,timeid='time')
            caseSCM.set_attribute('rtadv',1)

        var = 'qtadv'
        att = 'qtadv'
        if att in self.attlist and self.attributes[att] == 1:
            caseSCM.add_variable(var,dataout[var].data,lev=lev,levtype=levtype,levid='lev',time=time,timeid='time')
            print 'suppose qvadv=qtadv'
            caseSCM.add_variable('qvadv',dataout[var].data,lev=lev,levtype=levtype,levid='lev',time=time,timeid='time')
            caseSCM.set_attribute('qvadv',1)
            print 'Warning: Should compute rvadv and rtadv from qtadv - Supposed equal yet'
            caseSCM.add_variable('rvadv',dataout[var].data,lev=lev,levtype=levtype,levid='lev',time=time,timeid='time')
            caseSCM.set_attribute('rvadv',1)
            caseSCM.add_variable('rtadv',dataout[var].data,lev=lev,levtype=levtype,levid='lev',time=time,timeid='time')
            caseSCM.set_attribute('rtadv',1)

        var = 'rvadv'
        att = 'rvadv'
        if att in self.attlist and self.attributes[att] == 1:
            caseSCM.add_variable(var,dataout[var].data,lev=lev,levtype=levtype,levid='lev',time=time,timeid='time')
            print 'suppose rtadv=rvadv'
            caseSCM.add_variable('rtadv',dataout[var].data,lev=lev,levtype=levtype,levid='lev',time=time,timeid='time')
            caseSCM.set_attribute('rtadv',1)
            print 'Compute qvadv from rvadv'
            rvadv = caseSCM.variables['rvadv'].data
            rv = caseSCM.variables['rv'].data
            qvadv =  thermo.advrt2advqt(rt=rv,advrt=advrv)
            caseSCM.add_variable('qvadv',qvadv,lev=lev,levtype=levtype,levid='lev',time=time,timeid='time')
            caseSCM.set_attribute('qvadv',1)
            print 'suppose qtadv=qvadv'
            caseSCM.add_variable('qtadv',qvadv,lev=lev,levtype=levtype,levid='lev',time=time,timeid='time')
            caseSCM.set_attribute('qtadv',1)

        var = 'rtadv'
        att = 'rtadv'
        if att in self.attlist and self.attributes[att] == 1:
            caseSCM.add_variable(var,dataout[var].data,lev=lev,levtype=levtype,levid='lev',time=time,timeid='time')
            print 'suppose rvadv=rtadv'
            caseSCM.add_variable('rvadv',dataout[var].data,lev=lev,levtype=levtype,levid='lev',time=time,timeid='time')
            caseSCM.set_attribute('rvadv',1)
            print 'Compute qtadv from rtadv'
            rtadv = caseSCM.variables['rtadv'].data
            rt = caseSCM.variables['rt'].data
            qtadv =  thermo.advrt2advqt(rt=rt,advrt=rtadv)
            caseSCM.add_variable('qtadv',qtadv,lev=lev,levtype=levtype,levid='lev',time=time,timeid='time')
            caseSCM.set_attribute('qtadv',1)
            print 'suppose qvadv=qtadv'
            caseSCM.add_variable('qvadv',qtadv,lev=lev,levtype=levtype,levid='lev',time=time,timeid='time')
            caseSCM.set_attribute('qvadv',1)

        #---- Vertical velocity forcing
        var = 'w'
        att = 'forc_w'
        if att in self.attlist and self.attributes[att] == 1:
            caseSCM.add_variable(var,dataout[var].data,lev=lev,levtype=levtype,levid='lev',time=time,timeid='time')

        var = 'omega'
        att = 'forc_omega'
        if att in self.attlist and self.attributes[att] == 1:
            caseSCM.add_variable(var,dataout[var].data,lev=lev,levtype=levtype,levid='lev',time=time,timeid='time')

        #---- Geostrophic wind forcing
        att = 'forc_geo'
        if att in self.attlist and self.attributes[att] == 1:
            caseSCM.add_variable('ug',dataout['ug'].data,lev=lev,levtype=levtype,levid='lev',time=time,timeid='time')
            caseSCM.add_variable('vg',dataout['vg'].data,lev=lev,levtype=levtype,levid='lev',time=time,timeid='time')

        #---- Wind nudging

        var = 'u_nudging'
        att = 'nudging_u'
        if att in self.attlist and self.attributes[att] > 0:
            caseSCM.add_variable(var,dataout[var].data,lev=lev,levtype=levtype,levid='lev',time=time,timeid='time')

        var = 'v_nudging'
        att = 'nudging_v'
        if att in self.attlist and self.attributes[att] > 0:
            caseSCM.add_variable(var,dataout[var].data,lev=lev,levtype=levtype,levid='lev',time=time,timeid='time')

        #---- Temperature nudging

        var = 'temp_nudging'
        att = 'nudging_t'
        if att in self.attlist and self.attributes[att] > 0:
            caseSCM.add_variable(var,dataout[var].data,lev=lev,levtype=levtype,levid='lev',time=time,timeid='time')
            print 'Should compute theta_nudging from temp_nudging - not done yet'
            sys.exit()

        var = 'theta_nudging'
        att = 'nudging_th'
        if att in self.attlist and self.attributes[att] > 0:
            caseSCM.add_variable(var,dataout[var].data,lev=lev,levtype=levtype,levid='lev',time=time,timeid='time')
            print 'compute temp_nudging from theta_nudging'
            pressure = caseSCM.variables['pressure_forc'].data
            thnudg = caseSCM.variables['theta_nudging'].data
            tnudg = thermo.theta2t(p=pressure,theta=thnudg)
            caseSCM.add_variable('temp_nudging',tnudg,lev=lev,levtype=levtype,levid='lev',time=time,timeid='time')
            caseSCM.set_attribute('nudging_t',self.attributes[att])
            caseSCM.set_attribute('p_nudging_t',self.attributes['p_nudging_th'])

        #---- Humidity nudging

        var = 'qv_nudging'
        att = 'nudging_qv'
        if att in self.attlist and self.attributes[att] > 0:
            caseSCM.add_variable(var,dataout[var].data,lev=lev,levtype=levtype,levid='lev',time=time,timeid='time')
            print 'Should compute qt/rv/rt_nudging from qt_nudging - not done yet'
            sys.exit()

        var = 'qt_nudging'
        att = 'nudging_qt'
        if att in self.attlist and self.attributes[att] > 0:
            caseSCM.add_variable(var,dataout[var].data,lev=lev,levtype=levtype,levid='lev',time=time,timeid='time')
            print 'Should compute qv/rv/rt_nudging from qt_nudging - not done yet'
            sys.exit()

        var = 'rv_nudging'
        att = 'nudging_rv'
        if att in self.attlist and self.attributes[att] > 0:
            caseSCM.add_variable(var,dataout[var].data,lev=lev,levtype=levtype,levid='lev',time=time,timeid='time')
            print 'Should compute qv/qt/rt_nudging from rv_nudging - not done yet'
            sys.exit()

        var = 'rt_nudging'
        att = 'nudging_rt'
        if att in self.attlist and self.attributes[att] > 0:
            caseSCM.add_variable(var,dataout[var].data,lev=lev,levtype=levtype,levid='lev',time=time,timeid='time')
            print 'Should compute qv/qt/rv_nudging from rt_nudging - not done yet'
            sys.exit()

        #---- Surface forcing
        att = 'surfaceForcing'
        if att in self.attlist and self.attributes[att] == 'surfaceFlux':
            caseSCM.add_variable('sfc_sens_flx',dataout['sfc_sens_flx'].data,time=time,timeid='time')
            caseSCM.add_variable('sfc_lat_flx', dataout['sfc_lat_flx'].data, time=time,timeid='time')

        if att in self.attlist and self.attributes[att] == 'ts':
            caseSCM.add_variable('ts',dataout['ts'].data,time=time,timeid='time')

        att = 'surfaceForcingWind'
        if att in self.attlist and self.attributes[att] == 'ustar':
            caseSCM.add_variable('ustar',dataout['ustar'].data,time=time,timeid='time')

        ###########################
        # Final
        ###########################

        return caseSCM

