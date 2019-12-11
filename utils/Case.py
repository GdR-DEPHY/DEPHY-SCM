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

from Axis import Axis
from Variable import Variable, read as readvar, interpol

from variables_attributes import attributes as var_attributes
from attributes import known_attributes, required_attributes

import thermo

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
            print 'ERROR: Variable {0} is already defined'.format(varid)
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
                #print 'Pressure level type is not coded yet'
                #sys.exit()
                if levid is None:
                    levAxis = Axis('lev_{0}'.format(varid),lev,name='{0} for variable {1}'.format(levtype,varid),units='Pa')
                else:
                    levAxis = Axis(levid,lev,name='{0}'.format(levtype),units='Pa')
            else:
                print 'ERROR: levtype unexpected:', levtype
                print 'ERROR: levtype should be defined and in altitude, pressure:'
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
                print 'ERROR: nt=None and nlev=None unexpected'
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
                tmp = readvar(var,f)
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
            print 'ERROR: Output level should be given'
            sys.exit()

        if levtype == 'altitude':
            levout = Axis('lev',lev,name='Altitude',units='m')
        elif levtype == 'pressure':
            levout = Axis('lev',lev,name='pressure',units='Pa')
            print 'ERROR: Pressure level type is not coded yet'
            sys.exit()
        else:
            print 'ERROR: levtype unexpected:', levtype
            sys.exit()

        if time is None:
            print 'ERROR: Output time should be given'
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

        for var in ['height','pressure','u','v','temp','theta','thetal','qv','qt','rv','rt','ql','qi','tke']:
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
                elif 'thetal' in dataout.keys():
                    z = dataout['thetal'].level.data
                    thetal = dataout['thetal'].data[0,:,0,0]
                    print 'compute pressure from altitude, liquid potential temperature and surface pressure, assuming theta=thetal'
                    pressure = thermo.z2p(theta=thetal,z=z,ps=ps)
                elif 'temp' in dataout.keys():
                    z = dataout['temp'].level.data
                    temp = dataout['temp'].data[0,:,0,0]
                    print 'compute pressure from altitude, temperature and surface pressure'
                    pressure = thermo.z2p(temp=temp,z=z,ps=ps)                    
                else:
                    print 'ERROR: theta, thetal or temp should be defined'
                    sys.exit()
                pressure = np.reshape(pressure,(1,nlev,1,1))
                caseSCM.add_variable(var,pressure,lev=lev,levtype=levtype,levid='lev')
            elif var == 'theta':
                pressure = caseSCM.variables['pressure'].data[0,:,0,0]
                if 'temp' in dataout.keys():
                    temp = dataout['temp'].data[0,:,0,0]
                    print 'compute potential temperature from pressure and temperature'
                    theta = thermo.t2theta(p=pressure,temp=temp)
                elif 'thetal' in dataout.keys():
                    print 'assume theta=thetal'
                    theta = dataout['thetal'].data[0,:,0,0]
                else:
                    print 'ERROR: At least temp or thetal should be given'
                    sys.exit()                    
                caseSCM.add_variable(var,theta,lev=lev,levtype=levtype,levid='lev')
            elif var == 'temp':
                pressure = caseSCM.variables['pressure'].data[0,:,0,0]
                if 'theta' in dataout.keys():
                    theta = dataout['theta'].data[0,:,0,0]
                    print 'compute temperature from pressure and potential temperature'
                    temp = thermo.theta2t(p=pressure,theta=theta)
                elif 'thetal' in dataout.keys():
                    thetal = dataout['thetal'].data[0,:,0,0]
                    print 'compute temperature from pressure and liquid potential temperature (No liquid water considered)'
                    temp = thermo.theta2t(p=pressure,theta=thetal)
                else:
                    print 'ERROR: At least theta or thetal should be given'
                    sys.exit()
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
                    print 'ERROR: Either qt, rv or rt should be defined'
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
                    print 'ERROR: Either qv, rv or rt should be defined'
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
                    print 'ERROR: Either qv, qt or rt should be defined'
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
                    print 'ERROR: Either qv, qt or rv should be defined'
                    sys.exit()                
            elif var in ['ql','qi','tke']:
                caseSCM.add_variable(var,dataout['u'].data*0,lev=lev,levtype=levtype,levid='lev')
            else:
                print 'ERROR: Case unexpected: variable {0} have to be defined'.format(var)
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
            print 'Compute rvadv from qvadv'
            qvadv = caseSCM.variables['qvadv'].data
            qv = caseSCM.variables['qv'].data
            rvadv =  thermo.advqt2advrt(qt=qv,advqt=qvadv)
            caseSCM.add_variable('rvadv',rvadv,lev=lev,levtype=levtype,levid='lev',time=time,timeid='time')
            caseSCM.set_attribute('rvadv',1)
            print 'Suppose rtadv=rvadv'
            caseSCM.add_variable('rtadv',rvadv,lev=lev,levtype=levtype,levid='lev',time=time,timeid='time')
            caseSCM.set_attribute('rtadv',1)

        var = 'qtadv'
        att = 'qtadv'
        if att in self.attlist and self.attributes[att] == 1:
            caseSCM.add_variable(var,dataout[var].data,lev=lev,levtype=levtype,levid='lev',time=time,timeid='time')
            print 'suppose qvadv=qtadv'
            caseSCM.add_variable('qvadv',dataout[var].data,lev=lev,levtype=levtype,levid='lev',time=time,timeid='time')
            caseSCM.set_attribute('qvadv',1)
            print 'Compute rvadv from qvadv'
            qtadv = caseSCM.variables['qtadv'].data
            qt = caseSCM.variables['qt'].data
            rtadv =  thermo.advqt2advrt(qt=qt,advqt=qtadv)
            caseSCM.add_variable('rtadv',rtadv,lev=lev,levtype=levtype,levid='lev',time=time,timeid='time')
            caseSCM.set_attribute('rtadv',1)
            print 'Suppose rvadv=rtadv'
            caseSCM.add_variable('rvadv',rtadv,lev=lev,levtype=levtype,levid='lev',time=time,timeid='time')
            caseSCM.set_attribute('rvadv',1)

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
            qvadv =  thermo.advrt2advqt(rt=rv,advrt=rvadv)
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

        lflag = False

        var = 'temp_nudging'
        att = 'nudging_t'
        if att in self.attlist and self.attributes[att] > 0:
            lflag=True      
            print '-'*10
            print 'temp_nudging is given'
            caseSCM.add_variable(var,dataout[var].data,lev=lev,levtype=levtype,levid='lev',time=time,timeid='time')
            print 'compute theta_nudging from temp_nudging'
            pressure = caseSCM.variables['pressure_forc'].data
            tnudg = caseSCM.variables['temp_nudging'].data
            thnudg = thermo.t2theta(p=pressure,temp=tnudg)
            caseSCM.add_variable('theta_nudging',thnudg,lev=lev,levtype=levtype,levid='lev',time=time,timeid='time')
            caseSCM.set_attribute('nudging_th',self.attributes[att])
            if 'p_nudging_t' in self.attributes.keys():
                caseSCM.set_attribute('p_nudging_th',self.attributes['p_nudging_t'])
            if 'z_nudging_t' in self.attributes.keys():
                caseSCM.set_attribute('z_nudging_th',self.attributes['z_nudging_t'])

        var = 'theta_nudging'
        att = 'nudging_th'
        if att in self.attlist and self.attributes[att] > 0:
            if lflag:
                print 'Error: Several nudging variable for temperature are given, which might yield to inconsistencies'
                sys.exit()
            else:
                lflag=True   
            print '-'*10
            print 'theta_nudging is given'
            caseSCM.add_variable(var,dataout[var].data,lev=lev,levtype=levtype,levid='lev',time=time,timeid='time')
            print 'compute temp_nudging from theta_nudging'
            pressure = caseSCM.variables['pressure_forc'].data
            thnudg = caseSCM.variables['theta_nudging'].data
            tnudg = thermo.theta2t(p=pressure,theta=thnudg)
            caseSCM.add_variable('temp_nudging',tnudg,lev=lev,levtype=levtype,levid='lev',time=time,timeid='time')
            caseSCM.set_attribute('nudging_t',self.attributes[att])
            if 'p_nudging_th' in self.attributes.keys():
                caseSCM.set_attribute('p_nudging_t',self.attributes['p_nudging_th'])
            if 'z_nudging_th' in self.attributes.keys():
                caseSCM.set_attribute('z_nudging_t',self.attributes['z_nudging_th'])

        var = 'thetal_nudging'
        att = 'nudging_thl'
        if att in self.attlist and self.attributes[att] > 0:
            if lflag:
                print 'Error: Several nudging variable for temperature are given, which might yield to inconsistencies'
                sys.exit()
            else:
                lflag=True   
            print '-'*10
            print 'thetal_nudging is given'
            caseSCM.add_variable(var,dataout[var].data,lev=lev,levtype=levtype,levid='lev',time=time,timeid='time')
            print 'suppose theta_nudging=thetal_nudging'
            caseSCM.add_variable('theta_nudging',dataout[var].data,lev=lev,levtype=levtype,levid='lev',time=time,timeid='time')
            caseSCM.set_attribute('nudging_th',self.attributes[att])
            print 'compute temp_nudging from thetal_nudging assuming no liquid water'
            pressure = caseSCM.variables['pressure_forc'].data
            thlnudg = caseSCM.variables['thetal_nudging'].data
            tnudg = thermo.theta2t(p=pressure,theta=thlnudg)
            caseSCM.add_variable('temp_nudging',tnudg,lev=lev,levtype=levtype,levid='lev',time=time,timeid='time')
            caseSCM.set_attribute('nudging_t',self.attributes[att])
            if 'p_nudging_th' in self.attributes.keys():
                caseSCM.set_attribute('p_nudging_th',self.attributes['p_nudging_thl'])
                caseSCM.set_attribute('p_nudging_t',self.attributes['p_nudging_thl'])
            if 'z_nudging_th' in self.attributes.keys():
                caseSCM.set_attribute('z_nudging_th',self.attributes['z_nudging_thl'])
                caseSCM.set_attribute('z_nudging_t',self.attributes['z_nudging_thl'])

        #---- Humidity nudging

        lflag = False

        var = 'qv_nudging'
        att = 'nudging_qv'
        if att in self.attlist and self.attributes[att] > 0:
            lflag = True
            print '-'*10
            print 'qv_nudging is given'            
            caseSCM.add_variable(var,dataout[var].data,lev=lev,levtype=levtype,levid='lev',time=time,timeid='time')
            print 'Suppose qt_nudging=qv_nudging'
            caseSCM.add_variable('qt_nudging',dataout[var].data,lev=lev,levtype=levtype,levid='lev',time=time,timeid='time')
            caseSCM.set_attribute('nudging_qt',self.attributes[att])
            qvnudg = caseSCM.variables['qv_nudging'].data
            print 'compute rv_nudging from qv_nudging and suppose rt_nudging=rv_nudging'
            rvnudg = thermo.qt2rt(qvnudg)
            caseSCM.add_variable('rv_nudging',rvnudg,lev=lev,levtype=levtype,levid='lev',time=time,timeid='time')
            caseSCM.set_attribute('nudging_rv',self.attributes[att])
            caseSCM.add_variable('rt_nudging',rvnudg,lev=lev,levtype=levtype,levid='lev',time=time,timeid='time')
            caseSCM.set_attribute('nudging_rt',self.attributes[att])
            if 'p_nudging_qv' in self.attributes.keys():
                caseSCM.set_attribute('p_nudging_qt',self.attributes['p_nudging_qv'])
                caseSCM.set_attribute('p_nudging_rv',self.attributes['p_nudging_qv'])
                caseSCM.set_attribute('p_nudging_rt',self.attributes['p_nudging_qv'])
            if 'z_nudging_qv' in self.attributes.keys():
                caseSCM.set_attribute('z_nudging_qt',self.attributes['z_nudging_qv'])
                caseSCM.set_attribute('z_nudging_rv',self.attributes['z_nudging_qv'])
                caseSCM.set_attribute('z_nudging_rt',self.attributes['z_nudging_qv'])

        var = 'qt_nudging'
        att = 'nudging_qt'
        if att in self.attlist and self.attributes[att] > 0:
            if lflag:
                print 'Error: Several nudging variable for humidity are given, which might yield to inconsistencies'
                sys.exit()
            else:
                lflag=True
            print '-'*10
            print 'qt_nudging is given'            
            caseSCM.add_variable(var,dataout[var].data,lev=lev,levtype=levtype,levid='lev',time=time,timeid='time')
            print 'Suppose qv_nudging=qt_nudging'
            caseSCM.add_variable('qv_nudging',dataout[var].data,lev=lev,levtype=levtype,levid='lev',time=time,timeid='time')
            caseSCM.set_attribute('nudging_qv',self.attributes[att])
            qtnudg = caseSCM.variables['qt_nudging'].data
            print 'compute rt_nudging from qt_nudging and suppose rv_nudging=rt_nudging'
            rtnudg = thermo.qt2rt(qtnudg)
            caseSCM.add_variable('rt_nudging',rtnudg,lev=lev,levtype=levtype,levid='lev',time=time,timeid='time')
            caseSCM.set_attribute('nudging_rt',self.attributes[att])
            caseSCM.add_variable('rv_nudging',rtnudg,lev=lev,levtype=levtype,levid='lev',time=time,timeid='time')
            caseSCM.set_attribute('nudging_rv',self.attributes[att])
            if 'p_nudging_qt' in self.attributes.keys():
                caseSCM.set_attribute('p_nudging_qv',self.attributes['p_nudging_qt'])
                caseSCM.set_attribute('p_nudging_rv',self.attributes['p_nudging_qt'])
                caseSCM.set_attribute('p_nudging_rt',self.attributes['p_nudging_qt'])
            if 'z_nudging_qt' in self.attributes.keys():
                caseSCM.set_attribute('z_nudging_qv',self.attributes['z_nudging_qt'])
                caseSCM.set_attribute('z_nudging_rv',self.attributes['z_nudging_qt'])
                caseSCM.set_attribute('z_nudging_rt',self.attributes['z_nudging_qt'])

        var = 'rv_nudging'
        att = 'nudging_rv'
        if att in self.attlist and self.attributes[att] > 0:
            if lflag:
                print 'Error: Several nudging variable for humidity are given, which might yield to inconsistencies'
                sys.exit()
            else:
                lflag=True            
            print '-'*10
            print 'rv_nudging is given'
            caseSCM.add_variable(var,dataout[var].data,lev=lev,levtype=levtype,levid='lev',time=time,timeid='time')
            print 'Suppose rt_nudging=rv_nudging'
            caseSCM.add_variable('rt_nudging',dataout[var].data,lev=lev,levtype=levtype,levid='lev',time=time,timeid='time')
            caseSCM.set_attribute('nudging_rt',self.attributes[att])
            rvnudg = caseSCM.variables['rv_nudging'].data
            print 'compute qv_nudging from rv_nudging and suppose qt_nudging=qv_nudging'
            qvnudg = thermo.rt2qt(rvnudg)
            caseSCM.add_variable('qv_nudging',qvnudg,lev=lev,levtype=levtype,levid='lev',time=time,timeid='time')
            caseSCM.set_attribute('nudging_qv',self.attributes[att])
            caseSCM.add_variable('qt_nudging',qvnudg,lev=lev,levtype=levtype,levid='lev',time=time,timeid='time')   
            caseSCM.set_attribute('nudging_qt',self.attributes[att])
            if 'p_nudging_rv' in self.attributes.keys():
                caseSCM.set_attribute('p_nudging_qv',self.attributes['p_nudging_rv'])
                caseSCM.set_attribute('p_nudging_qt',self.attributes['p_nudging_rv'])
                caseSCM.set_attribute('p_nudging_rt',self.attributes['p_nudging_rv'])
            if 'z_nudging_rv' in self.attributes.keys():
                caseSCM.set_attribute('z_nudging_qv',self.attributes['z_nudging_rv'])
                caseSCM.set_attribute('z_nudging_qt',self.attributes['z_nudging_rv'])
                caseSCM.set_attribute('z_nudging_rt',self.attributes['z_nudging_rv'])

        var = 'rt_nudging'
        att = 'nudging_rt'
        if att in self.attlist and self.attributes[att] > 0:
            if lflag:
                print 'Error: Several nudging variable for humidity are given, which might yield to inconsistencies'
                sys.exit()
            else:
                lflag=True            
            print '-'*10
            print 'rt_nudging is given'
            caseSCM.add_variable(var,dataout[var].data,lev=lev,levtype=levtype,levid='lev',time=time,timeid='time')
            print 'Suppose rv_nudging=rt_nudging'
            caseSCM.add_variable('rv_nudging',dataout[var].data,lev=lev,levtype=levtype,levid='lev',time=time,timeid='time')
            caseSCM.set_attribute('nudging_rv',self.attributes[att])
            rtnudg = caseSCM.variables['rt_nudging'].data
            print 'compute qt_nudging from rt_nudging and suppose qv_nudging=qt_nudging'
            qtnudg = thermo.rt2qt(rtnudg)
            caseSCM.add_variable('qt_nudging',qtnudg,lev=lev,levtype=levtype,levid='lev',time=time,timeid='time')
            caseSCM.set_attribute('nudging_qt',self.attributes[att])
            caseSCM.add_variable('qv_nudging',qvnudg,lev=lev,levtype=levtype,levid='lev',time=time,timeid='time') 
            caseSCM.set_attribute('nudging_qv',self.attributes[att])
            if 'p_nudging_rt' in self.attributes.keys():
                caseSCM.set_attribute('p_nudging_qv',self.attributes['p_nudging_rt'])
                caseSCM.set_attribute('p_nudging_qt',self.attributes['p_nudging_rt'])
                caseSCM.set_attribute('p_nudging_rv',self.attributes['p_nudging_rt'])
            if 'z_nudging_rt' in self.attributes.keys():
                caseSCM.set_attribute('z_nudging_qv',self.attributes['z_nudging_rt'])
                caseSCM.set_attribute('z_nudging_qt',self.attributes['z_nudging_rt'])
                caseSCM.set_attribute('z_nudging_rv',self.attributes['z_nudging_rt'])

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

