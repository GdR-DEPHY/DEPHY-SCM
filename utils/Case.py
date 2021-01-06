#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on 27 November 2019

@author: Romain Roehrig
"""

import os
import sys
import time

from datetime import datetime 

import netCDF4 as nc

import numpy as np

from Axis import Axis
from Variable import Variable, read as readvar, interpol

from variables_attributes import attributes as var_attributes
from attributes import known_attributes, required_attributes

import thermo

class Case:

    def __init__(self,caseid,lat=None,lon=None,startDate=None,endDate=None,surfaceType='ocean',zorog=0.,z0=None):

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

        # Surface type
        self.surfaceType = surfaceType

        # Altitude above sea level (m)
        self.zorog = zorog

        # Roughness length
        self.z0 = z0

        # Initial time axis
        self.t0 = 0
        if not(startDate is None):
            self.tunits = 'seconds since {0}-{1:0>2}-{2:0>2} {3:0>2}:{4:0>2}:{5:0>2}'.format(startDate[0:4],startDate[4:6],startDate[6:8],startDate[8:10],startDate[10:12],startDate[12:14])

            self.t0Axis = Axis('t0',[self.t0,],name='Initial time',units=self.tunits,calendar='gregorian')

        # Convert startDate and endDate 

        self.tstart = None
        if startDate is not None:
            d = datetime(int(startDate[0:4]),int(startDate[4:6]),int(startDate[6:8]),int(startDate[8:10]),int(startDate[10:12]),int(startDate[12:14]))
            self.tstart = nc.date2num(d,self.tunits,calendar='gregorian') # should be 0
        self.tend = None
        if endDate is not None:
            d = datetime(int(endDate[0:4]),int(endDate[4:6]),int(endDate[6:8]),int(endDate[8:10]),int(endDate[10:12]),int(endDate[12:14]))
            self.tend = nc.date2num(d,self.tunits,calendar='gregorian')

        # Variables
        self.varlist = []
        self.variables = {}

        # Attributes
        self.attlist = ['case','title','reference','author','version','format_version','modifications','script','comment',
                'startDate','endDate',
                'surfaceType',
                'zorog'
                ]
        self.attributes = {
                'case': self.id,
                'title': "",
                'reference': "",
                'author': "",
                'version': "Created on " + time.ctime(time.time()),
                'format_version': "DEPHY SCM format version 0",
                'modifications': "",
                'script': "",
                'comment': "",
                'startDate': self.startDate,
                'endDate': self.endDate,
                'surfaceType': self.surfaceType,
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
        self.tunits = 'seconds since {0}-{1:0>2}-{2:0>2} {3:0>2}:{4:0>2}:{5:0>2}'.format(startDate[0:4],startDate[4:6],startDate[6:8],startDate[8:10],startDate[10:12],startDate[12:14])
        self.t0Axis = Axis('t0',[self.t0,],name='Initial time',units=self.tunits,calendar='gregorian')

        d = datetime(int(startDate[0:4]),int(startDate[4:6]),int(startDate[6:8]),int(startDate[8:10]),int(startDate[10:12]),int(startDate[12:14]))
        self.tstart = nc.date2num(d,self.tunits,calendar='gregorian') # should be 0
        d = datetime(int(endDate[0:4]),int(endDate[4:6]),int(endDate[6:8]),int(endDate[8:10]),int(endDate[10:12]),int(endDate[12:14]))
        self.tend = nc.date2num(d,self.tunits,calendar='gregorian')

    def set_latlon(self,lat,lon):

        self.lat = lat
        self.lon = lon

        self.latAxis = Axis('lat',[lat,],name='Latitude', units='degrees_north')
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

    def set_z0(self,z0):

        self.attlist.append('z0')
        self.attributes['z0'] = z0

    def add_variable(self,varid,vardata,lev=None,levtype=None,levid=None,time=None,timeid=None,name=None,units=None):
        """ Add a variable to a Case object.
            
            Required arguments:
            varid   -- string for the variable id. Should be in ...
            vardata -- input data as a list or a numpy array

            Optional (keyword) arguments:
            lev     -- input data for the level axis as a list, a numpy array or an Axis object (default None)
            levtype -- string describing the type of the level axis: None (default), 'altitude' or 'pressure'
            levid   -- string for the level axis id. The default is None, which implies a generic id (default None)
            time    -- input data for the time axis, as a list, a numpy array or an Axis object (default None)
            timeid  -- string for the time axis id. The default is None, which implies a generic id (default None)
            name    -- string of the name attribute of the variable (to be use as long_name in a netCDF file) (default None)
            units   -- string of the units attribute of the variable (default None)
        """

        # if variable is already defined, stop
        # if not, add it to case variables dictionnary

        if varid in self.varlist:
            print 'ERROR: Variable {0} is already defined'.format(varid)
            sys.exit()

        self.varlist.append(varid)
        #print varid, lev, time

        ######################
        # Prepare level axis, if needed
        if lev is None:
            levAxis=None
            nlev = None
        elif isinstance(lev,Axis):
            levAxis = lev
            nlev, = lev.data.shape
        else:
            nlev, = np.array(lev).shape # In case lev is given as a list
            if levtype == 'altitude':
                levunits = 'm'
            elif levtype == 'pressure':
                levunits = 'Pa'
            else:
                print 'ERROR: levtype unexpected:', levtype
                print 'ERROR: levtype should be defined and in altitude, pressure:'
                sys.exit()

            if levid is None:
                levAxis = Axis('lev_{0}'.format(varid),lev,name='{0} for variable {1}'.format(levtype,varid),units=levunits)
            else:
                levAxis = Axis(levid,lev,name='{0}'.format(levtype),units=levunits)

        ######################
        # Prepare time axis
        if time is None:
            timeAxis = None
            nt = None
        elif isinstance(time,Axis):
            timeAxis = time
            nt, = time.data.shape
        else:
            nt, = np.array(time).shape # In case time is given as a list
            # time is supposed to be given in seconds since beginning
            if timeid is None:
                timeAxis = Axis('time_{0}'.format(varid),time,
                                name='Forcing time for variable {0}'.format(varid),
                                units=self.tunits)
            else:
                timeAxis = Axis(timeid,time,name='Forcing time',units=self.tunits)

        ######################
        # Get variable attributes
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
        # Create variable
        if levAxis is None and timeAxis is None:
            print 'ERROR: level and time axes are None. Unexpected'
            sys.exit()
        else:
            if timeAxis is None:
                tmp = np.reshape(vardata,(nlev,1,1))
            elif levAxis is None:
                tmp = np.reshape(vardata,(nt,1,1))
            else:
                tmp = np.reshape(vardata,(nt,nlev,1,1))

        self.variables[varid] = Variable(varid,name=varname,units=varunits,data=tmp,level=levAxis,time=timeAxis,lat=self.latAxis,lon=self.lonAxis,plotcoef=plotcoef,plotunits=plotunits)


    def add_init_variable(self,varid,vardata,**kwargs): 
        """ Add an initial state variable to a case object.
            Prepare time axis for such initial state variable and possibly reshape input data to conform with DEPHY format.
            
            Required argument:
            varid   -- id of the initial variable. Should be in ...
            vardata -- input data as a numeric value (int or float), a list or a numpy array

            See add_variable function for optional arguments.
            Note that, for all variable except ps:
            - a level axis is required (lev optional argument).
            - a levtype is required (levtype optional argument).
        """

        # Get time axis for initial state variables
        kwargs['time'] = self.t0Axis

        if varid in ['ps']:
            # Put the expected shape of the input data
            tmp = np.reshape(vardata,(1,1,1))
        else:
            # Check if lev optional argument is given
            if not(kwargs.has_key('lev')):
                print 'ERROR: level axis should be given for variable', varid
                sys.exit()

            nlev, = np.array(kwargs['lev']).shape # In case lev is given as a list

            # Put the expected shape of the input data
            tmp = np.reshape(vardata,(1,nlev,1,1))

        # add initial variable to the Case object
        self.add_variable(varid,tmp,**kwargs)

    def add_init_ps(self,vardata,**kwargs):
        """Add initial state variable for surface pressure to a Case object.
           
           Required argument:
           vardata -- input data as an integer or a float.

           See add_variable function for optional arguments.
        """

        self.add_init_variable('ps',vardata,**kwargs)

    def add_init_height(self,vardata,**kwargs):
        """Add initial state variable for height to a Case object.
           Required argument:
           vardata -- input data as a list or a numpy array.

           See add_variable function for optional arguments.
           Note that:
           - a level axis is required (lev optional argument).
           - a levtype is required (levtype optional argument).
        """

        self.add_init_variable('height',vardata,**kwargs)

    def add_init_pressure(self,vardata,**kwargs):
        """Add initial state variable for pressure to a Case object.
           Required argument:
           vardata -- input data as a list or a numpy array.

           See add_variable function for optional arguments.
           Note that:
           - a level axis is required (lev optional argument).
           - a levtype is required (levtype optional argument).
        """

        self.add_init_variable('pressure',vardata,**kwargs)

    def add_init_temp(self,vardata,**kwargs):
        """Add initial state variable for temperature to a Case object.
           Required argument:
           vardata -- input data as a list or a numpy array.

           See add_variable function for optional arguments.
           Note that:
           - a level axis is required (lev optional argument).
           - a levtype is required (levtype optional argument).
        """

        self.add_init_variable('temp',vardata,**kwargs)

    def add_init_theta(self,vardata,**kwargs):
        """Add initial state variable for potential temperature to a Case object.
           Required argument:
           vardata -- input data as a list or a numpy array.

           See add_variable function for optional arguments.
           Note that:
           - a level axis is required (lev optional argument).
           - a levtype is required (levtype optional argument).
        """

        self.add_init_variable('theta',vardata,**kwargs)

    def add_init_thetal(self,vardata,**kwargs):
        """Add initial state variable for liquid water potential temperature to a Case object.
           Required argument:
           vardata -- input data as a list or a numpy array.

           See add_variable function for optional arguments.
           Note that:
           - a level axis is required (lev optional argument).
           - a levtype is required (levtype optional argument).
        """

        self.add_init_variable('thetal',vardata,**kwargs)

    def add_init_qv(self,vardata,**kwargs):
        """Add initial state variable for specific humidity to a Case object.
           Required argument:
           vardata -- input data as a list or a numpy array.

           See add_variable function for optional arguments.
           Note that:
           - a level axis is required (lev optional argument).
           - a levtype is required (levtype optional argument).
        """

        self.add_init_variable('qv',vardata,**kwargs)

    def add_init_qt(self,vardata,**kwargs):
        """Add initial state variable for total water to a Case object.
           Required argument:
           vardata -- input data as a list or a numpy array.

           See add_variable function for optional arguments.
           Note that:
           - a level axis is required (lev optional argument).
           - a levtype is required (levtype optional argument).
        """

        self.add_init_variable('qt',vardata,**kwargs)

    def add_init_rv(self,vardata,**kwargs):
        """Add initial state variable for water vapor mixing ratio to a Case object.
           Required argument:
           vardata -- input data as a list or a numpy array.

           See add_variable function for optional arguments.
           Note that:
           - a level axis is required (lev optional argument).
           - a levtype is required (levtype optional argument).
        """

        self.add_init_variable('rv',vardata,**kwargs)

    def add_init_rt(self,vardata,**kwargs):
        """Add initial state variable for total water mixing ratio to a Case object.
           Required argument:
           vardata -- input data as a list or a numpy array.

           See add_variable function for optional arguments.
           Note that:
           - a level axis is required (lev optional argument).
           - a levtype is required (levtype optional argument).
        """

        self.add_init_variable('rt',vardata,**kwargs)

    def add_init_wind(self,u=None,v=None,ulev=None,vlev=None,**kwargs):
        """Add initial state variable for total water to a Case object.
           Required arguments:
           u -- input data for zonal wind as a list or a numpy array.
           v -- input data for meridional wind as a list or a numpy array.

           Optional arguments:
           lev -- level axis for both u and v as a list or a numpy array (default None)
           ulev -- level axis for u as a list or a numpy array (default None)
           vlev -- level axis for v as a list or a numpy array (default None)
           levtype -- type of vertical axis (pressure or altitude)

           Either lev or ulev/vlev should be provided

           See add_variable function for optional arguments.
        """

        if u is None or v is None:
            print 'ERROR: you must provide both zonal and meridional wind'
            sys.exit()

        if not(kwargs.has_key('lev')) and ulev is None and vlev is None:
            print 'ERROR: you must provide a vertical axis either with lev or with both ulev/vlev'
            sys.exit()

        if ulev is not None:
            kwargs['lev'] = ulev
        self.add_init_variable('u',u,**kwargs)

        if vlev is not None:
            kwargs['lev'] = vlev
        self.add_init_variable('v',v,**kwargs)

    def add_init_tke(self,vardata,**kwargs):
        """Add initial state variable for turbulent kinetic energy to a Case object.
           Required argument:
           vardata -- input data as a list or a numpy array.

           See add_variable function for optional arguments.
           Note that:
           - a level axis is required (lev optional argument).
           - a levtype is required (levtype optional argument).
        """

        self.add_init_variable('tke',vardata,**kwargs)

    def add_forcing_variable(self,varid,vardata,**kwargs): 
        """ Add a forcing variable to a case object.
            
            Required argument:
            varid   -- id of the forcing variable. Should be in ...
            vardata -- input data as a numeric value (int or float), a list or a numpy array

            See add_variable function for optional arguments.
            Note that, for all variable except ps_forc:
            - a level axis is required (lev optional argument).
            - a levtype is required (levtype optional argument).

            If time is not provided, forcing is assumed constant in time
        """

        # Prepare time axis
        if kwargs.has_key('time'):
            lconstant = False
            nt, = np.array(kwargs['time']).shape
        else: # forcing is constant in time
            lconstant = True
            kwargs['time'] = [self.tstart,self.tend]
            nt = 2

        if varid in ['ps_forc','sfc_sens_flx','sfc_lat_flx','ustar','ts']:
            # Put the expected shape of the input data
            if lconstant: 
                tmp = np.zeros((nt,1,1),dtype=np.float32)
                tmp[0,0,0] = vardata
                tmp[1,0,0] = vardata
            else:
                tmp = np.reshape(vardata,(nt,1,1))
        else:
            # Check if lev optional argument is given
            if not(kwargs.has_key('lev')):
                print 'ERROR: level axis should be given for variable', varid
                sys.exit()

            nlev, = np.array(kwargs['lev']).shape # In case lev is given as a list

            # Put the expected shape of the input data
            if lconstant:
                tmp = np.zeros((nt,nlev,1,1),dtype=np.float32)
                tmp[0,:,0,0] = vardata[:]
                tmp[1,:,0,0] = vardata[:]
            else:
                tmp = np.reshape(vardata,(nt,nlev,1,1))

        # add initial variable to the Case object
        self.add_variable(varid,tmp,**kwargs)

    def add_surface_pressure_forcing(self,data,**kwargs):
        """Add a surface pressure forcing to a Case object
           Required argument:
           data -- input data as a list or a numpy array.

           See add_variable function for optional arguments.

           If time is not provided, forcing is assumed constant in time.           
        """

        self.add_forcing_variable('ps_forc',data,**kwargs)

    def add_pressure_forcing(self,data,**kwargs):
        """Add forcing pressure levels to a Case object
           Required argument:
           data -- input data as a list or a numpy array.

           See add_variable function for optional arguments.
           Note that:
           - a level axis is required (lev optional argument).
           - a levtype is required (levtype optional argument).           

           If time is not provided, forcing is assumed constant in time.           
        """

        self.add_forcing_variable('pressure_forc',data,**kwargs)

    def add_height_forcing(self,data,**kwargs):
        """Add forcing height levels to a Case object
           Required argument:
           data -- input data as a list or a numpy array.

           See add_variable function for optional arguments.
           Note that:
           - a level axis is required (lev optional argument).
           - a levtype is required (levtype optional argument).           

           If time is not provided, forcing is assumed constant in time.           
        """

        self.add_forcing_variable('height_forc',data,**kwargs)

    def add_geostrophic_wind(self,ug=None,vg=None,uglev=None,vglev=None,**kwargs):
        """Add a geostrophic wind forcing to a Case object.
           Required argument:
           ug -- input data for geostrophic zonal wind as a list or a numpy array.
           vg -- input data for geostrophic meridional wind as a list or a numpy array.

           Optional arguments:
           lev     -- level axis for both u and v as a list or a numpy array (default None)
           uglev   -- level axis for u as a list or a numpy array (default None)
           vglev   -- level axis for v as a list or a numpy array (default None)
           levtype -- type of vertical axis (pressure or altitude)

           Either lev or ulev/vlev should be provided

           See add_variable function for optional arguments.

           If time is not provided, forcing is assumed constant in time
        """

        if ug is None or vg is None:
            print 'ERROR: you must provide both zonal and meridional geostrophic wind'
            sys.exit()

        if not(kwargs.has_key('lev')) and uglev is None and vglev is None:
            print 'ERROR: you must provide a vertical axis either with lev or with both uglev/vglev'
            sys.exit()

        self.set_attribute('forc_geo',1)

        if uglev is not None:
            kwargs['lev'] = uglev
        self.add_forcing_variable('ug',ug,**kwargs)

        if vglev is not None:
            kwargs['lev'] = vglev        
        self.add_forcing_variable('vg',vg,**kwargs)

    def add_vertical_velocity(self,w=None,omega=None,**kwargs):
        """Add a potential temperature advection to a Case object.
           Argument:
           w     -- input vertical velocity in m s-1 as a list or a numpy array (default None).
           omega -- input pressure vertical velocity in Pa s-1 as a list or a numpy array (default None).

           Either w or omega should be given

           See add_variable function for optional arguments.
           Note that:
           - a level axis is required (lev optional argument).
           - a levtype is required (levtype optional argument).

           If time is not provided, forcing is assumed constant in time.
        """

        if w is None and omega is None:
            print 'ERROR: you must provide either w or omega'
            sys.exit()

        if w is not None and omega is not None:
            print 'ERROR: you cannot provide both w and omega at the same time'
            sys.exit()

        if w is not None:
            self.set_attribute('forc_w',1)
            self.add_forcing_variable('w',w,**kwargs)

        if omega is not None:
            self.set_attribute('forc_omega',1)
            self.add_forcing_variable('omega',omega,**kwargs)

    def add_temp_advection(self,data,include_rad=False,**kwargs):
        """Add a temperature advection to a Case object.
           Required argument:
           data -- input data as a list or a numpy array.

           Optional (keyword) argument:
           include_rad -- boolean indicated whether the radiative tendency 
                          is included in the advection (default False)

           See add_variable function for optional arguments.
           Note that:
           - a level axis is required (lev optional argument).
           - a levtype is required (levtype optional argument).

           If time is not provided, forcing is assumed constant in time.
        """

        self.set_attribute('adv_temp',1)
        if include_rad:
            self.set_attribute('rad_temp','adv')

        self.add_forcing_variable('temp_adv',data,**kwargs)

    def add_theta_advection(self,data,include_rad=False,**kwargs):
        """Add a potential temperature advection to a Case object.
           Required argument:
           data -- input data as a list or a numpy array.

           Optional (keyword) argument:
           include_rad -- boolean indicated whether the radiative tendency 
                          is included in the advection (default False)

           See add_variable function for optional arguments.
           Note that:
           - a level axis is required (lev optional argument).
           - a levtype is required (levtype optional argument).

           If time is not provided, forcing is assumed constant in time.
        """

        self.set_attribute('adv_theta',1)
        if include_rad:
            self.set_attribute('rad_theta','adv')

        self.add_forcing_variable('theta_adv',data,**kwargs)

    def add_thetal_advection(self,data,include_rad=False,**kwargs):
        """Add a liquid potential temperature advection to a Case object.
           Required argument:
           data -- input data as a list or a numpy array.

           Optional (keyword) argument:
           include_rad -- boolean indicated whether the radiative tendency 
                          is included in the advection (default False)

           See add_variable function for optional arguments.
           Note that:
           - a level axis is required (lev optional argument).
           - a levtype is required (levtype optional argument).

           If time is not provided, forcing is assumed constant in time.
        """

        self.set_attribute('adv_thetal',1)
        if include_rad:
            self.set_attribute('rad_thetal','adv')

        self.add_forcing_variable('thetal_adv',data,**kwargs)

    def add_qv_advection(self,data,**kwargs):
        """Add a specific humidity advection to a Case object.
           Required argument:
           data -- input data as a list or a numpy array.

           See add_variable function for optional arguments.
           Note that:
           - a level axis is required (lev optional argument).
           - a levtype is required (levtype optional argument).

           If time is not provided, forcing is assumed constant in time.
        """

        self.set_attribute('adv_qv',1)

        self.add_forcing_variable('qv_adv',data,**kwargs)

    def add_qt_advection(self,data,**kwargs):
        """Add a total water advection to a Case object.
           Required argument:
           data -- input data as a list or a numpy array.

           See add_variable function for optional arguments.
           Note that:
           - a level axis is required (lev optional argument).
           - a levtype is required (levtype optional argument).

           If time is not provided, forcing is assumed constant in time.
        """

        self.set_attribute('adv_qt',1)

        self.add_forcing_variable('qt_adv',data,**kwargs)

    def add_rv_advection(self,data,**kwargs):
        """Add a water vapor mixing ratio advection to a Case object.
           Required argument:
           data -- input data as a list or a numpy array.

           See add_variable function for optional arguments.
           Note that:
           - a level axis is required (lev optional argument).
           - a levtype is required (levtype optional argument).

           If time is not provided, forcing is assumed constant in time.
        """

        self.set_attribute('adv_rv',1)

        self.add_forcing_variable('rv_adv',data,**kwargs)

    def add_rt_advection(self,data,**kwargs):
        """Add a total water mixing ratio advection to a Case object.
           Required argument:
           data -- input data as a list or a numpy array.

           See add_variable function for optional arguments.
           Note that:
           - a level axis is required (lev optional argument).
           - a levtype is required (levtype optional argument).

           If time is not provided, forcing is assumed constant in time.
        """

        self.set_attribute('adv_rt',1)

        self.add_forcing_variable('rt_adv',data,**kwargs)

    def add_nudging(self,varid,data,timescale=None,z_nudging=None,p_nudging=None,**kwargs):
        """Add a nudging forcing to a Case object.
           Required argument:
           varid     -- id of the variable to be nudged as a string
           data      -- input data as a list or a numpy array.
           timescale -- nudging timescale in seconds (integer or float)
           z_nudging -- altitude above which nudging is applied (integer or float)
           p_nudging -- pressure altitude under which nudging is applied (integer or float)

           Either z_nudging or p_nudging must be defined.

           See add_variable function for optional arguments.
           Note that:
           - a level axis is required (lev optional argument).
           - a levtype is required (levtype optional argument).

           If time is not provided, forcing is assumed constant in time.
        """

        if timescale is None:
            print 'ERROR: you must provide a nudging timescale for variable {0}'.format(varid)
            sys.exit()

        self.set_attribute('nudging_{0}'.format(varid),float(timescale))

        if z_nudging is None and p_nudging is None:
            print 'WARNING: {0} will be nudged over the whole atmosphere'
            self.set_attribute('z_nudging_{0}'.format(varid),0)

        if z_nudging is not None:
            self.set_attribute('z_nudging_{0}'.format(varid),float(z_nudging))
        else:
            self.set_attribute('p_nudging_{0}'.format(varid),float(p_nudging))

        var = '{0}_nudging'.format(varid)
        self.add_forcing_variable(var,data,**kwargs)


    def add_wind_nudging(self,unudg=None,vnudg=None,ulev=None,vlev=None,**kwargs):
        """Add a wind nudging forcing to a Case object.
           Required argument:
           unudg -- input data for zonal wind as a list or a numpy array.
           vnudg -- input data for meridional wind as a list or a numpy array.

           Optional arguments:
           lev     -- level axis for both unudg and vnudg as a list or a numpy array (default None)
           ulev    -- level axis for unudg as a list or a numpy array (default None)
           vlev    -- level axis for vnudg as a list or a numpy array (default None)
           levtype -- type of vertical axis (pressure or altitude)

           Either lev or ulev/vlev should be provided

           See add_variable function for optional arguments.

           If time is not provided, forcing is assumed constant in time
        """

        if unudg is None or vnudg is None:
            print 'ERROR: you must provide both zonal and meridional nudging wind'
            sys.exit()

        if not(kwargs.has_key('lev')) and ulev is None and vlev is None:
            print 'ERROR: you must provide a vertical axis either with lev or with both ulev/vglev'
            sys.exit()

        if ulev is not None:
            kwargs['lev'] = ulev
        self.add_nudging('u',unudg,**kwargs)

        if vlev is not None:
            kwargs['lev'] = vlev        
        self.add_nudging('v',vnudg,**kwargs)

    def add_temp_nudging(self,data,**kwargs):
        """Add a temperature nudging forcing to a Case object.
           Required argument:
           data -- input data as a list or a numpy array.

           See add_nudging for other required arguments

           See add_variable function for optional arguments.
           Note that:
           - a level axis is required (lev optional argument).
           - a levtype is required (levtype optional argument).

           If time is not provided, forcing is assumed constant in time.
        """

        self.add_nudging('temp',data,**kwargs)

    def add_theta_nudging(self,data,**kwargs):
        """Add a potential temperature nudging forcing to a Case object.
           Required argument:
           data -- input data as a list or a numpy array.

           See add_nudging for other required arguments

           See add_variable function for optional arguments.
           Note that:
           - a level axis is required (lev optional argument).
           - a levtype is required (levtype optional argument).

           If time is not provided, forcing is assumed constant in time.
        """

        self.add_nudging('theta',data,**kwargs)

    def add_thetal_nudging(self,data,**kwargs):
        """Add a liquid water potential temperature nudging forcing to a Case object.
           Required argument:
           data -- input data as a list or a numpy array.

           See add_nudging for other required arguments

           See add_variable function for optional arguments.
           Note that:
           - a level axis is required (lev optional argument).
           - a levtype is required (levtype optional argument).

           If time is not provided, forcing is assumed constant in time.
        """

        self.add_nudging('thetal',data,**kwargs)

    def add_qv_nudging(self,data,**kwargs):
        """Add a specific humidity nudging forcing to a Case object.
           Required argument:
           data -- input data as a list or a numpy array.

           See add_nudging for other required arguments

           See add_variable function for optional arguments.
           Note that:
           - a level axis is required (lev optional argument).
           - a levtype is required (levtype optional argument).

           If time is not provided, forcing is assumed constant in time.
        """

        self.add_nudging('qv',data,**kwargs)

    def add_qt_nudging(self,data,**kwargs):
        """Add a total water nudging forcing to a Case object.
           Required argument:
           data -- input data as a list or a numpy array.

           See add_nudging for other required arguments

           See add_variable function for optional arguments.
           Note that:
           - a level axis is required (lev optional argument).
           - a levtype is required (levtype optional argument).

           If time is not provided, forcing is assumed constant in time.
        """

        self.add_nudging('qt',data,**kwargs)

    def add_rv_nudging(self,data,**kwargs):
        """Add a water vapor mixing ratio nudging forcing to a Case object.
           Required argument:
           data -- input data as a list or a numpy array.

           See add_nudging for other required arguments

           See add_variable function for optional arguments.
           Note that:
           - a level axis is required (lev optional argument).
           - a levtype is required (levtype optional argument).

           If time is not provided, forcing is assumed constant in time.
        """

        self.add_nudging('rv',data,**kwargs)

    def add_rt_nudging(self,data,**kwargs):
        """Add a total water mixing ratio nudging forcing to a Case object.
           Required argument:
           data -- input data as a list or a numpy array.

           See add_nudging for other required arguments

           See add_variable function for optional arguments.
           Note that:
           - a level axis is required (lev optional argument).
           - a levtype is required (levtype optional argument).

           If time is not provided, forcing is assumed constant in time.
        """

        self.add_nudging('rt',data,**kwargs)

    def add_ozone(self,data,**kwargs):
        """Add an ozone forcing to a Case object.

           Required argument:
           data -- input data as a list or a numpy array.

           If time is not provided, forcing is assumed constant in time.
        """

        self.add_forcing_variable('o3',data,**kwargs)

    def deactivate_radiation(self):
        """Deactivate radiation in a Case object
           
           No argument required.
        """

        for var in ['temp','theta','thetal']:
            if var in self.varlist:
                self.set_attribute('rad_{0}'.format(var),"adv")

    def add_surface_temp(self,data,**kwargs):
        """Add a surface temperature forcing to a Case object.
           This function does not imply that the surface forcing type is ts.
           This function can be used to add a useful surface temperature in case surface fluxes are prescribed.

           Required argument:
           data -- input data as a list or a numpy array.

           If time is not provided, forcing is assumed constant in time.
        """

        self.add_forcing_variable('ts',data,**kwargs)

    def add_forcing_ts(self,data,z0=None,**kwargs):
        """Add a surface temperature forcing to a Case object.
           This function sets a surface temperature forcing as the case surface forcing.

           Required argument:
           data -- input data as a list or a numpy array.

           If time is not provided, forcing is assumed constant in time.
        """

        self.set_attribute('surfaceForcing','ts')

        if z0 is not None:
            self.set_attribute('surfaceForcingWind','z0')
            self.set_attribute('z0',z0)

        self.add_forcing_variable('ts',data,**kwargs)

    def add_surface_fluxes(self,sens=None,lat=None,time_sens=None,time_lat=None,forc_wind=None,z0=None,ustar=None,time_ustar=None,**kwargs):
        """Add a surface flux forcing to a Case object.

           Required argument:
           sens -- input data for surface sensible heat flux as a numeric, a list or a numpy array.
           lat -- input data for surface latent heat flux as a numeric, a list or a numpy array.

           Optional (keyword) argument:
           time       -- time axis for both sensible and latent heat fluxes
           time_sens  -- time axis for sensibile heat flux
           time_lat   -- time axis for latent heat flux
           forc_wind  -- type of surface wind forcing as a string: 'z0' or 'ustar' (default 'z0')
           z0         -- numeric value for surface roughness (default None)
           ustar      -- input data for surface friction velocity as a numeric, a list or a numpy array (default None)
           time_ustar -- time axis for ustar (default None)

           See add_variable function for other optional arguments.

           If time is not provided, surface fluxes are assumed constant in time.
           If time_ustar is not provided, friction velocity is assumed constant in time.
        """

        if sens is None or lat is None:
            print 'ERROR: you must provide both sensible and latent heat fluxes'
            sys.exit()

        self.set_attribute('surfaceForcing','surfaceFlux')

        if forc_wind == 'z0':
            self.set_attribute("surfaceForcingWind","z0")
            if z0 is None:
                print 'ERROR: z0 must be provided'
                sys.exit()
            self.set_z0(z0)
        elif forc_wind == 'ustar':
            self.set_attribute("surfaceForcingWind","ustar")
            if ustar is None:
                print 'ERROR: ustar must be provided'
                sys.exit()
            self.add_forcing_variable('ustar',ustar,time=time_ustar)
        elif forc_wind is None:
            print 'ERROR: you must specify the surface wind forcing (z0 or ustar)'
        else:
            print 'ERROR: surface wind forcing unexpected:', forc_wind
            print 'ERROR: it should be either z0 or ustar'
            sys.exit()

        if time_sens is not None:
            kwargs['time'] = time_sens
            self.add_forcing_variable('sfc_sens_flx',sens,**kwargs)
        else:
            self.add_forcing_variable('sfc_sens_flx',sens,**kwargs)

        if time_lat is not None:
            kwargs['time'] = time_lat
            self.add_forcing_variable('sfc_lat_flx',lat,**kwargs)
        else:
            self.add_forcing_variable('sfc_lat_flx',lat,**kwargs)

    def set_betaevap(self,beta=1.):
        """Activate a beta model for surface evaporation in a Case object
           
           Optional (keyword) argument:
           beta -- beta value of the beta model (default: 1.)
        """

        self.set_attribute("surfaceForcingMoisture","betaevap")
        self.set_attribute("betaevap",float(beta))

    def deactivate_surface_evaporation(self):
        """Deactivate surface evoporation in a Case object
           
           No argument required.
        """

        self.set_attribute("surfaceForcingMoisture","betaevap")
        self.set_attribute("betaevap",0.)

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
                print 'Reading', var
                tmp = readvar(var,f)
                if tmp.level is None:
                    self.add_variable(var,tmp.data,time=tmp.time.data)
                elif tmp.level.units == 'm':
                    self.add_variable(var,tmp.data,time=tmp.time.data,lev=tmp.level.data,levtype='altitude')
                elif tmp.level.units == 'Pa':
                    self.add_variable(var,tmp.data,time=tmp.time.data,lev=tmp.level.data,levtype='pressure')                    
                if verbose:
                    tmp.info()

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

    def plot(self,rep_images='./images/',timeunits=None,levunits=None):

        if not(os.path.exists(rep_images)):
            os.makedirs(rep_images)

        for var in self.varlist:
            self.variables[var].plot(rep_images=rep_images,timeunits=timeunits,levunits=levunits)

    def plot_compare(self,cc,rep_images='./images/',label1=None,label2=None,timeunits=None,levunits=None):

        if not(os.path.exists(rep_images)):
            os.makedirs(rep_images)

        for var in self.varlist:
            if var in cc.varlist: 
                if len(self.variables[var].data.shape) <= 3  or self.variables[var].data.shape[0] == 1:
                    self.variables[var].plot(rep_images=rep_images,
                            var2=cc.variables[var],
                            label=label1,label2=label2,
                            timeunits=timeunits,levunits=levunits)


    def convert2SCM(self,time=None,lev=None,levtype=None):

        lvert = True
        ltime = True

        if lev is None:
            print 'No vertical interpolation'
            lvert = False
            lev = self.variables['temp'].level.data
            if self.variables['temp'].level.units == 'm':
                levtype = 'altitude'
                levout = Axis('lev',self.variables['temp'].level.data,name='Altitude',units='m')
            elif self.variables['temp'].level.units == 'Pa':
                levtype = 'pressure'
                levout = Axis('lev',self.variables['temp'].level.data,name='pressure',units='Pa')
            else:
                print 'ERROR: Level type undefined for level units:', self.variables['temp'].level.units
                sys.exit()

        else:
            if levtype == 'altitude':
                levout = Axis('lev',lev,name='Altitude',units='m')
            elif levtype == 'pressure':
                levout = Axis('lev',lev,name='pressure',units='Pa')
                print 'ERROR: Pressure level type is not coded yet for interpolation'
                sys.exit()
            else:
                print 'ERROR: levtype unexpected:', levtype
                sys.exit()

        if time is None:
            print 'No time interpolation'
            ltime = False
            try:
                time = self.variables['pressure_forc'].time.data
            except:
                try:
                    time = self.variables['height_forc'].time.data
                except:
                    print 'ERROR: cannot define time axis'
                    raise

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
                if lvert and ltime:
                    dataout[var] = interpol(self.variables[var],levout=levout,timeout=timeout)
                elif lvert or ltime:
                    print 'ERROR: case unexpected for interpolation'
                    print 'ERROR: lvert:', lvert
                    print 'ERROR: ltime:', ltime
                    sys.exit()
                else:
                    dataout[var] = self.variables[var]
            else:
                dataout[var] = self.variables[var]

        ###########################
        # Initial state
        ###########################
 
        caseSCM.add_init_ps(dataout['ps'].data)

        nt0,nlev,nlat,nlon = dataout['u'].data.shape

        for var in ['height','pressure','u','v','temp','theta','thetal','qv','qt','rv','rt','rl','ri','ql','qi','tke']:
            if var in dataout.keys():
                caseSCM.add_init_variable(var,dataout[var].data,lev=lev,levtype=levtype,levid='lev')
            elif var == 'height':
                caseSCM.add_init_variable(var,dataout['u'].level.data,lev=lev,levtype=levtype,levid='lev')
            elif var == 'pressure':
                if var in dataout.keys():
                    caseSCM.add_init_variable(var,dataout[var].data,lev=lev,levtype=levtype,levid='lev')
                else:
                    ps = dataout['ps'].data[0,0,0]
                    if 'theta' in dataout.keys():
                        z = dataout['theta'].level.data
                        theta = dataout['theta'].data[0,:,0,0]
                        if 'qv' in dataout.keys():
                            qv = dataout['qv'].data[0,:,0,0]
                            print 'compute pressure from altitude, potential temperature, qv and surface pressure'
                            pressure = thermo.z2p(theta=theta,z=z,ps=ps,qv=qv)
                        elif 'qt' in dataout.keys():
                            qt = dataout['qt'].data[0,:,0,0]
                            print 'compute pressure from altitude, potential temperature, qt and surface pressure'
                            pressure = thermo.z2p(theta=theta,z=z,ps=ps,qv=qt)
                        else:
                            print 'compute pressure from altitude, potential temperature and surface pressure'
                            pressure = thermo.z2p(theta=theta,z=z,ps=ps)

                    elif 'thetal' in dataout.keys():
                        z = dataout['thetal'].level.data
                        thetal = dataout['thetal'].data[0,:,0,0]
                        if 'qv' in dataout.keys():
                            qv = dataout['qv'].data[0,:,0,0]
                            print 'compute pressure from altitude, liquid potential temperature, qv and surface pressure'
                            pressure = thermo.z2p(theta=thetal,z=z,ps=ps,qv=qv)
                        elif 'qt' in dataout.keys():
                            qt = dataout['qt'].data[0,:,0,0]
                            print 'compute pressure from altitude, liquid potential temperature, qt and surface pressure'
                            pressure = thermo.z2p(theta=thetal,z=z,ps=ps,qv=qt)
                        else:
                            print 'compute pressure from altitude, potential temperature and surface pressure'
                            pressure = thermo.z2p(theta=thetal,z=z,ps=ps)

                    elif 'temp' in dataout.keys():
                        z = dataout['temp'].level.data
                        temp = dataout['temp'].data[0,:,0,0]
                        if 'qv' in dataout.keys():
                            qv = dataout['qv'].data[0,:,0,0]
                            print 'compute pressure from altitude, temperature, qv and surface pressure'
                            pressure = thermo.z2p(temp=temp,z=z,ps=ps,qv=qv)
                        elif 'qt' in dataout.keys():
                            qt = dataout['qt'].data[0,:,0,0]
                            print 'compute pressure from altitude, temperature, qt and surface pressure'
                            pressure = thermo.z2p(temp=temp,z=z,ps=ps,qv=qt)
                        else:
                            print 'compute pressure from altitude, temperature and surface pressure'
                            pressure = thermo.z2p(temp=temp,z=z,ps=ps)

                    else:
                        print 'ERROR: theta, thetal or temp should be defined'
                        sys.exit()
                    pressure = np.reshape(pressure,(1,nlev,1,1))
                    caseSCM.add_init_variable(var,pressure,lev=lev,levtype=levtype,levid='lev')
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
                caseSCM.add_init_variable(var,theta,lev=lev,levtype=levtype,levid='lev')
            elif var == 'thetal':
                pressure = caseSCM.variables['pressure'].data[0,:,0,0]
                if 'temp' in dataout.keys():
                    temp = dataout['temp'].data[0,:,0,0]
                    print 'compute potential temperature from pressure and temperature and assume thetal=theta'
                    thetal = thermo.t2theta(p=pressure,temp=temp)
                elif 'theta' in dataout.keys():
                    print 'assume thetal=theta'
                    thetal = dataout['theta'].data[0,:,0,0]
                else:
                    print 'ERROR: At least temp or thetal should be given'
                    sys.exit()                    
                caseSCM.add_init_variable(var,thetal,lev=lev,levtype=levtype,levid='lev')
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
                caseSCM.add_init_variable(var,temp,lev=lev,levtype=levtype,levid='lev')
            elif var == 'qv':
                if 'qt' in dataout.keys():
                    print 'assume qv=qt'
                    caseSCM.add_init_variable(var,dataout['qt'].data,lev=lev,levtype=levtype,levid='lev')
                elif 'rv'in dataout.keys():
                    rv = dataout['rv'].data[0,:,0,0]
                    print 'compute qv from rv'
                    qv = thermo.rt2qt(rv)
                    caseSCM.add_init_variable(var,qv,lev=lev,levtype=levtype,levid='lev')
                elif 'rt' in dataout.keys():
                    rt = dataout['rt'].data[0,:,0,0]
                    print 'compute qt from rt and assume qv=qt'
                    qv = thermo.rt2qt(rt)
                    caseSCM.add_init_variable(var,qv,lev=lev,levtype=levtype,levid='lev')
                else:
                    print 'ERROR: Either qt, rv or rt should be defined'
                    sys.exit()
            elif var == 'qt':
                if 'qv' in dataout.keys():
                    print 'assume qt=qv'
                    caseSCM.add_init_variable(var,dataout['qv'].data,lev=lev,levtype=levtype,levid='lev')
                elif 'rv'in dataout.keys():
                    rv = dataout['rv'].data[0,:,0,0]
                    print 'compute qv from rv and assume qt=qv'
                    qt = thermo.rt2qt(rv)
                    caseSCM.add_init_variable(var,qt,lev=lev,levtype=levtype,levid='lev')
                elif 'rt' in dataout.keys():
                    rt = dataout['rt'].data[0,:,0,0]
                    print 'compute qt from rt'
                    qt = thermo.rt2qt(rt)
                    caseSCM.add_init_variable(var,qt,lev=lev,levtype=levtype,levid='lev')
                else:
                    print 'ERROR: Either qv, rv or rt should be defined'
                    sys.exit()
            elif var == 'rv':
                if 'qv' in dataout.keys():
                    qv = dataout['qv'].data[0,:,0,0]
                    print 'compute rv from qv'
                    rv = thermo.qt2rt(qv)
                    caseSCM.add_init_variable(var,rv,lev=lev,levtype=levtype,levid='lev')                    
                elif 'qt'in dataout.keys():
                    qt = dataout['qt'].data[0,:,0,0]
                    print 'compute rt from qt and assume rv=rt'
                    rt = thermo.qt2rt(qt)
                    caseSCM.add_init_variable(var,rt,lev=lev,levtype=levtype,levid='lev')                     
                elif 'rt' in dataout.keys():
                    print 'assume rv=rt'
                    caseSCM.add_init_variable(var,dataout['rt'].data,lev=lev,levtype=levtype,levid='lev')
                else:
                    print 'ERROR: Either qv, qt or rt should be defined'
                    sys.exit()
            elif var == 'rt':
                if 'qv' in dataout.keys():
                    qv = dataout['qv'].data[0,:,0,0]
                    print 'compute rv from qv and assume rt=rv'
                    rv = thermo.qt2rt(qv)
                    caseSCM.add_init_variable(var,rv,lev=lev,levtype=levtype,levid='lev')                    
                elif 'qt'in dataout.keys():
                    qt = dataout['qt'].data[0,:,0,0]
                    print 'compute rt from qt'
                    rt = thermo.qt2rt(qt)
                    caseSCM.add_init_variable(var,rt,lev=lev,levtype=levtype,levid='lev')                     
                elif 'rv' in dataout.keys():
                    print 'assume rt=rv'
                    caseSCM.add_init_variable(var,dataout['rv'].data,lev=lev,levtype=levtype,levid='lev')
                else:
                    print 'ERROR: Either qv, qt or rv should be defined'
                    sys.exit()                
            elif var in ['rl','ri','ql','qi','tke']:
                caseSCM.add_init_variable(var,dataout['u'].data*0,lev=lev,levtype=levtype,levid='lev')
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
        var = 'temp_adv'
        att = 'adv_temp'
        if att in self.attlist and self.attributes[att] == 1:
            caseSCM.add_variable(var,dataout[var].data,lev=lev,levtype=levtype,levid='lev',time=time,timeid='time')
            print 'compute theta_adv from temp_adv'
            pressure = caseSCM.variables['pressure_forc'].data
            tadv = caseSCM.variables['temp_adv'].data
            thadv = thermo.t2theta(p=pressure,temp=tadv)
            caseSCM.add_variable('theta_adv',thadv,lev=lev,levtype=levtype,levid='lev',time=time,timeid='time')
            caseSCM.set_attribute('adv_theta',1)
            print 'assume thetal_adv=theta_adv'
            caseSCM.add_variable('thetal_adv',thadv,lev=lev,levtype=levtype,levid='lev',time=time,timeid='time')
            caseSCM.set_attribute('adv_thetal',1)

        var = 'theta_adv'
        att = 'adv_theta'
        if att in self.attlist and self.attributes[att] == 1:
            caseSCM.add_variable(var,dataout[var].data,lev=lev,levtype=levtype,levid='lev',time=time,timeid='time')
            print 'compute temp_adv from theta_adv'
            pressure = caseSCM.variables['pressure_forc'].data
            thadv = caseSCM.variables['theta_adv'].data
            tadv = thermo.theta2t(p=pressure,theta=thadv)
            caseSCM.add_variable('temp_adv',tadv,lev=lev,levtype=levtype,levid='lev',time=time,timeid='time')
            caseSCM.set_attribute('adv_temp',1)
            print 'assumee thetal_adv=theta_adv'
            caseSCM.add_variable('thetal_adv',thadv,lev=lev,levtype=levtype,levid='lev',time=time,timeid='time')
            caseSCM.set_attribute('adv_thetal',1)

        #---- Temperature radiative tendency
        var = 'temp_rad'
        att = 'rad_temp'
        if att in self.attlist and self.attributes[att] == 1:
            caseSCM.add_variable(var,dataout[var].data,lev=lev,levtype=levtype,levid='lev',time=time,timeid='time')
            print 'compute theta_rad from temp_rad'
            pressure = caseSCM.variables['pressure_forc'].data
            trad = caseSCM.variables['temp_rad'].data
            thrad = thermo.t2theta(p=pressure,temp=trad)
            caseSCM.add_variable('theta_rad',thrad,lev=lev,levtype=levtype,levid='lev',time=time,timeid='time')
            caseSCM.set_attribute('rad_theta',1)
            print 'assume thetal_rad=theta_rad'
            caseSCM.add_variable('thetal_rad',thrad,lev=lev,levtype=levtype,levid='lev',time=time,timeid='time')
            caseSCM.set_attribute('rad_thetal',1)

        if att in self.attlist and self.attributes[att] == "adv":
            caseSCM.set_attribute('rad_theta',"adv")
            caseSCM.set_attribute('rad_thetal',"adv")

        var = 'theta_rad'
        att = 'rad_theta'
        if att in self.attlist and self.attributes[att] == 1:
            caseSCM.add_variable(var,dataout[var].data,lev=lev,levtype=levtype,levid='lev',time=time,timeid='time')
            print 'compute temp_rad from theta_rad'
            pressure = caseSCM.variables['pressure_forc'].data
            thrad = caseSCM.variables['theta_rad'].data
            trad = thermo.theta2t(p=pressure,theta=thrad)
            caseSCM.add_variable('temp_rad',trad,lev=lev,levtype=levtype,levid='lev',time=time,timeid='time')
            caseSCM.set_attribute('rad_temp',1)
            print 'assume thetal_rad=theta_rad'
            caseSCM.add_variable('thetal_rad',thrad,lev=lev,levtype=levtype,levid='lev',time=time,timeid='time')
            caseSCM.set_attribute('rad_thetal',1)

        if att in self.attlist and self.attributes[att] == "adv":
            caseSCM.set_attribute('rad_temp',"adv")
            caseSCM.set_attribute('rad_thetal',"adv")

        var = 'thetal_rad'
        att = 'rad_thetal'
        if att in self.attlist and self.attributes[att] == 1:
            caseSCM.add_variable(var,dataout[var].data,lev=lev,levtype=levtype,levid='lev',time=time,timeid='time')
            print 'compute temp_rad from thetal_rad assuming thetal_rad=theta_rad'
            pressure = caseSCM.variables['pressure_forc'].data
            thlrad = caseSCM.variables['thetal_rad'].data
            trad = thermo.theta2t(p=pressure,theta=thlrad)
            caseSCM.add_variable('temp_rad',trad,lev=lev,levtype=levtype,levid='lev',time=time,timeid='time')
            caseSCM.set_attribute('rad_temp',1)
            print 'assume theta_rad=thetal_rad'
            caseSCM.add_variable('theta_rad',thlrad,lev=lev,levtype=levtype,levid='lev',time=time,timeid='time')
            caseSCM.set_attribute('rad_theta',1)

        if att in self.attlist and self.attributes[att] == "adv":
            caseSCM.set_attribute('rad_temp',"adv")
            caseSCM.set_attribute('rad_theta',"adv")

        #---- Large-scale humidity advection
        var = 'qv_adv'
        att = 'adv_qv'
        if att in self.attlist and self.attributes[att] == 1:
            caseSCM.add_variable(var,dataout[var].data,lev=lev,levtype=levtype,levid='lev',time=time,timeid='time')
            print 'assume qt_adv=qv_adv'
            caseSCM.add_variable('qt_adv',dataout[var].data,lev=lev,levtype=levtype,levid='lev',time=time,timeid='time')
            caseSCM.set_attribute('adv_qt',1)
            print 'Compute rv_adv from qv_adv'
            qvadv = caseSCM.variables['qv_adv'].data
            qv = caseSCM.variables['qv'].data
            rvadv =  thermo.advqt2advrt(qt=qv,advqt=qvadv)
            caseSCM.add_variable('rv_adv',rvadv,lev=lev,levtype=levtype,levid='lev',time=time,timeid='time')
            caseSCM.set_attribute('adv_rv',1)
            print 'assume rt_adv=rv_adv'
            caseSCM.add_variable('rt_adv',rvadv,lev=lev,levtype=levtype,levid='lev',time=time,timeid='time')
            caseSCM.set_attribute('adv_rt',1)

        var = 'qt_adv'
        att = 'adv_qt'
        if att in self.attlist and self.attributes[att] == 1:
            caseSCM.add_variable(var,dataout[var].data,lev=lev,levtype=levtype,levid='lev',time=time,timeid='time')
            print 'assume qv_adv=qt_adv'
            caseSCM.add_variable('qv_adv',dataout[var].data,lev=lev,levtype=levtype,levid='lev',time=time,timeid='time')
            caseSCM.set_attribute('adv_qv',1)
            print 'Compute rv_adv from qv_adv'
            qtadv = caseSCM.variables['qt_adv'].data
            qt = caseSCM.variables['qt'].data
            rtadv =  thermo.advqt2advrt(qt=qt,advqt=qtadv)
            caseSCM.add_variable('rt_adv',rtadv,lev=lev,levtype=levtype,levid='lev',time=time,timeid='time')
            caseSCM.set_attribute('adv_rt',1)
            print 'assume rv_adv=rt_adv'
            caseSCM.add_variable('rv_adv',rtadv,lev=lev,levtype=levtype,levid='lev',time=time,timeid='time')
            caseSCM.set_attribute('adv_rv',1)

        var = 'rv_adv'
        att = 'adv_rv'
        if att in self.attlist and self.attributes[att] == 1:
            caseSCM.add_variable(var,dataout[var].data,lev=lev,levtype=levtype,levid='lev',time=time,timeid='time')
            print 'assume rt_adv=rv_adv'
            caseSCM.add_variable('rt_adv',dataout[var].data,lev=lev,levtype=levtype,levid='lev',time=time,timeid='time')
            caseSCM.set_attribute('adv_rt',1)
            print 'Compute qv_adv from rv_adv'
            rvadv = caseSCM.variables['rv_adv'].data
            rv = caseSCM.variables['rv'].data
            qvadv =  thermo.advrt2advqt(rt=rv,advrt=rvadv)
            caseSCM.add_variable('qv_adv',qvadv,lev=lev,levtype=levtype,levid='lev',time=time,timeid='time')
            caseSCM.set_attribute('adv_qv',1)
            print 'assume qt_adv=qv_adv'
            caseSCM.add_variable('qt_adv',qvadv,lev=lev,levtype=levtype,levid='lev',time=time,timeid='time')
            caseSCM.set_attribute('adv_qt',1)

        var = 'rt_adv'
        att = 'adv_rt'
        if att in self.attlist and self.attributes[att] == 1:
            caseSCM.add_variable(var,dataout[var].data,lev=lev,levtype=levtype,levid='lev',time=time,timeid='time')
            print 'assume rv_adv=rt_adv'
            caseSCM.add_variable('rv_adv',dataout[var].data,lev=lev,levtype=levtype,levid='lev',time=time,timeid='time')
            caseSCM.set_attribute('adv_rv',1)
            print 'Compute qt_adv from rt_adv'
            rtadv = caseSCM.variables['rt_adv'].data
            rt = caseSCM.variables['rt'].data
            qtadv =  thermo.advrt2advqt(rt=rt,advrt=rtadv)
            caseSCM.add_variable('qt_adv',qtadv,lev=lev,levtype=levtype,levid='lev',time=time,timeid='time')
            caseSCM.set_attribute('adv_qt',1)
            print 'assume qv_adv=qt_adv'
            caseSCM.add_variable('qv_adv',qtadv,lev=lev,levtype=levtype,levid='lev',time=time,timeid='time')
            caseSCM.set_attribute('adv_qv',1)

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
            if not('z_nudging_u' in self.attributes.keys()):
                height = np.squeeze(caseSCM.variables['height'].data)
                pressure = np.squeeze(caseSCM.variables['pressure'].data)
                zlev = thermo.plev2zlev(self.attributes['p_nudging_u'],height,pressure)
                caseSCM.set_attribute('z_nudging_u',zlev)
            if not('p_nudging_u' in self.attributes.keys()):
                height = np.squeeze(caseSCM.variables['height'].data)
                pressure = np.squeeze(caseSCM.variables['pressure'].data)
                plev = thermo.zlev2plev(self.attributes['z_nudging_u'],height,pressure)
                caseSCM.set_attribute('p_nudging_u',plev)

        var = 'v_nudging'
        att = 'nudging_v'
        if att in self.attlist and self.attributes[att] > 0:
            caseSCM.add_variable(var,dataout[var].data,lev=lev,levtype=levtype,levid='lev',time=time,timeid='time')
            if not('z_nudging_v' in self.attributes.keys()):
                height = np.squeeze(caseSCM.variables['height'].data)
                pressure = np.squeeze(caseSCM.variables['pressure'].data)
                zlev = thermo.plev2zlev(self.attributes['p_nudging_v'],height,pressure)
                caseSCM.set_attribute('z_nudging_v',zlev)
            if not('p_nudging_v' in self.attributes.keys()):
                height = np.squeeze(caseSCM.variables['height'].data)
                pressure = np.squeeze(caseSCM.variables['pressure'].data)
                plev = thermo.zlev2plev(self.attributes['z_nudging_v'],height,pressure)
                caseSCM.set_attribute('p_nudging_v',plev)

        #---- Temperature nudging

        lflag = False

        var = 'temp_nudging'
        att = 'nudging_temp'
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
            print 'assume thetal_nudging=theta_nudging'
            caseSCM.add_variable('thetal_nudging',thnudg,lev=lev,levtype=levtype,levid='lev',time=time,timeid='time')
            caseSCM.set_attribute('nudging_theta',self.attributes[att])
            caseSCM.set_attribute('nudging_thetal',self.attributes[att])
            if 'p_nudging_temp' in self.attributes.keys():
                caseSCM.set_attribute('p_nudging_theta',self.attributes['p_nudging_temp'])
                if not('z_nudging_temp' in self.attributes.keys()):
                    height = np.squeeze(caseSCM.variables['height'].data)
                    pressure = np.squeeze(caseSCM.variables['pressure'].data)
                    zlev = int(thermo.plev2zlev(self.attributes['p_nudging_temp'],height,pressure))
                    caseSCM.set_attribute('z_nudging_thetal',zlev)
                    caseSCM.set_attribute('z_nudging_theta',zlev)
                    caseSCM.set_attribute('z_nudging_temp',zlev)                  
            if 'z_nudging_temp' in self.attributes.keys():
                caseSCM.set_attribute('z_nudging_theta',self.attributes['z_nudging_temp'])
                if not('p_nudging_temp' in self.attributes.keys()):
                    height = np.squeeze(caseSCM.variables['height'].data)
                    pressure = np.squeeze(caseSCM.variables['pressure'].data)
                    plev = int(thermo.zlev2plev(self.attributes['z_nudging_temp'],height,pressure))
                    caseSCM.set_attribute('p_nudging_thetal',plev)
                    caseSCM.set_attribute('p_nudging_theta',plev)
                    caseSCM.set_attribute('p_nudging_temp',plev) 

        var = 'theta_nudging'
        att = 'nudging_theta'
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
            caseSCM.set_attribute('nudging_temp',self.attributes[att])
            if 'p_nudging_theta' in self.attributes.keys():
                caseSCM.set_attribute('p_nudging_temp',self.attributes['p_nudging_theta'])
                if not('z_nudging_theta' in self.attributes.keys()):
                    height = np.squeeze(caseSCM.variables['height'].data)
                    pressure = np.squeeze(caseSCM.variables['pressure'].data)
                    zlev = int(thermo.plev2zlev(self.attributes['p_nudging_theta'],height,pressure))
                    caseSCM.set_attribute('z_nudging_thetal',zlev)
                    caseSCM.set_attribute('z_nudging_theta',zlev)
                    caseSCM.set_attribute('z_nudging_temp',zlev)  
            if 'z_nudging_theta' in self.attributes.keys():
                caseSCM.set_attribute('z_nudging_temp',self.attributes['z_nudging_theta'])
                if not('p_nudging_theta' in self.attributes.keys()):
                    height = np.squeeze(caseSCM.variables['height'].data)
                    pressure = np.squeeze(caseSCM.variables['pressure'].data)
                    plev = int(thermo.zlev2plev(self.attributes['z_nudging_theta'],height,pressure))
                    caseSCM.set_attribute('p_nudging_thetal',plev)
                    caseSCM.set_attribute('p_nudging_theta',plev)
                    caseSCM.set_attribute('p_nudging_temp',plev) 

        var = 'thetal_nudging'
        att = 'nudging_thetal'
        if att in self.attlist and self.attributes[att] > 0:
            if lflag:
                print 'Error: Several nudging variable for temperature are given, which might yield to inconsistencies'
                sys.exit()
            else:
                lflag=True   
            print '-'*10
            print 'thetal_nudging is given'
            caseSCM.add_variable(var,dataout[var].data,lev=lev,levtype=levtype,levid='lev',time=time,timeid='time')
            print 'assume theta_nudging=thetal_nudging'
            caseSCM.add_variable('theta_nudging',dataout[var].data,lev=lev,levtype=levtype,levid='lev',time=time,timeid='time')
            caseSCM.set_attribute('nudging_theta',self.attributes[att])
            print 'compute temp_nudging from thetal_nudging assuming no liquid water'
            pressure = caseSCM.variables['pressure_forc'].data
            thlnudg = caseSCM.variables['thetal_nudging'].data
            tnudg = thermo.theta2t(p=pressure,theta=thlnudg)
            caseSCM.add_variable('temp_nudging',tnudg,lev=lev,levtype=levtype,levid='lev',time=time,timeid='time')
            caseSCM.set_attribute('nudging_temp',self.attributes[att])
            if 'p_nudging_thetal' in self.attributes.keys():
                caseSCM.set_attribute('p_nudging_theta',self.attributes['p_nudging_thetal'])
                caseSCM.set_attribute('p_nudging_temp',self.attributes['p_nudging_thetal'])
                if not('z_nudging_thetal' in self.attributes.keys()):
                    height = np.squeeze(caseSCM.variables['height'].data)
                    pressure = np.squeeze(caseSCM.variables['pressure'].data)
                    zlev = int(thermo.plev2zlev(self.attributes['p_nudging_thetal'],height,pressure))
                    caseSCM.set_attribute('z_nudging_thetal',zlev)
                    caseSCM.set_attribute('z_nudging_theta',zlev)
                    caseSCM.set_attribute('z_nudging_temp',zlev)                
            if 'z_nudging_thetal' in self.attributes.keys():
                caseSCM.set_attribute('z_nudging_theta',self.attributes['z_nudging_thetal'])
                caseSCM.set_attribute('z_nudging_temp',self.attributes['z_nudging_thetal'])
                if not('p_nudging_thetal' in self.attributes.keys()):
                    height = np.squeeze(caseSCM.variables['height'].data)
                    pressure = np.squeeze(caseSCM.variables['pressure'].data)
                    plev = int(thermo.zlev2plev(self.attributes['z_nudging_thetal'],height,pressure))
                    caseSCM.set_attribute('p_nudging_thetal',plev)
                    caseSCM.set_attribute('p_nudging_theta',plev)
                    caseSCM.set_attribute('p_nudging_temp',plev) 

        #---- Humidity nudging

        lflag = False

        var = 'qv_nudging'
        att = 'nudging_qv'
        if att in self.attlist and self.attributes[att] > 0:
            lflag = True
            print '-'*10
            print 'qv_nudging is given'            
            caseSCM.add_variable(var,dataout[var].data,lev=lev,levtype=levtype,levid='lev',time=time,timeid='time')
            print 'assume qt_nudging=qv_nudging'
            caseSCM.add_variable('qt_nudging',dataout[var].data,lev=lev,levtype=levtype,levid='lev',time=time,timeid='time')
            caseSCM.set_attribute('nudging_qt',self.attributes[att])
            qvnudg = caseSCM.variables['qv_nudging'].data
            print 'compute rv_nudging from qv_nudging and assume rt_nudging=rv_nudging'
            rvnudg = thermo.qt2rt(qvnudg)
            caseSCM.add_variable('rv_nudging',rvnudg,lev=lev,levtype=levtype,levid='lev',time=time,timeid='time')
            caseSCM.set_attribute('nudging_rv',self.attributes[att])
            caseSCM.add_variable('rt_nudging',rvnudg,lev=lev,levtype=levtype,levid='lev',time=time,timeid='time')
            caseSCM.set_attribute('nudging_rt',self.attributes[att])
            if 'p_nudging_qv' in self.attributes.keys():
                caseSCM.set_attribute('p_nudging_qt',self.attributes['p_nudging_qv'])
                caseSCM.set_attribute('p_nudging_rv',self.attributes['p_nudging_qv'])
                caseSCM.set_attribute('p_nudging_rt',self.attributes['p_nudging_qv'])
                if not('z_nudging_qv' in self.attributes.keys()):
                    height = np.squeeze(caseSCM.variables['height'].data)
                    pressure = np.squeeze(caseSCM.variables['pressure'].data)
                    zlev = int(thermo.plev2zlev(self.attributes['p_nudging_qv'],height,pressure))
                    caseSCM.set_attribute('z_nudging_qv',zlev)
                    caseSCM.set_attribute('z_nudging_qt',zlev)
                    caseSCM.set_attribute('z_nudging_rv',zlev)
                    caseSCM.set_attribute('z_nudging_rt',zlev)                  
            if 'z_nudging_qv' in self.attributes.keys():
                caseSCM.set_attribute('z_nudging_qt',self.attributes['z_nudging_qv'])
                caseSCM.set_attribute('z_nudging_rv',self.attributes['z_nudging_qv'])
                caseSCM.set_attribute('z_nudging_rt',self.attributes['z_nudging_qv'])
                if not('p_nudging_qv' in self.attributes.keys()):
                    height = np.squeeze(caseSCM.variables['height'].data)
                    pressure = np.squeeze(caseSCM.variables['pressure'].data)
                    plev = int(thermo.zlev2plev(self.attributes['z_nudging_qv'],height,pressure))
                    caseSCM.set_attribute('p_nudging_qv',plev)
                    caseSCM.set_attribute('p_nudging_qt',plev)
                    caseSCM.set_attribute('p_nudging_rv',plev)
                    caseSCM.set_attribute('p_nudging_rt',plev)

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
            print 'assume qv_nudging=qt_nudging'
            caseSCM.add_variable('qv_nudging',dataout[var].data,lev=lev,levtype=levtype,levid='lev',time=time,timeid='time')
            caseSCM.set_attribute('nudging_qv',self.attributes[att])
            qtnudg = caseSCM.variables['qt_nudging'].data
            print 'compute rt_nudging from qt_nudging and assume rv_nudging=rt_nudging'
            rtnudg = thermo.qt2rt(qtnudg)
            caseSCM.add_variable('rt_nudging',rtnudg,lev=lev,levtype=levtype,levid='lev',time=time,timeid='time')
            caseSCM.set_attribute('nudging_rt',self.attributes[att])
            caseSCM.add_variable('rv_nudging',rtnudg,lev=lev,levtype=levtype,levid='lev',time=time,timeid='time')
            caseSCM.set_attribute('nudging_rv',self.attributes[att])
            if 'p_nudging_qt' in self.attributes.keys():
                caseSCM.set_attribute('p_nudging_qv',self.attributes['p_nudging_qt'])
                caseSCM.set_attribute('p_nudging_rv',self.attributes['p_nudging_qt'])
                caseSCM.set_attribute('p_nudging_rt',self.attributes['p_nudging_qt'])
                if not('z_nudging_qt' in self.attributes.keys()):
                    height = np.squeeze(caseSCM.variables['height'].data)
                    pressure = np.squeeze(caseSCM.variables['pressure'].data)
                    zlev = int(thermo.plev2zlev(self.attributes['p_nudging_qt'],height,pressure))
                    caseSCM.set_attribute('z_nudging_qv',zlev)
                    caseSCM.set_attribute('z_nudging_qt',zlev)
                    caseSCM.set_attribute('z_nudging_rv',zlev)
                    caseSCM.set_attribute('z_nudging_rt',zlev)                  
            if 'z_nudging_qt' in self.attributes.keys():
                caseSCM.set_attribute('z_nudging_qv',self.attributes['z_nudging_qt'])
                caseSCM.set_attribute('z_nudging_rv',self.attributes['z_nudging_qt'])
                caseSCM.set_attribute('z_nudging_rt',self.attributes['z_nudging_qt'])
                if not('p_nudging_qt' in self.attributes.keys()):
                    height = np.squeeze(caseSCM.variables['height'].data)
                    pressure = np.squeeze(caseSCM.variables['pressure'].data)
                    plev = int(thermo.zlev2plev(self.attributes['z_nudging_qt'],height,pressure))
                    caseSCM.set_attribute('p_nudging_qv',plev)
                    caseSCM.set_attribute('p_nudging_qt',plev)
                    caseSCM.set_attribute('p_nudging_rv',plev)
                    caseSCM.set_attribute('p_nudging_rt',plev)

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
            print 'assume rt_nudging=rv_nudging'
            caseSCM.add_variable('rt_nudging',dataout[var].data,lev=lev,levtype=levtype,levid='lev',time=time,timeid='time')
            caseSCM.set_attribute('nudging_rt',self.attributes[att])
            rvnudg = caseSCM.variables['rv_nudging'].data
            print 'compute qv_nudging from rv_nudging and assume qt_nudging=qv_nudging'
            qvnudg = thermo.rt2qt(rvnudg)
            caseSCM.add_variable('qv_nudging',qvnudg,lev=lev,levtype=levtype,levid='lev',time=time,timeid='time')
            caseSCM.set_attribute('nudging_qv',self.attributes[att])
            caseSCM.add_variable('qt_nudging',qvnudg,lev=lev,levtype=levtype,levid='lev',time=time,timeid='time')   
            caseSCM.set_attribute('nudging_qt',self.attributes[att])
            if 'p_nudging_rv' in self.attributes.keys():
                caseSCM.set_attribute('p_nudging_qv',self.attributes['p_nudging_rv'])
                caseSCM.set_attribute('p_nudging_qt',self.attributes['p_nudging_rv'])
                caseSCM.set_attribute('p_nudging_rt',self.attributes['p_nudging_rv'])
                if not('z_nudging_rv' in self.attributes.keys()):
                    height = np.squeeze(caseSCM.variables['height'].data)
                    pressure = np.squeeze(caseSCM.variables['pressure'].data)
                    zlev = int(thermo.plev2zlev(self.attributes['p_nudging_rv'],height,pressure))
                    caseSCM.set_attribute('z_nudging_qv',zlev)
                    caseSCM.set_attribute('z_nudging_qt',zlev)
                    caseSCM.set_attribute('z_nudging_rv',zlev)
                    caseSCM.set_attribute('z_nudging_rt',zlev)                  
            if 'z_nudging_rv' in self.attributes.keys():
                caseSCM.set_attribute('z_nudging_qv',self.attributes['z_nudging_rv'])
                caseSCM.set_attribute('z_nudging_qt',self.attributes['z_nudging_rv'])
                caseSCM.set_attribute('z_nudging_rt',self.attributes['z_nudging_rv'])
                if not('p_nudging_rv' in self.attributes.keys()):
                    height = np.squeeze(caseSCM.variables['height'].data)
                    pressure = np.squeeze(caseSCM.variables['pressure'].data)
                    plev = int(thermo.zlev2plev(self.attributes['z_nudging_rv'],height,pressure))
                    caseSCM.set_attribute('p_nudging_qv',plev)
                    caseSCM.set_attribute('p_nudging_qt',plev)
                    caseSCM.set_attribute('p_nudging_rv',plev)
                    caseSCM.set_attribute('p_nudging_rt',plev)

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
            print 'assume rv_nudging=rt_nudging'
            caseSCM.add_variable('rv_nudging',dataout[var].data,lev=lev,levtype=levtype,levid='lev',time=time,timeid='time')
            caseSCM.set_attribute('nudging_rv',self.attributes[att])
            rtnudg = caseSCM.variables['rt_nudging'].data
            print 'compute qt_nudging from rt_nudging and assume qv_nudging=qt_nudging'
            qtnudg = thermo.rt2qt(rtnudg)
            caseSCM.add_variable('qt_nudging',qtnudg,lev=lev,levtype=levtype,levid='lev',time=time,timeid='time')
            caseSCM.set_attribute('nudging_qt',self.attributes[att])
            caseSCM.add_variable('qv_nudging',qvnudg,lev=lev,levtype=levtype,levid='lev',time=time,timeid='time') 
            caseSCM.set_attribute('nudging_qv',self.attributes[att])
            if 'p_nudging_rt' in self.attributes.keys():
                caseSCM.set_attribute('p_nudging_qv',self.attributes['p_nudging_rt'])
                caseSCM.set_attribute('p_nudging_qt',self.attributes['p_nudging_rt'])
                caseSCM.set_attribute('p_nudging_rv',self.attributes['p_nudging_rt'])
                if not('z_nudging_rt' in self.attributes.keys()):
                    height = np.squeeze(caseSCM.variables['height'].data)
                    pressure = np.squeeze(caseSCM.variables['pressure'].data)
                    zlev = int(thermo.plev2zlev(self.attributes['p_nudging_rt'],height,pressure))
                    caseSCM.set_attribute('z_nudging_qv',zlev)
                    caseSCM.set_attribute('z_nudging_qt',zlev)
                    caseSCM.set_attribute('z_nudging_rv',zlev)
                    caseSCM.set_attribute('z_nudging_rt',zlev)                  
            if 'z_nudging_rt' in self.attributes.keys():
                caseSCM.set_attribute('z_nudging_qv',self.attributes['z_nudging_rt'])
                caseSCM.set_attribute('z_nudging_qt',self.attributes['z_nudging_rt'])
                caseSCM.set_attribute('z_nudging_rv',self.attributes['z_nudging_rt'])
                if not('p_nudging_rt' in self.attributes.keys()):
                    height = np.squeeze(caseSCM.variables['height'].data)
                    pressure = np.squeeze(caseSCM.variables['pressure'].data)
                    plev = int(thermo.zlev2plev(self.attributes['z_nudging_rt'],height,pressure))
                    caseSCM.set_attribute('p_nudging_qv',plev)
                    caseSCM.set_attribute('p_nudging_qt',plev)
                    caseSCM.set_attribute('p_nudging_rv',plev)
                    caseSCM.set_attribute('p_nudging_rt',plev)

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

