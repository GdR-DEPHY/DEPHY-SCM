#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on 27 November 2019

@author: Romain Roehrig
"""

import os
import sys
import time
import copy

from datetime import datetime 

import netCDF4 as nc

import numpy as np

from Axis import Axis
from Variable import Variable, read as readvar, interpol

from variables_attributes import attributes as var_attributes
from attributes import known_attributes, required_attributes

import thermo
import constants as CC

lwarnings = False

startDate0 = datetime(1979,1,1,0,0,0)
endDate0 = datetime(1979,1,1,0,0,0)

class Case:

    def __init__(self,caseid,lat=None,lon=None,startDate=startDate0,endDate=endDate0,surfaceType='ocean',zorog=0.):

        self.id = caseid

        self.set_dates(startDate,endDate)

        # Latitude (degrees_noth) and Longitude (degrees_east)
        self.lat = lat
        self.lon = lon

        # Surface type
        self.surface_type = surfaceType

        # Variables
        self.var_init_list = []
        self.var_forcing_list = []
        self.variables = {}

        # Attributes
        self.attlist = ['case','title','reference','author','version','format_version','modifications','script','comment',
                'start_date','end_date',
                'radiation',
                'surface_type','surface_forcing_temp','surface_forcing_moisture','surface_forcing_wind'
                ]
        self.attributes = {
                'case': self.id,
                'title': "",
                'reference': "",
                'author': "",
                'version': "Created on " + time.ctime(time.time()),
                'format_version': "DEPHY SCM format version 1",
                'modifications': "",
                'script': "",
                'comment': "",
                'start_date': self.start_date.strftime('%Y-%m-%d %H:%M:%S'),
                'end_date': self.end_date.strftime('%Y-%m-%d %H:%M:%S'),
                'radiation': 'on',
                'surface_type': self.surface_type,
                'surface_forcing_temp': 'none',
                'surface_forcing_moisture': 'none',
                'surface_forcing_wind': 'none',
                }

        for att in set(required_attributes).difference(set(self.attlist)):
            self.attlist.append(att)
            if att in ['surface_type','surface_forcing_temp','surface_forcing_moisture','surface_forcing_wind','radiation']:
                self.attributes[att] = ""
            else:
                self.attributes[att] = 0

        # Add latitude and longitude variables
        if not(lat is None):
            self.add_latitude(lat)
        if not(lon is None):
            self.add_longitude(lon)

        # Altitude above sea level (m)
        self.add_orography(zorog)


    def set_dates(self,startDate,endDate):

        # Start en end dates of the simulation
        if isinstance(startDate,datetime):
            self.start_date = startDate
        elif len(startDate) == 14:
            self.start_date = datetime.strptime(startDate,'%Y%m%d%H%M%S')
        else:
            self.start_date = datetime.strptime(startDate,'%Y-%m-%d %H:%M:%S')

        if isinstance(endDate,datetime):
            self.end_date = endDate
        elif len(endDate) == 14:
            self.end_date = datetime.strptime(endDate,'%Y%m%d%H%M%S')
        else:
            self.end_date = datetime.strptime(endDate,'%Y-%m-%d %H:%M:%S')
        
        # Initial time axis
        self.t0 = 0
        if not(startDate is None):
            self.tunits = self.start_date.strftime('seconds since %Y-%m-%d %H:%M:%S')
            self.t0Axis = Axis('t0',[self.t0,],name='initial_time',units=self.tunits,calendar='gregorian')

        # Convert startDate and endDate for time axes 
        self.tstart = None
        if startDate is not None:
            self.tstart = nc.date2num(self.start_date,self.tunits,calendar='gregorian') # should be 0
        self.tend = None
        if endDate is not None:
            self.tend = nc.date2num(self.end_date,self.tunits,calendar='gregorian')

    def set_latlon(self,lat,lon):

        self.lat = lat
        self.lon = lon

        self.add_latitude(lat)
        self.add_latitude(lon)

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

#    def set_z0(self,z0):
#
#        self.attlist.append('z0')
#        self.attributes['z0'] = z0

    def add_variable(self, varid, vardata, name=None, units=None,
            lev=None, levtype=None, levid=None,
            height=None, pressure=None,
            time=None, timeid=None):
        """Add a variable to a Case object.
            
        Required arguments:
        varid   -- string for the variable id. Should be in ...
        vardata -- input data as a list or a numpy array

        Optional (keyword) arguments:
        lev     -- input data for the level axis as a list, a numpy array or an Axis object (default None)
        levtype -- string describing the type of the level axis: None (default), 'altitude' or 'pressure'
        levid   -- string for the level axis id. The default is None, which implies a generic id (default None)
        height  -- variable object describing height for current variable
        pressure-- variable object describing pressure for current variable
        time    -- input data for the time axis, as a list, a numpy array or an Axis object (default None)
        timeid  -- string for the time axis id. The default is None, which implies a generic id (default None)
        name    -- string of the name attribute of the variable (to be use as long_name in a netCDF file) (default None)
        units   -- string of the units attribute of the variable (default None)
        """

        # if variable is already defined, stop
        # if not, add it to case variables dictionnary

        if varid in self.var_init_list + self.var_forcing_list:
            if lwarnings:
                print 'WARNING: Variable {0} is already defined'.format(varid)
                print 'WARNING: It will be overwritten'
            #sys.exit()
        else:
            if varid in ['ps','zh','pa','ua','va','ta','theta','thetal','qv','qt','rv','rt','rl','ri','ql','qi','tke','ts']:
                self.var_init_list.append(varid)
            else:
                self.var_forcing_list.append(varid)
        #print varid, name, lev, time

        ######################
        # Prepare time axis
        if time is None:
            timeAxis = None
            nt = None
            print 'ERROR'
            raise ValueError
        elif isinstance(time,Axis):
            timeAxis = time
            nt, = time.data.shape
        else:
            nt, = np.array(time).shape # In case time is given as a list
            # time is supposed to be given in seconds since beginning
            if timeid is None:
                timeAxis = Axis('time_{0}'.format(varid),time,
                                name='forcing_time_for_{0}'.format(varid),
                                units=self.tunits)
            else:
                timeAxis = Axis(timeid,time,name='forcing_time',units=self.tunits)

        ######################
        # Prepare level axis, if needed
        if lev is None:
            levAxis=None
            nlev = None
        elif isinstance(lev,Axis):
            levAxis = lev
            nlev, = lev.data.shape
            levdata = lev.data
        else:
            levdata = np.array(lev)
            nlev, = np.array(lev).shape # In case lev is given as a list
            if levtype == 'altitude':
                levunits = 'm'
                lev_name = 'height_for_{0}'.format(varid)
            elif levtype == 'pressure':
                levunits = 'Pa'
                lev_name = 'air_pressure_for_{0}'.format(varid)
            else:
                print 'ERROR: levtype unexpected:', levtype
                print 'ERROR: levtype should be defined and in altitude, pressure:'
                sys.exit()

            if levid is None:
                levAxis = Axis('lev_{0}'.format(varid),lev,name=lev_name,units=levunits)
            else:
                levAxis = Axis(levid,lev,name='{0}'.format(levtype),units=levunits)

        height_id = None
        height_units = None

        pressure_id = None
        pressure_units = None

        if levAxis is not None:
            if levtype == 'altitude':
                if height is None:
                    height = np.tile(levdata,(nt,1))
                    height_id = 'zh_{0}'.format(varid)
                    height_units = 'm'
            elif levtype == 'pressure':
                if pressure is None:
                    pressure = np.tile(levdata,(nt,1))
                    pressure_id = 'pa_{0}'.format(varid)
                    pressure_units = 'Pa'

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

        try:
            plotcoef = var_attributes[varid]['plotcoef']
            plotunits = var_attributes[varid]['plotunits']
        except KeyError:
            plotcoef = 1.
            plotunits = None
        except:
            raise

        ######################
        # Create variable
        if levAxis is None and timeAxis is None:
            print 'ERROR: level and time axes are None. Unexpected'
            sys.exit()
        else:
            if timeAxis is None:
                tmp = np.reshape(vardata,(nlev,))
            elif levAxis is None:
                tmp = np.reshape(vardata,(nt,))
            else:
                tmp = np.reshape(vardata,(nt,nlev,))

        self.variables[varid] = Variable(varid, data=tmp, name=varname, units=varunits,
                height=height, height_id=height_id, height_units=height_units,
                pressure=pressure, pressure_id=pressure_id, pressure_units=pressure_units,
                level=levAxis,time=timeAxis,
                plotcoef=plotcoef,plotunits=plotunits)


    def add_init_variable(self,varid,vardata,coordinate=False,**kwargs): 
        """Add an initial state variable to a case object.
            
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

        if varid in ['ps','ts']:
            # Put the expected shape of the input data
            tmp = np.reshape(vardata,(1,))
        else:
            # Check if lev optional argument is given
            if not(kwargs.has_key('lev')):
                print 'ERROR: level axis should be given for variable', varid
                sys.exit()

            if isinstance(kwargs['lev'],Axis):
                nlev, = kwargs['lev'].data.shape
            else:
                nlev, = np.array(kwargs['lev']).shape # In case lev is given as a list

            # Put the expected shape of the input data
            tmp = np.reshape(vardata,(1,nlev,))

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

        self.add_init_variable('zh',vardata,**kwargs)

    def add_init_pressure(self,vardata,**kwargs):
        """Add initial state variable for pressure to a Case object.
           
        Required argument:
        vardata -- input data as a list or a numpy array.

        See add_variable function for optional arguments.
        Note that:
        - a level axis is required (lev optional argument).
        - a levtype is required (levtype optional argument).
        """

        self.add_init_variable('pa',vardata,**kwargs)

    def add_init_temp(self,vardata,**kwargs):
        """Add initial state variable for temperature to a Case object.
           
        Required argument:
        vardata -- input data as a list or a numpy array.

        See add_variable function for optional arguments.
        Note that:
        - a level axis is required (lev optional argument).
        - a levtype is required (levtype optional argument).
        """

        self.add_init_variable('ta',vardata,**kwargs)

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
        self.add_init_variable('ua',u,**kwargs)

        if vlev is not None:
            kwargs['lev'] = vlev
        self.add_init_variable('va',v,**kwargs)

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

    def add_init_ts(self,vardata,**kwargs):
        """Add initial state variable for surface temperature to a Case object.
           
        Required argument:
        vardata -- input data as an integer or a float.

        See add_variable function for optional arguments.
        """

        self.add_init_variable('ts',vardata,**kwargs)

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

        if varid in ['ps_forc','hfss','hfls','ustar','ts_forc','tskin','orog','lat','lon','z0','z0h','z0q','beta']:
            # Put the expected shape of the input data
            if lconstant: 
                tmp = np.zeros((nt),dtype=np.float32)
                tmp[0] = vardata
                tmp[1] = vardata
            else:
                tmp = np.reshape(vardata,(nt,))
        else:
            # Check if lev optional argument is given
            if not(kwargs.has_key('lev')):
                print 'ERROR: level axis should be given for variable', varid
                sys.exit()

            nlev, = np.array(kwargs['lev']).shape # In case lev is given as a list

            # Put the expected shape of the input data
            if lconstant:
                tmp = np.zeros((nt,nlev),dtype=np.float32)
                tmp[0,:] = vardata[:]
                tmp[1,:] = vardata[:]
            else:
                tmp = np.reshape(vardata,(nt,nlev))

        # add initial variable to the Case object
        self.add_variable(varid,tmp,**kwargs)

    def add_latitude(self,data,**kwargs):
        """Add latitude to a Case object
           
        Required argument:
        data -- input data as a list or a numpy array.

        See add_variable function for optional arguments.

        If time is not provided, forcing is assumed constant in time.           
        """

        self.add_forcing_variable('lat',data,**kwargs)

    def add_longitude(self,data,**kwargs):
        """Add longitude to a Case object
           
        Required argument:
        data -- input data as a list or a numpy array.

        See add_variable function for optional arguments.

        If time is not provided, forcing is assumed constant in time.           
        """

        self.add_forcing_variable('lon',data,**kwargs)

    def add_orography(self,data,**kwargs):
        """Add orography to a Case object
        
        Required argument:
        data -- input data as a list or a numpy array.

        See add_variable function for optional arguments.

        If time is not provided, forcing is assumed constant in time.           
        """

        self.add_forcing_variable('orog',data,**kwargs)

    def add_z0(self,data,**kwargs):
        """Add roughness length for wind to a Case object
        
        Required argument:
        data -- input data as a list or a numpy array.

        See add_variable function for optional arguments.

        If time is not provided, forcing is assumed constant in time.           
        """

        self.add_forcing_variable('z0',data,**kwargs)

    def add_z0h(self,data,**kwargs):
        """Add roughness length for heat to a Case object
        
        Required argument:
        data -- input data as a list or a numpy array.

        See add_variable function for optional arguments.

        If time is not provided, forcing is assumed constant in time.           
        """

        self.add_forcing_variable('z0h',data,**kwargs)

    def add_z0q(self,data,**kwargs):
        """Add roughness length for moisture to a Case object
        
        Required argument:
        data -- input data as a list or a numpy array.

        See add_variable function for optional arguments.

        If time is not provided, forcing is assumed constant in time.           
        """

        self.add_forcing_variable('z0q',data,**kwargs)

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

        self.add_forcing_variable('pa_forc',data,**kwargs)

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

        self.add_forcing_variable('zh_forc',data,**kwargs)

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
            self.set_attribute('forc_wa',1)
            self.add_forcing_variable('wa',w,**kwargs)

        if omega is not None:
            self.set_attribute('forc_wap',1)
            self.add_forcing_variable('wap',omega,**kwargs)

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

        self.set_attribute('adv_ta',1)
        if include_rad:
            self.set_attribute('radiation','off')

        self.add_forcing_variable('tnta_adv',data,**kwargs)

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
            self.set_attribute('radiation','off')

        self.add_forcing_variable('tntheta_adv',data,**kwargs)

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
            self.set_attribute('radiation','off')

        self.add_forcing_variable('tnthetal_adv',data,**kwargs)

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

        self.add_forcing_variable('tnqv_adv',data,**kwargs)

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

        self.add_forcing_variable('tnqt_adv',data,**kwargs)

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

        self.add_forcing_variable('tnrv_adv',data,**kwargs)

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

        self.add_forcing_variable('tnrt_adv',data,**kwargs)

    def add_nudging(self,varid,data,timescale=None,z_nudging=None,p_nudging=None,**kwargs):
        """Add a nudging forcing to a Case object.
        
        Required argument:
        varid     -- id of the variable to be nudged as a string
        data      -- input data as a list or a numpy array.
        timescale -- nudging timescale in seconds (integer or float)
        z_nudging -- altitude above which nudging is applied (integer or float)
        p_nudging -- pressure altitude under which nudging is applied (integer or float)

        Either z_nudging, p_nudging or nudging_profile must be defined.

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
            self.set_attribute('za_nudging_{0}'.format(varid),0)

        if z_nudging is not None:
            self.set_attribute('zh_nudging_{0}'.format(varid),float(z_nudging))
        else:
            self.set_attribute('pa_nudging_{0}'.format(varid),float(p_nudging))

        var = '{0}_nud'.format(varid)
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
        self.add_nudging('ua',unudg,**kwargs)

        if vlev is not None:
            kwargs['lev'] = vlev        
        self.add_nudging('va',vnudg,**kwargs)

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

        self.add_nudging('ta',data,**kwargs)

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

    def add_ozone(self,ro3=None,o3=None,**kwargs):
        """Add an ozone forcing to a Case object.

        
        One of the following argument is required:
        ro3 -- ozone mixing ratio as a list or a numpy array.
        o3 -- ozone mole fraction in air as a list or a numpy array.

        If time is not provided, forcing is assumed constant in time.
        """

        if ro3 is not None:
            self.add_forcing_variable('o3',ro3*CC.Md/CC.Mo3,**kwargs)
        elif o3 is not None:
            self.add_forcing_variable('o3',o3,**kwargs)
        else:
            raise ValueError('Either ozone mixing ratio or ozone mole fraction in air must be provided')

    def activate_radiation(self):
        """Activate radiation in a Case object
           
        No argument required.
        """

        self.set_attribute('radiation'.format(var),"on")

    def deactivate_radiation(self):
        """Deactivate radiation in a Case object
           
        No argument required.
        """

        self.set_attribute('radiation',"off")

    def add_surface_temp(self,data,**kwargs):
        """Add a surface temperature forcing to a Case object.
        
        This function does not imply that the surface forcing type is ts.
        This function can be used to add a useful surface temperature in case surface fluxes are prescribed.

        Required argument:
        data -- input data as a list or a numpy array.

        If time is not provided, forcing is assumed constant in time.
        """

        self.add_forcing_variable('ts_forc',data,**kwargs)

    def add_surface_skin_temp(self,data,**kwargs):
        """Add a surface skin temperature forcing to a Case object.
        
        This function does not imply that the surface forcing type is ts.
        This function can be used to add a useful surface temperature in case surface fluxes are prescribed.

        Required argument:
        data -- input data as a list or a numpy array.

        If time is not provided, forcing is assumed constant in time.
        """

        self.add_forcing_variable('tskin',data,**kwargs)

    def add_forcing_ts(self,data,z0=None,**kwargs):
        """Add a surface temperature forcing to a Case object.
        
        This function sets a surface temperature forcing as the case surface forcing.

        Required argument:
        data -- input data as a list or a numpy array.

        If time is not provided, forcing is assumed constant in time.
        """

        self.set_attribute('surface_forcing_temp','ts')
        self.add_forcing_variable('ts_forc',data,**kwargs)

        if z0 is not None:
            self.set_attribute('surface_forcing_wind','z0')
            self.add_forcing_variable('z0',z0)
        

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

        self.set_attribute('surface_forcing_temp','surface_flux')
        self.set_attribute('surface_forcing_moisture','surface_flux')

        if forc_wind == 'z0':
            self.set_attribute("surface_forcing_wind","z0")
            if z0 is None:
                print 'ERROR: z0 must be provided'
                sys.exit()
            self.add_forcing_variable('z0',z0)
            #self.set_z0(z0)
        elif forc_wind == 'ustar':
            self.set_attribute("surface_forcing_wind","ustar")
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
            self.add_forcing_variable('hfss',sens,**kwargs)
        else:
            self.add_forcing_variable('hfss',sens,**kwargs)

        if time_lat is not None:
            kwargs['time'] = time_lat
            self.add_forcing_variable('hfls',lat,**kwargs)
        else:
            self.add_forcing_variable('hfls',lat,**kwargs)

    def set_betaevap(self,beta=1.):
        """Activate a beta model for surface evaporation in a Case object
           
        Optional (keyword) argument:
        beta -- beta value of the beta model (default: 1.)
        """

        self.set_attribute("surface_forcing_moisture","beta")
        self.add_forcing_variable('beta',beta)

    def deactivate_surface_evaporation(self):
        """Deactivate surface evoporation in a Case object
           
        No argument required.
        """

        self.set_attribute("surface_forcing_moisture","beta")
        self.add_forcing_variable('beta',0.)

    def info(self):
        """Print information about a case object.

        No argument required.
        """

        for att in known_attributes:
            if att in self.attlist:
                print '{0}: {1}'.format(att,self.attributes[att])
     
        print "######################"
        print "# Variable information"
        for var in self.var_init_list + self.var_forcing_list:
            self.variables[var].info()

    def write(self,fileout,verbose=False):
        """Write case object into a netCDF file

        Required argument:
        fileout -- netCDF file name as a string

        Time axes are first written, then level axes, 
        then altitude/pressure variables 
        and finally the variables themselves.
        """

        g = nc.Dataset(fileout,'w',format='NETCDF3_CLASSIC')

        # Writing first only time axes
        for var in var_attributes.keys():
            if var in self.var_init_list + self.var_forcing_list:
                if verbose: print 'writing time axis for', var
                self.variables[var].write(g,
                        write_time_axes=True, write_level_axes=False,
                        write_data=False, write_vertical=False)

        # Writing first only level axes
        for var in var_attributes.keys():
            if var in self.var_init_list + self.var_forcing_list:
                if verbose: print 'writing level axis for', var
                self.variables[var].write(g,
                        write_time_axes=False, write_level_axes=True,
                        write_data=False, write_vertical=False)

        # Writing then only vertical variables
        for var in var_attributes.keys():
            if var in self.var_init_list + self.var_forcing_list:
                if verbose: print 'writing vertical variable for', var
                self.variables[var].write(g, 
                        write_time_axes=False, write_level_axes=False,
                        write_data=False, write_vertical=True)

        # Finally write data
        for var in var_attributes.keys():
            if var in self.var_init_list + self.var_forcing_list:
                if verbose: print 'writing data for', var
                self.variables[var].write(g,
                        write_time_axes=False, write_level_axes=False,
                        write_data=True, write_vertical=False)

        # Write global attributes
        for att in known_attributes:
            if att in self.attlist:
                g.setncattr(att,self.attributes[att])

        g.close()

    def read(self,filein,verbose=False):
        """Read a netCDF file to store data and information into a case object.

        Required argument:
        filein -- netCDF file name as a string

        Note that the case object should be initialize first.
        """

        f = nc.Dataset(filein,'r')

        self.set_dates(f.start_date,f.end_date)
        #self.set_latlon(f['lat'][0],f['lon'][0])

        for var in f.variables:
            if (var not in f.dimensions and var[0:6] != 'bounds' and var[0:3] not in ['zh_','pa_'])\
                    or var in ['zh_forc','pa_forc']:
                if verbose: print 'Reading', var
                tmp = readvar(var,f)
                if tmp.level is None:
                    self.add_variable(var, tmp.data,
                            time=tmp.time)
                elif tmp.level.units == 'm':
                    self.add_variable(var, tmp.data,
                            time=tmp.time,
                            lev=tmp.level, levtype='altitude')
                elif tmp.level.units == 'Pa':
                    self.add_variable(var, tmp.data,
                            time=tmp.time,
                            lev=tmp.level, levtype='pressure')                    
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

        for var in self.var_init_list + self.var_forcing_list:
            self.variables[var].plot(rep_images=rep_images,timeunits=timeunits,levunits=levunits)

    def plot_compare(self,cc,rep_images='./images/',label1=None,label2=None,timeunits=None,levunits=None):

        if not(os.path.exists(rep_images)):
            os.makedirs(rep_images)

        for var in self.var_init_list + self.var_forcing_list:
            if var in cc.var_init_list + cc.var_forcing_list: 
                if len(self.variables[var].data.shape) <= 3  or self.variables[var].data.shape[0] == 1:
                    self.variables[var].plot(rep_images=rep_images,
                            var2=cc.variables[var],
                            label=label1,label2=label2,
                            timeunits=timeunits,levunits=levunits)

#    def compute_height(self):
#        """Compute height for a case object
#
#        No input argument.
#
#        Return a Variable object for the height
#        """
#
#        pressure = self.variables['pa']
#        if 'ta' in self.variables.keys():
#            temp = self.variables['ta'].data
#            height = thermo.p2z(temp=temp,p=pressure)
#        elif 'theta' in self.variables.keys():
#            theta = self.variables['theta'].data
#            height = thermo.p2z(theta=theta,p=pressure)
#        elif 'thetal' in self.variables.keys():
#            thetal = self.variables['thetal'].data
#            height = thermo.p2z(theta=thetal,p=pressure)
#        else:
#            print 'ERROR: ta, theta or thetal should be provided'
#            raise ValueError
#
#        return Variable('zh', data=height, 
#                name=var_attributes['zh']['name'], units=var_attributes['zh']['units'],
#                #height=height, height_id=height_id, height_units=height_units,
#                #pressure=pressure, pressure_id=pressure_id, pressure_units=pressure_units,
#                level=self.variables['zh'].level,time=self.variables['zh'].time)
#
#    def compute_pressure(self):
#        """Compute pressure for a case object
#
#        No input argument.
#
#        Return a Variable object for the pressure
#        """
#
#        ps = self.variables['ps'].data[0]
#
#        if 'qv' in self.variables.keys():
#            humidity = self.variables['qv'].data[0,:]
#        elif 'qt' in self.variables.keys():
#            humidity = self.variables['qt'].data[0,:]
#        else:
#            humidity = None
#
#        if 'theta' in self.variables.keys():
#            z = self.variables['theta'].level.data
#            theta = self.variables['theta'].data[0,:]
#            pressure = thermo.z2p(theta=theta,z=z,ps=ps,qv=humidity)
#        elif 'thetal' in self.variables.keys():
#            z = self.variables['thetal'].level.data
#            thetal = self.variables['thetal'].data[0,:]
#            pressure = thermo.z2p(theta=thetal,z=z,ps=ps,qv=humidity)
#        elif 'temp' in self.variables.keys():
#            z = self.variables['ta'].level.data
#            temp = self.variables['ta'].data[0,:,0,0]
#            pressure = thermo.z2p(temp=temp,z=z,ps=ps,qv=humidity)
#        else:
#            print 'ERROR: theta, thetal or temp should be defined'
#            sys.exit()
#
#        return Variable('pa', data=pressure, 
#                name=var_attributes['pa']['name'], units=var_attributes['pa']['units'],
#                #height=height, height_id=height_id, height_units=height_units,
#                #pressure=pressure, pressure_id=pressure_id, pressure_units=pressure_units,
#                level=self.variables['pa'].level,time=self.variables['pa'].time)


    def compute_theta(self,pressure=None):

        if 'ta' in self.var_init_list:
            if pressure is None:
                print 'ERROR: pressure should be None to theta from ta'
                raise ValueError
            else:
                 print 'compute potential temperature from pressure and temperature'
                 theta = thermo.t2theta(p=pressure.data[0,:], temp=self.variables['ta'].data[0,:])
        elif 'thetal' in self.var_init_list:
            print 'assume theta=thetal'
            theta = self.variables['thetal'].data[0,:]
        else:
            print 'ERROR: At least ta or thetal should be given to compute theta'
            raise ValueError

        return theta

    def compute_thetal(self,pressure=None):

        if 'ta' in self.var_init_list:
            if pressure is None:
                print 'ERROR: pressure should be None to compute thetal from ta'
                raise ValueError
            else:
                print 'compute potential temperature from pressure and temperature, assuming thetal=theta'
                thetal = thermo.t2theta(p=pressure.data[0,:], temp=self.variables['ta'].data[0,:])
        elif 'theta' in self.var_init_list:
            print 'assume thetal=theta'
            thetal = self.variables['theta'].data[0,:]
        else:
            print 'ERROR: At least ta or theta should be given to compute thetal'
            raise ValueError

        return thetal

    def compute_temp(self,pressure=None):

        if pressure is None:
            print 'ERROR: pressure should be None to compute ta from theta or thetahl'
            raise ValueError

        if 'theta' in self.var_init_list:
            print 'compute temperature from pressure and potential temperature'
            temp = thermo.theta2t(p=pressure.data[0,:], theta=self.variables['theta'].data[0,:])
        elif 'thetal' in self.var_init_list:
            print 'compute temperature from pressure and liquid potential temperature (No liquid water considered)'
            temp = thermo.theta2t(p=pressure.data[0,:], theta=self.variables['thetal'].data[0,:])
        else:
            print 'ERROR: At least theta or thetal should be given'
            raise ValueError

        return temp

    def compute_qv(self):

        if 'qt' in self.var_init_list:
            print 'assume qv=qt'
            qv = self.variables['qt'].data[0,:]
        elif 'rv' in self.var_init_list:
            print 'compute qv from rv'
            qv = thermo.rt2qt(self.variables['rv'].data[0,:])
        elif 'rt' in self.var_init_list:
            print 'compute qt from rt and assume qv=qt'
            qv = thermo.rt2qt(self.variables['rt'].data[0,:])
        else:
            print 'ERROR: Either qt, rv or rt should be defined to compute qv'
            raise ValueError

        return qv

    def compute_qt(self):

        if 'qv' in self.var_init_list:
            print 'assume qt=qv'
            qt = self.variables['qv'].data[0,:]
        elif 'rv' in self.var_init_list:
            print 'compute qv from rv and assume qt=qv'
            qt = thermo.rt2qt(self.variables['rv'].data[0,:])
        elif 'rt' in self.var_init_list:
            print 'compute qt from rt'
            qt = thermo.rt2qt(self.variables['rt'].data[0,:])
        else:
            print 'ERROR: Either qv, rv or rt should be defined to compute qt'
            raise ValueError

        return qt

    def compute_rv(self):

        if 'qv' in self.var_init_list:
            print 'compute rv from qv'
            rv = thermo.qt2rt(self.variables['qv'].data[0,:])
        elif 'qt' in self.var_init_list:
            print 'compute rt from qt and assume rv=rt'
            rv = thermo.qt2rt(self.variables['qt'].data[0,:])
        elif 'rt' in self.var_init_list:
            print 'assume rv=rt'
            rv = selv.variables['rt'].data[0,:]
        else:
            print 'ERROR: Either qv, qt or rt should be defined to compute rv'
            raise ValueError

        return rv

    def compute_rt(self):

        if 'qv' in self.var_init_list:
            print 'compute rt from qt, assuming qt=av'
            rt = thermo.qt2rt(self.variables['qv'].data[0,:])
        elif 'qt' in self.var_init_list:
            print 'compute rt from qt'
            rt = thermo.qt2rt(self.variables['qt'].data[0,:])
        elif 'rv' in self.var_init_list:
            print 'assume rt=rv'
            rt = selv.variables['rv'].data[0,:]
        else:
            print 'ERROR: Either qv, qt or rv should be defined to compute rt'
            raise ValueError

        return rt

    def compute_tnta_adv(self):

        pressure = self.variables['pa_forc'].data

        if 'tntheta_adv' in self.var_forcing_list:
            print 'compute tnta_adv from tntheta_adv'
            thadv = self.variables['tntheta_adv'].data
            tadv = thermo.theta2t(p=pressure,theta=thadv)
        elif 'tnthetal_adv' in self.var_forcing_list:
            print 'assume tnthetal_adv=tntheta_adv and compute tnta_adv from tntheta_adv'
            thladv = self.variables['tnthetal_adv'].data
            tadv = thermo.theta2t(p=pressure,theta=thladv)
        else:
            print 'ERROR: To compute tnta_adv, tntheta_adv or tnthetal_adv must be known'
            raise ValueError

        return tadv

    def compute_tntheta_adv(self):

        pressure = self.variables['pa_forc'].data

        if 'tnta_adv' in self.var_forcing_list:
            print 'compute tntheta_adv from tnta_adv'
            tadv = self.variables['tnta_adv'].data
            thadv = thermo.t2theta(p=pressure,temp=tadv)
        elif 'tnthetal_adv' in self.var_forcing_list:
            print 'assume tntheta_adv=tnthetal_adv'
            thadv = self.variables['tnthetal_adv'].data
        else:
            print 'ERROR: To compute tntheta_adv, tnta_adv or tnthetal_adv must be known'
            raise ValueError

        return thadv

    def compute_tnthetal_adv(self):

        pressure = self.variables['pa_forc'].data

        if 'tnta_adv' in self.var_forcing_list:
            print 'compute tnthetal_adv from tnta_adv assuming tnthetal_adv=tntheta_adv'
            tadv = self.variables['tnta_adv'].data
            thladv = thermo.t2theta(p=pressure,temp=tadv)
        elif 'tntheta_adv' in self.var_forcing_list:
            print 'assume tnthetal_adv=tntheta_adv'
            thladv = self.variables['tntheta_adv'].data
        else:
            print 'ERROR: To compute tnthetal_adv, tnta_adv or tntheta_adv must be known'
            raise ValueError

        return thladv

    def compute_tnta_rad(self):

        pressure = self.variables['pa_forc'].data

        if 'tntheta_rad' in self.var_forcing_list:
            print 'compute tnta_rad from tntheta_rad'
            thrad = self.variables['tntheta_rad'].data
            trad = thermo.theta2t(p=pressure,theta=thrad)
        elif 'tnthetal_rad' in self.var_forcing_list:
            print 'assume tnthetal_rad=tntheta_rad and compute tnta_rad from tntheta_rad'
            thlrad = self.variables['tnthetal_rad'].data
            trad = thermo.theta2t(p=pressure,theta=thlrad)
        else:
            print 'ERROR: To compute tnta_rad, tntheta_rad or tnthetal_rad must be known'
            raise ValueError

        return trad

    def compute_tntheta_rad(self):

        pressure = self.variables['pa_forc'].data

        if 'tnta_rad' in self.var_forcing_list:
            print 'compute tntheta_rad from tnta_rad'
            trad = self.variables['tnta_rad'].data
            thrad = thermo.t2theta(p=pressure,temp=trad)
        elif 'tnthetal_rad' in self.var_forcing_list:
            print 'assume tntheta_rad=tnthetal_rad'
            thrad = self.variables['tnthetal_rad'].data
        else:
            print 'ERROR: To compute tntheta_rad, tnta_rad or tnthetal_rad must be known'
            raise ValueError

        return thrad

    def compute_tnthetal_rad(self):

        pressure = self.variables['pa_forc'].data

        if 'tnta_rad' in self.var_forcing_list:
            print 'compute tnthetal_rad from tnta_rad assuming tnthetal_rad=tntheta_rad'
            trad = self.variables['tnta_rad'].data
            thlrad = thermo.t2theta(p=pressure,temp=trad)
        elif 'tntheta_rad' in self.var_forcing_list:
            print 'assume tnthetal_rad=tntheta_rad'
            thlrad = self.variables['tntheta_rad'].data
        else:
            print 'ERROR: To compute tnthetal_rad, tnta_rad or tntheta_rad must be known'
            raise ValueError

        return thlrad

    def compute_tnqv_adv(self):

        if 'tnqt_adv' in self.var_forcing_list:
            print 'assume tnqv_adv=tnqt_adv'
            qvadv = self.variables['tnqt_adv'].data
        elif 'tnrv_adv' in self.var_forcing_list:
            print 'compute tnqv_adv from tnrv_adv using initial rv profile'
            rvadv = self.variables['tnrv_adv'].data
            rv = self.variables['rv'].data
            qvadv =  thermo.advrt2advqt(rt=rv,advrt=rvadv)
        elif 'tnrt_adv' in self.var_forcing_list:
            print 'compute tnqv_adv from tnrv_adv assuming tnrv_adv=tnrt_adv and using initial rt profile'
            rtadv = self.variables['tnrt_adv'].data
            rt = self.variables['rt'].data
            qvadv =  thermo.advrt2advqt(rt=rt,advrt=rtadv)
        else:
            print 'ERROR: To compute tnqv_adv, tnqt_adv, tnrv_adv or tnrt_adv must be known'
            raise ValueError

        return qvadv

    def compute_tnqt_adv(self):

        if 'tnqv_adv' in self.var_forcing_list:
            print 'assume tnqt_adv=tnqv_adv'
            qtadv = self.variables['tnqv_adv'].data
        elif 'tnrv_adv' in self.var_forcing_list:
            print 'compute tnqt_adv from tnrt_adv assuming tnrt_adv=tnrv_adv using initial rv profile'
            rtadv = self.variables['tnrv_adv'].data
            rt = self.variables['rv'].data
            qtadv =  thermo.advrt2advqt(rt=rt,advrt=rtadv)
        elif 'tnrt_adv' in self.var_forcing_list:
            print 'compute tnqt_adv from tnrt_adv using initial rt profile'
            rtadv = self.variables['tnrt_adv'].data
            rt = self.variables['rt'].data
            qtadv =  thermo.advrt2advqt(rt=rt,advrt=rtadv)
        else:
            print 'ERROR: To compute tnqt_adv, tnqv_adv, tnrv_adv or tnrt_adv must be known'
            raise ValueError

        return qtadv

    def compute_tnrv_adv(self):

        if 'tnrt_adv' in self.var_forcing_list:
            print 'assume tnrv_adv=tnrt_adv'
            rvadv = self.variables['tnrt_adv'].data
        elif 'tnqv_adv' in self.var_forcing_list:
            print 'compute tnrv_adv from tnqv_adv using initial qv profile'
            qvadv = self.variables['tnqv_adv'].data
            qv = self.variables['qv'].data
            rvadv =  thermo.advqt2advrt(qt=qv,advqt=qvadv)
        elif 'tnqt_adv' in self.var_forcing_list:
            print 'compute tnrv_adv from tnqv_adv assuming tnqv_adv=tnqt_adv and using initial qt profile'
            qvadv = self.variables['tnqt_adv'].data
            qv = self.variables['qt'].data
            rvadv =  thermo.advqt2advrt(qt=qv,advqt=qvadv)
        else:
            print 'ERROR: To compute tnrv_adv, tnqv_adv, tnqt_adv or tnrt_adv must be known'
            raise ValueError

        return rvadv

    def compute_tnrt_adv(self):

        if 'tnrv_adv' in self.var_forcing_list:
            print 'assume tnrt_adv=tnrv_adv'
            rtadv = self.variables['tnrv_adv'].data
        elif 'tnqv_adv' in self.var_forcing_list:
            print 'compute tnrt_adv from tnqt_adv assuming tnqt_adv=tnqv_adv and using initial qv profile'
            qtadv = self.variables['tnqv_adv'].data
            qt = self.variables['qv'].data
            rtadv =  thermo.advqt2advrt(qt=qt,advqt=qtadv)
        elif 'tnqt_adv' in self.var_forcing_list:
            print 'compute tnrt_adv from tnqt_adv using initial qt profile'
            qtadv = self.variables['tnqt_adv'].data
            qt = self.variables['qt'].data
            rtadv =  thermo.advqt2advrt(qt=qt,advqt=qtadv)
        else:
            print 'ERROR: To compute tnrt_adv, tnqv_adv, tnqt_adv or tnrv_adv must be known'
            raise ValueError

        return rtadv

    def compute_ta_nud(self):

        pressure = self.variables['pa_forc'].data

        if 'theta_nud' in self.var_forcing_list:
            print 'compute ta_nud from theta_nud'
            thnud = self.variables['theta_nud'].data
            tnud = thermo.theta2t(p=pressure,theta=thnud)
        elif 'thetal_nud' in self.var_forcing_list:
            print 'compute ta_nud from theta_nud assuming theta_nud=thetal_nud'
            thnud = self.variables['thetal_nud'].data
            tnud = thermo.theta2t(p=pressure,theta=thnud)
        else:
            print 'ERROR: To compute ta_nud, theta_nud or thetal_nud must be known'
            raise ValueError

        return tnud

    def compute_ta_nud(self):

        pressure = self.variables['pa_forc'].data

        if 'theta_nud' in self.var_forcing_list:
            print 'compute ta_nud from theta_nud'
            thnud = self.variables['theta_nud'].data
            tnud = thermo.theta2t(p=pressure,theta=thnud)
        elif 'thetal_nud' in self.var_forcing_list:
            print 'compute ta_nud from theta_nud assuming theta_nud=thetal_nud'
            thnud = self.variables['thetal_nud'].data
            tnud = thermo.theta2t(p=pressure,theta=thnud)
        else:
            print 'ERROR: To compute ta_nud, theta_nud or thetal_nud must be known'
            raise ValueError

        return tnud

    def compute_theta_nud(self):

        pressure = self.variables['pa_forc'].data

        if 'ta_nud' in self.var_forcing_list:
            print 'compute theta_nud from ta_nud'
            tnud = self.variables['ta_nud'].data
            thnud = thermo.t2theta(p=pressure,temp=tnud)
        elif 'thetal_nud' in self.var_forcing_list:
            print 'assume theta_nud=thetal_nud'
            thnud = self.variables['thetal_nud'].data
        else:
            print 'ERROR: To compute theta_nud, ta_nud or thetal_nud must be known'
            raise ValueError

        return thnud

    def compute_thetal_nud(self):

        pressure = self.variables['pa_forc'].data

        if 'ta_nud' in self.var_forcing_list:
            print 'compute thetal_nud from ta_nud assuming thetal_nud=theta_nud'
            tnud = self.variables['ta_nud'].data
            thlnud = thermo.t2theta(p=pressure,temp=tnud)
        elif 'theta_nud' in self.var_forcing_list:
            print 'assume thetal_nud=theta_nud'
            thlnud = self.variables['theta_nud'].data
        else:
            print 'ERROR: To compute thetal_nud, ta_nud or theta_nud must be known'
            raise ValueError

        return thlnud

    def compute_qv_nud(self):

        if 'qt_nud' in self.var_forcing_list:
            print 'assume qv_nud=qt_nud'
            qvnud = self.variables['qt_nud'].data
        elif 'rv_nud' in self.var_forcing_list:
            print 'compute qv_nud from rv_nud'
            rvnud = self.variables['rv_nud'].data
            qvnud =  thermo.rt2qt(rvnud)
        elif 'rt_nud' in self.var_forcing_list:
            print 'compute qv_nud from rv_nud assuming rv_nud=rt_nud'
            rvnud = self.variables['rt_nud'].data
            qvnud =  thermo.rt2qt(rvnud)
        else:
            print 'ERROR: To compute qv_nud, qt_nud, rv_nud or rt_nud must be known'
            raise ValueError

        return qvnud

    def compute_qt_nud(self):

        if 'qv_nud' in self.var_forcing_list:
            print 'assume qt_nud=qv_nud'
            qtnud = self.variables['qv_nud'].data
        elif 'rv_nud' in self.var_forcing_list:
            print 'compute qt_nud from rt_nud assuming rt_nud=rv_nud'
            rtnud = self.variables['rv_nud'].data
            qtnud =  thermo.rt2qt(rtnud)
        elif 'rt_nud' in self.var_forcing_list:
            print 'compute qt_nud from rt_nud'
            rtnud = self.variables['rt_nud'].data
            qtnud =  thermo.rt2qt(rtnud)
        else:
            print 'ERROR: To compute qt_nud, qv_nud, rv_nud or rt_nud must be known'
            raise ValueError

        return qtnud

    def compute_rv_nud(self):

        if 'qv_nud' in self.var_forcing_list:
            print 'compute rv_nud from qv_nud'
            qvnud = self.variables['qv_nud'].data
            rvnud =  thermo.qt2rt(qvnud)
        elif 'qt_nud' in self.var_forcing_list:
            print 'compute rv_nud from qv_nud assuming qv_nud=qt_nud'
            qvnud = self.variables['qt_nud'].data
            rvnud =  thermo.qt2rt(qvnud)
        elif 'rt_nud' in self.var_forcing_list:
            print 'assume rv_nud=rt_nud'
            rvnud = self.variables['rt_nud'].data
        else:
            print 'ERROR: To compute rv_nud, qv_nud, qt_nud or rt_nud must be known'
            raise ValueError

        return rvnud

    def compute_rt_nud(self):

        if 'qv_nud' in self.var_forcing_list:
            print 'compute rt_nud from qt_nud assuming qt_nud=qv_nud'
            qtnud = self.variables['qv_nud'].data
            rtnud =  thermo.qt2rt(qtnud)
        elif 'qt_nud' in self.var_forcing_list:
            print 'compute rt_nud from qt_nud'
            qtnud = self.variables['qt_nud'].data
            rtnud =  thermo.qt2rt(qtnud)
        elif 'rv_nud' in self.var_forcing_list:
            print 'assume rt_nud=rv_nud'
            rtnud = self.variables['rv_nud'].data
        else:
            print 'ERROR: To compute rt_nud, qv_nud, qt_nud or rv_nud must be known'
            raise ValueError

        return rtnud




    def interpolate(self,time=None,lev=None,levtype=None,usetemp=True,usetheta=True,usethetal=True):

        ###########################
        # Init new case structure
        ###########################

        newcase = Case(self.id,
                startDate=self.start_date, endDate=self.end_date)
        for att in self.attlist:
            newcase.set_attribute(att,self.attributes[att])

        ###########################
        # Interpolation
        ###########################

        dataout = {}

        if time is None:
            print 'WARNING: No time interpolation'
            for var in self.var_init_list + self.var_forcing_list:
                VV = self.variables[var]
                dataout[var] = Variable(var, data=VV.data, name=VV.name, units=VV.units,
                        height=VV.height, pressure=VV.pressure,
                        level=VV.level, time=VV.time,
                        plotcoef=VV.plotcoef, plotunits=VV.plotunits)
                if VV.time is not self.t0Axis:
                    dataout[var].time.id = 'time'
        else:
            timeout = Axis('time',time,name='forcing_time',units=self.tunits, calendar='gregorian')
            for var in self.var_init_list + self.var_forcing_list:
                VV = self.variables[var]
                if VV.time is not self.t0Axis:
                    dataout[var] = VV.interpol_time(time=timeout)
                else:
                    dataout[var] = Variable(var, data=VV.data, name=VV.name, units=VV.units,
                        height=VV.height, pressure=VV.pressure,
                        level=VV.level, time=VV.time,
                        plotcoef=VV.plotcoef, plotunits=VV.plotunits)

        if lev is None:

            print 'WARNING: No vertical interpolation'

            for var in self.var_init_list:
                VV = dataout[var]

                if VV.level is not None:
                    if VV.level.units == 'm':
                        levout = Axis('lev',VV.level.data,name='height',units='m')
                    elif VV.level.units == 'Pa':
                        levout = Axis('lev',VV.level.data,name='air_pressure',units='Pa')
                    else:
                        print 'ERROR: Level type undefined for level units:', VV.level.units
                        raise ValueError

                    height = None
                    pressure = None

                    if VV.height is not None:
                        height = Variable('zh', data=VV.height.data, units=VV.height.units, name='height',
                                level=levout, time=self.t0Axis)
                    if VV.pressure is not None:
                        pressure = Variable('pa', data=VV.pressure.data, units=VV.pressure.units, name='air_pressure',
                                level=levout, time=self.t0Axis)

                        
                    dataout[var] = Variable(var, data=VV.data, name=VV.name, units=VV.units,
                            height=height, pressure=pressure,
                            level=VV.level, time=self.t0Axis,
                            plotcoef=VV.plotcoef, plotunits=VV.plotunits)
                else:
                    dataout[var] = Variable(var, data=VV.data, name=VV.name, units=VV.units,
                            time=self.t0Axis,
                            plotcoef=VV.plotcoef, plotunits=VV.plotunits)

            for var in self.var_forcing_list:
                VV = dataout[var]

                if VV.level is not None:
                    if VV.level.units == 'm':
                        levout = Axis('lev',VV.level.data,name='height',units='m')
                    elif VV.level.units == 'Pa':
                        levout = Axis('lev',VV.level.data,name='air_pressure',units='Pa')
                    else:
                        print 'ERROR: Level type undefined for level units:', VV.level.units
                        raise ValueError

                    height = None
                    pressure = None

                    if VV.height is not None:
                        height = Variable('zh_forc', data=VV.height.data, units=VV.height.units, name='height_forcing',
                                level=levout, time=VV.time)
                    if VV.pressure is not None:
                        pressure = Variable('pa_forc', data=VV.pressure.data, units=VV.pressure.units, name='air_pressure_forcing',
                                level=levout, time=VV.time)

                        
                    dataout[var] = Variable(var, data=VV.data, name=VV.name, units=VV.units,
                            height=height, pressure=pressure,
                            level=VV.level, time=VV.time,
                            plotcoef=VV.plotcoef, plotunits=VV.plotunits)
                else:
                    dataout[var] = Variable(var, data=VV.data, name=VV.name, units=VV.units,
                            time=VV.time,
                            plotcoef=VV.plotcoef, plotunits=VV.plotunits)

        else:

            if levtype == 'altitude':

                levout = Axis('lev', lev, name='height', units='m')
                for var in self.var_init_list:
                    VV = dataout[var]
                    dataout[var] = VV.interpol_vert(height=lev)
                    if dataout[var].level is not None:
                        dataout[var].set_level(lev=levout)
                        dataout[var].height.set_level(lev=levout)
                        dataout[var].height.id = 'zh'
                        dataout[var].height.name = 'height'
                        dataout[var].set_coordinates('time','zh','lat','lon')
                        dataout[var].height.set_coordinates('time','zh','lat','lon')
                    #dataout[var].info()

                for var in self.var_forcing_list:
                    VV = dataout[var]
                    dataout[var] = VV.interpol_vert(height=lev)
                    if dataout[var].level is not None:
                        dataout[var].set_level(lev=levout)
                        dataout[var].height.set_level(lev=levout)
                        dataout[var].height.id = 'zh_forc'
                        dataout[var].height.name = 'height_forcing'
                        dataout[var].set_coordinates('time','zh_forc','lat','lon')
                        dataout[var].height.set_coordinates('time','zh_forc','lat','lon')
                    #dataout[var].info()

            elif levtype == 'pressure':

                levout = Axis('lev', lev, name='air_pressure', units='Pa')
                print 'ERROR: Pressure level type is not coded yet for interpolation'
                raise ValueError

            else:

                print 'ERROR: levtype unexpected:', levtype
                raise ValueError

        for var in self.var_init_list:
            newcase.var_init_list.append(var)
            newcase.variables[var] = dataout[var]

        for var in self.var_forcing_list:
            newcase.var_forcing_list.append(var)
            newcase.variables[var] = dataout[var]

        return newcase

    def add_missing_init_variables(self):

        _, nlev = self.variables['ua'].data.shape

        ##################################
        # Special case of height/pressure
        ##################################

        ps = self.variables['ps'].data[0]

        var = None
        for v in ['ta','theta','thetal']:
            if v in self.var_init_list:
                var = v
        if var is None:
            print 'ERROR: ta, theta or thetal must be given'
            raise ValueError

        VV = self.variables[var]
        levAxis = VV.level

        if VV.height is None and VV.pressure is None:

            print 'ERROR: height and pressure are None for {0}'.format(var)
            raise ValueError

        elif VV.pressure is None:

            height = VV.height

            kwargs = {}
            kwargs[var] = VV.data[0,:]
            kwargs['z'] = VV.height.data[0,:]
            kwargs['ps'] = ps
            if 'qv' in self.var_init_list:
                kwargs['qv'] = self.variables['qv'].data[0,:]
            elif 'qt' in self.var_init_list:
                kwargs['qv'] = self.variables['qt'].data[0,:]
            elif 'rv' in self.var_init_list:
                kwargs['qv'] = thermo.rt2qt(self.variables['rv'].data[0,:])
            elif 'rt' in self.var_init_list:
                kwargs['qv'] = thermo.rt2qt(self.variables['rt'].data[0,:])
            else:
                kwargs['qv'] = None

            pressure = thermo.z2p(**kwargs)
            pressure = np.reshape(pressure,(1,nlev))
            pressure = Variable('pa', data=pressure, name='air_pressure', units='Pa',
                        height=VV.height, 
                        level=VV.level, time=VV.time,
                        plotcoef=var_attributes['pa']['plotcoef'],
                        plotunits=var_attributes['pa']['plotunits'])

            self.variables['pa'] = pressure
            self.var_init_list.append('pa')

            for var in self.var_init_list:
                if self.variables[var].time.id == 't0':
                    self.variables[var].pressure = pressure

            if 'zh' not in self.var_init_list:
                self.var_init_list.append('zh')
                self.variables['zh'] = height

        elif VV.height is None:

            pressure = VV.pressure

            if 'zh' in self.var_init_list:
                height = self.variables['zh'].data
            else:
                kwargs = {}
                kwargs[var] = VV.data[0,:]
                kwargs['p'] = pressure.data[0,:]
                if 'qv' in self.var_init_list:
                    kwargs['qv'] = self.variables['qv'].data[0,:]
                elif 'qt' in self.var_init_list:
                    kwargs['qv'] = self.variables['qt'].data[0,:]
                elif 'rv' in self.var_init_list:
                    kwargs['qv'] = thermo.rt2qt(self.variables['rv'].data[0,:])
                elif 'rt' in self.var_init_list:
                    kwargs['qv'] = thermo.rt2qt(self.variables['rt'].data[0,:])
                else:
                    kwargs['qv'] = None

                height = thermo.p2z(**kwargs)

            height = np.reshape(height,(1,nlev))
            height = Variable('zh', data=height, name='height', units='Pa',
                        pressure=VV.pressure, 
                        level=VV.level, time=VV.time,
                        plotcoef=var_attributes['zh']['plotcoef'],
                        plotunits=var_attributes['zh']['plotunits'])

            self.variables['zh'] = height
            self.var_init_list.append('zh')

            for var in self.var_init_list:
                if self.variables[var].time.id == 't0':
                    self.variables[var].height = height

            if 'pa' not in self.var_init_list:
                self.variables['pa'] = pressure
                self.var_init_list.append('pa')

        else:
            
            if lwarnings: print 'WARNING: Nothing to do. height and pressure already defined. Just pass'

        ##################################
        # Initial state variables
        ##################################

        for var in ['zh','pa','ua','va','ta','theta','thetal','qv','qt','rv','rt','rl','ri','ql','qi','tke']:
            if var in self.var_init_list:
                pass 
            elif var == 'theta':
                theta = self.compute_theta(pressure=pressure)
                self.add_init_variable(var, theta, lev=levAxis, height=height, pressure=pressure)
            elif var == 'thetal':
                thetal = self.compute_thetal(pressure=pressure)
                self.add_init_variable(var, thetal, lev=levAxis, height=height, pressure=pressure)
            elif var == 'ta':
                temp = self.compute_temp(pressure=pressure)
                self.add_init_variable(var, temp, lev=levAxis, height=height, pressure=pressure)
            elif var == 'qv':
                qv = self.compute_qv()
                self.add_init_variable(var, qv, lev=levAxis, height=height, pressure=pressure)
            elif var == 'qt':
                qt = self.compute_qt()
                self.add_init_variable(var, qt, lev=levAxis, height=height, pressure=pressure)
            elif var == 'rv':
                rv = self.compute_rv()
                self.add_init_variable(var, rv, lev=levAxis, height=height, pressure=pressure)
            elif var == 'rt':
                rt = self.compute_rt()
                self.add_init_variable(var, rt, lev=levAxis, height=height, pressure=pressure)
            elif var in ['rl','ri','ql','qi','tke']:
                self.add_init_variable(var, self.variables['ua'].data*0, lev=levAxis, height=height, pressure=pressure)
            else:
                print 'ERROR: Case unexpected: variable {0} have to be defined'.format(var)
                sys.exit()

    def add_missing_forcing_variables(self):

        if not(self.var_forcing_list):
            print 'WARNING: no forcing variable. Nothing to do'
            # TODO: add ps_forc
            return

        ##################################
        # Special case of height/pressure
        ##################################

        height = None
        pressure = None

        for var in self.var_forcing_list:
            VV = self.variables[var]
            if VV.time is not None:
                time = VV.time
            if VV.level is not None:
                level = VV.level
            if VV.height is not None:
                height = VV.height
                nt, nlev = VV.data.shape
            if VV.pressure is not None:
                pressure = VV.pressure
                nt, nlev = VV.data.shape

        #---- Surface pressure forcing
        var = 'ps_forc'
        if var not in self.var_forcing_list:
            tmp = np.zeros((nt,),dtype=np.float64) + self.variables['ps'].data[0]
            self.add_variable(var, tmp, time=time)

        #---- Height/pressure
        if height is None and pressure is None:

            print 'ERROR: height and pressure are None. Unexpected in add_missing_forcing_variables'
            raise ValueError

        elif pressure is None:

            print 'assume pa_forc is constant over time'
            pressure = np.tile(self.variables['pa'].data[0,:],(nt,1))
            pressure = Variable('pa_forc', data=pressure, name='air_pressure_forcing', units='Pa',
                        height=height, 
                        level=level, time=time,
                        plotcoef=var_attributes['pa_forc']['plotcoef'],
                        plotunits=var_attributes['pa_forc']['plotunits'])

            self.variables['pa_forc'] = pressure
            self.var_forcing_list.append('pa_forc')

            for var in self.var_forcing_list:
                if self.variables[var].time.id == 'time':
                    self.variables[var].pressure = pressure

            if 'zh_forc' not in self.var_forcing_list:
                self.var_forcing_list.append('zh_forc')
                self.variables['zh_forc'] = height

        elif VV.height is None:

            print 'assume zh_forc is constant over time'
            height = np.tile(self.variables['zh'].data[0,:],(nt,1))
            height = Variable('zh_forc', data=height, name='height_forcing', units='m',
                        pressure=pressure, 
                        level=level, time=time,
                        plotcoef=var_attributes['zh_forc']['plotcoef'],
                        plotunits=var_attributes['zh_forc']['plotunits'])

            self.variables['zh_forc'] = height
            self.var_forcing_list.append('zh_forc')

            for var in self.var_forcing_list:
                if self.variables[var].time.id == 'time':
                    self.variables[var].height = height

            if 'pa_forc' not in self.var_forcing_list:
                self.var_forcing_list.append('pa_forc')
                self.variables['pa_forc'] = pressure

        else:
            
            if lwarnings: print 'WARNING: Nothing to do. height and pressure already defined. Just pass' 

        ##################################
        # Forcing variables
        ##################################

        #---- Large-scale temperature advection
        atts = ['adv_ta','adv_theta','adv_thetal']
        flag = False
        for att in atts:
            if att in self.attlist and self.attributes[att] == 1:
                flag = True

        if flag:
            # large-scale advection of temperature is active. All temperature variables are added, if needed.
            if 'tnta_adv' not in self.var_forcing_list:
                tadv = self.compute_tnta_adv()
                self.add_variable('tnta_adv', tadv, 
                    lev=level, time=time,
                    height=height, pressure=pressure)
                self.set_attribute('adv_ta',1)
            if 'tntheta_adv' not in self.var_forcing_list:
                thadv = self.compute_tntheta_adv()
                self.add_variable('tntheta_adv', thadv, 
                    lev=level, time=time,
                    height=height, pressure=pressure)
                self.set_attribute('adv_theta',1)
            if 'tnthetal_adv' not in self.var_forcing_list:
                thladv = self.compute_tnthetal_adv()
                self.add_variable('tnthetal_adv', thladv, 
                    lev=level, time=time,
                    height=height, pressure=pressure)
                self.set_attribute('adv_thetal',1)


        #---- Temperature radiative tendency
        att = 'radiation'
        if att in self.attlist and self.attributes[att] in ['tend',]:
            # radiative tendency is active. All temperature variables are added, if needed.

            if 'tnta_rad' not in self.var_forcing_list:
                trad = self.compute_tnta_rad()
                self.add_variable('tnta_rad', trad, 
                    lev=level, time=time,
                    height=height, pressure=pressure)
            if 'tntheta_rad' not in self.var_forcing_list:
                thrad = self.compute_tntheta_rad()
                self.add_variable('tntheta_rad', thrad, 
                    lev=level, time=time,
                    height=height, pressure=pressure)
            if 'tnthetal_rad' not in self.var_forcing_list:
                thlrad = self.compute_tnthetal_rad()
                self.add_variable('tnthetal_rad', thlrad, 
                    lev=level, time=time,
                    height=height, pressure=pressure)

        #---- Large-scale humidity advection
        atts = ['adv_qv','adv_qt','adv_rv','adv_rt']
        flag = False
        for att in atts:
            if att in self.attlist and self.attributes[att] == 1:
                flag = True 

        if flag: 
            # large-scale advection of humidity is active. All humidity variables are added, if needed.
            if 'tnqv_adv' not in self.var_forcing_list:
                qvadv = self.compute_tnqv_adv()
                self.add_variable('tnqv_adv', qvadv, 
                    lev=level, time=time,
                    height=height, pressure=pressure)
                self.set_attribute('adv_qv',1)
            if 'tnqt_adv' not in self.var_forcing_list:
                qtadv = self.compute_tnqt_adv()
                self.add_variable('tnqt_adv', qtadv, 
                    lev=level, time=time,
                    height=height, pressure=pressure)
                self.set_attribute('adv_qt',1)
            if 'tnrv_adv' not in self.var_forcing_list:
                rvadv = self.compute_tnrv_adv()
                self.add_variable('tnrv_adv', rvadv, 
                    lev=level, time=time,
                    height=height, pressure=pressure)
                self.set_attribute('adv_rv',1)
            if 'tnrt_adv' not in self.var_forcing_list:
                rtadv = self.compute_tnrt_adv()
                self.add_variable('tnrt_adv', rtadv, 
                    lev=level, time=time,
                    height=height, pressure=pressure)
                self.set_attribute('adv_rt',1)

        #---- Wind nudging

        for var in ['ua','va']:
            att = 'nudging_{0}'.format(var)
            if att in self.attlist and self.attributes[att] > 0:
                height = np.squeeze(self.variables['zh'].data)
                pressure = np.squeeze(self.variables['pa'].data)
                # Add further description of nudging altitude/pressure
                if 'zh_nudging_{0}'.format(var) not in self.attributes.keys():
                    zlev = thermo.plev2zlev(self.attributes['pa_nudging_{0}'.format(var)],height,pressure)
                    self.set_attribute('zh_nudging_{0}'.format(var),zlev)
                if 'pa_nudging_{0}'.format(var) not in self.attributes.keys():
                    plev = thermo.zlev2plev(self.attributes['zh_nudging_{0}'.format(var)],height,pressure)
                    self.set_attribute('pa_nudging_{0}'.format(var),plev)

        #---- Temperature nudging

        atts = ['nudging_ta','nudging_theta','nudging_thetal']
        flag = False
        zlev = None
        plev = None
        for att in atts:
            if att in self.attlist and self.attributes[att] != 0:
                flag = True
                nudging_timescale = self.attributes[att]
                if self.attributes[att] > 0: # simple nudging profile described in global attributes
                    if 'zh_{0}'.format(att) in self.attributes:
                        zlev = self.attributes['zh_{0}'.format(att)]
                    if 'pa_{0}'.format(att) in self.attributes:
                        plev = self.attributes['pa_{0}'.format(att)]

        if flag:
            if zlev is not None or plev is not None: # simple nudging profile described in global attributes
                # update nudging height/pressure for all temperature variables
                height = np.squeeze(self.variables['zh'].data)
                pressure = np.squeeze(self.variables['pa'].data)
                if zlev is None:
                    zlev = int(thermo.plev2zlev(plev,height,pressure))
                if plev is None:
                    plev = int(thermo.zlev2plev(zlev,height,pressure))

                for var in ['ta','theta','thetal']:
                    self.set_attribute('zh_nudging_{0}'.format(var),zlev)
                    self.set_attribute('pa_nudging_{0}'.format(var),plev)
            else:
                print 'ERROR: nudging case not yet coded'
                raise ValueError

            # temperature nudging is active. All temperature variables are added, if needed.
            if 'ta_nud' not in self.var_forcing_list:
                tnud = self.compute_ta_nud()
                self.add_variable('ta_nud', tnud, 
                    lev=level, time=time,
                    height=height, pressure=pressure)
                self.set_attribute('nudging_ta',nudging_timescale)
            if 'theta_nud' not in self.var_forcing_list:
                thnud = self.compute_theta_nud()
                self.add_variable('theta_nud', thnud, 
                    lev=level, time=time,
                    height=height, pressure=pressure)
                self.set_attribute('nudging_theta',nudging_timescale)
            if 'thetal_nud' not in self.var_forcing_list:
                thlnud = self.compute_thetal_nud()
                self.add_variable('thetal_nud', thlnud, 
                    lev=level, time=time,
                    height=height, pressure=pressure)
                self.set_attribute('nudging_thetal',nudging_timescale)

        #---- Humidity nudging

        atts = ['nudging_qv','nudging_qt','nudging_rv','nudging_rt']
        flag = False
        zlev = None
        plev = None
        for att in atts:
            if att in self.attlist and self.attributes[att] != 0:
                flag = True
                nudging_timescale = self.attributes[att]
                if self.attributes[att] > 0: # simple nudging profile described in global attributes
                    if 'zh_{0}'.format(att) in self.attributes:
                        zlev = self.attributes['zh_{0}'.format(att)]
                    if 'pa_{0}'.format(att) in self.attributes:
                        plev = self.attributes['pa_{0}'.format(att)]

        if flag:
            if zlev is not None or plev is not None: # simple nudging profile described in global attributes
                # update nudging height/pressure for all humidity variables
                height = np.squeeze(self.variables['zh'].data)
                pressure = np.squeeze(self.variables['pa'].data)
                if zlev is None:
                    zlev = int(thermo.plev2zlev(plev,height,pressure))
                if plev is None:
                    plev = int(thermo.zlev2plev(zlev,height,pressure))

                for var in ['qv','qt','rv','rt']:
                    self.set_attribute('zh_nudging_{0}'.format(var),zlev)
                    self.set_attribute('pa_nudging_{0}'.format(var),plev)
            else:
                print 'ERROR: nudging case not yet coded'
                raise ValueError

            # humidity nudging is active. All temperature variables are added, if needed.
            if 'qv_nud' not in self.var_forcing_list:
                qvnud = self.compute_qv_nud()
                self.add_variable('qv_nud', qvnud, 
                    lev=level, time=time,
                    height=height, pressure=pressure)
                self.set_attribute('nudging_qv',nudging_timescale)
            if 'qt_nud' not in self.var_forcing_list:
                qtnud = self.compute_qt_nud()
                self.add_variable('qt_nud', qtnud, 
                    lev=level, time=time,
                    height=height, pressure=pressure)
                self.set_attribute('nudging_qt',nudging_timescale)
            if 'rv_nud' not in self.var_forcing_list:
                rvnud = self.compute_rv_nud()
                self.add_variable('rv_nud', rvnud, 
                    lev=level, time=time,
                    height=height, pressure=pressure)
                self.set_attribute('nudging_rv',nudging_timescale)
            if 'rt_nud' not in self.var_forcing_list:
                rtnud = self.compute_rt_nud()
                self.add_variable('rt_nud', rtnud, 
                    lev=level, time=time,
                    height=height, pressure=pressure)
                self.set_attribute('nudging_rt',nudging_timescale)


    def convert2SCM(self, time=None, lev=None, levtype=None,
            usetemp=True, usetheta=True, usethetal=True):

        ###########################################
        # Interpolation
        ###########################################

        print '#'*40
        print '#### Interpolate available variables'

        caseSCM = self.interpolate(time=time, lev=lev, levtype=levtype,
                usetemp=usetemp, usetheta=usetheta, usethetal=usethetal)

        ###########################################
        # Add missing variables for initial state
        ###########################################

        print '#'*40
        print '#### Add missing initial variables'

        caseSCM.add_missing_init_variables()

        ###########################################
        # Add missing forcing variables
        ###########################################

        print '#'*40
        print '#### Add missing forcing variables'

        caseSCM.add_missing_forcing_variables()

        ###########################################
        # Final
        ###########################################

        print '#'*40

        return caseSCM


    def convert2SCM_old(self,time=None,lev=None,levtype=None,usetemp=True,usetheta=True,usethetal=True):

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
                    if levtype == 'altitude' and var in ['pressure','pressure_forc']:
                        dataout[var] = interpol(self.variables[var],levout=levout,timeout=timeout,log=True)
                    else:
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
            print 'assume height_forc is constant over time'
            tmp = np.zeros((ntout,nlev,1,1),dtype=np.float32)
            for it in range(0,ntout):
                tmp[it,:,0,0] = lev[:]
            caseSCM.add_variable(var,tmp,lev=lev,levtype=levtype,levid='lev',time=time,timeid='time')

        #---- Forcing presssure
        var = 'pressure_forc'
        if var in self.varlist:
            caseSCM.add_variable(var,dataout[var].data,lev=lev,levtype=levtype,levid='lev',time=time,timeid='time')
        else:
            print 'assume pressure_forc is constant over time'
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

        ltemp = False
        ltheta = False
        lthetal = False

        var = 'temp_nudging'
        att = 'nudging_temp'
        if att in self.attlist and self.attributes[att] > 0:
            if not(usetemp) and (ltheta or lthetal):
                print 'Warning: Several nudging variable for temperature are given, which might yield to inconsistencies'
                #sys.exit()
            else:
                ltemp=True
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
            if not(usetheta) and (ltemp or lthetal):
                print 'Warning: Several nudging variable for temperature are given, which might yield to inconsistencies'
                #sys.exit()
            else:
                ltheta=True   
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
            if not(usethetal) and (ltemp or ltheta):
                print 'Warning: Several nudging variable for temperature are given, which might yield to inconsistencies'
                #sys.exit()
            else:
                lthetal=True   
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
            if 'ts' in dataout.keys():
                caseSCM.add_variable('ts', dataout['ts'].data, time=time,timeid='time')

        if att in self.attlist and self.attributes[att] == 'ts':
            caseSCM.add_variable('ts',dataout['ts'].data,time=time,timeid='time')

        att = 'surfaceForcingWind'
        if att in self.attlist and self.attributes[att] == 'ustar':
            caseSCM.add_variable('ustar',dataout['ustar'].data,time=time,timeid='time')

        ###########################
        # Final
        ###########################

        return caseSCM

