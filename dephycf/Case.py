#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 27 November 2019

@author: Romain Roehrig
"""

import os
import logging
logging.basicConfig(format='%(asctime)s - %(name)20s - %(levelname)s - %(message)s', level=logging.INFO)
logger = logging.getLogger(__name__)

import time
import copy

from datetime import datetime 

import netCDF4 as nc

import numpy as np

from .Axis import Axis
from .Variable import Variable, read as readvar, interpol

from .variables_attributes import attributes as var_attributes
from .attributes import known_attributes, required_attributes

from . import thermo
from . import constants as CC

# Default start and en dates
startDate0 = datetime(1979,1,1,0,0,0)
endDate0 = datetime(1979,1,1,0,0,0)

init_vars_1D = ['ps','ts','thetas']
forc_vars_1D = ['ps_forc','hfss','hfls','ustar',\
                'ts_forc','thetas_forc','tskin',\
                'orog','lat','lon','z0','z0h','z0q','beta','alb','emis']

class Case:

    def __init__(self, caseid,
            lat=None, lon=None,
            startDate=startDate0, endDate=endDate0,
            surfaceType='ocean', zorog=0.,
            forcing_scale=-1):

        self.id = caseid

        self.set_dates(startDate,endDate)

        # Latitude (degrees_noth) and Longitude (degrees_east)
        self.lat = lat
        self.lon = lon

        # Surface type
        self.surface_type = surfaceType

        # Forcing scale
        self.forcing_scale = forcing_scale

        # Variables
        self.var_init_list = []
        self.var_forcing_list = []
        self.variables = {}

        # Attributes
        self.attlist = ['case','title','reference','author','version','format_version','modifications','script','comment',
                'start_date','end_date',
                'forcing_scale',
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
                'forcing_scale': self.forcing_scale,
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
            logger.warning("Warning, attribute {0} is not known. It might not be written in the case output file.".format(attid))

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

###################################################################################################
#                  Generic removal/addition of a variable
###################################################################################################

    def remove_variable(self, varid):
        """Remove a variable of a Case object

        Require arguments:
        varid -- string for the variable id
        """

        if varid in self.var_init_list:
            self.var_init_list.remove(varid)

        if varid in self.var_forcing_list:
            self.var_forcing_list.remove(varid)

        del(self.variables[varid])


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
            logger.debug('Variable {0} is already defined. It will be overwritten'.format(varid))
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
            raise ValueError('time axis must not be None for variable {0}'.format(varid))
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
                logger.error('levtype unexpected: {0}'.format(levtype))
                logger.error('levtype should be defined for variable {0} and given as altitude or pressure'.format(varid))
                raise ValueError('levtype unexpected: {0}'.format(levtype))

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
                    self.set_attribute('forc_z',1)
            elif levtype == 'pressure':
                if pressure is None:
                    pressure = np.tile(levdata,(nt,1))
                    pressure_id = 'pa_{0}'.format(varid)
                    pressure_units = 'Pa'
                    self.set_attribute('forc_p',1)

        ######################
        # Get variable attributes
        if name is None:
            varname = var_attributes[varid]['name']
        else:
            varname = name

        if units is None:
            varunits = var_attributes[varid]['units']
        else:
            logger.warning('Warning: the framework expects SI units')
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
            logger.error('level and time axes are None. Case unexpected')
            raise ValueError('level and time axes are None. Case unexpected')
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

###################################################################################################
#                  Setting initial state
###################################################################################################

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

        if varid in init_vars_1D:
            # Put the expected shape of the input data
            tmp = np.reshape(vardata,(1,))
        else:
            # Check if lev optional argument is given
            if 'lev' not in kwargs:
                logger.error('level axis should be given for variable {0}'.format(varid))
                raise ValueError('level axis should be given for variable {0}'.format(varid))

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

    def add_init_ts(self,vardata,**kwargs):
        """Add initial state variable for surface temperature to a Case object.
           
        Required argument:
        vardata -- input data as an integer or a float.

        See add_variable function for optional arguments.
        """

        self.add_init_variable('ts',vardata,**kwargs)

    def add_init_thetas(self,vardata,**kwargs):
        """Add initial state variable for surface potential temperature to a Case object.
           
        Required argument:
        vardata -- input data as an integer or a float.

        See add_variable function for optional arguments.
        """

        self.add_init_variable('thetas',vardata,**kwargs)

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
        self.set_attribute('ini_ta',1)
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
        self.set_attribute('ini_theta',1)
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
        self.set_attribute('ini_thetal',1)
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
        self.set_attribute('ini_qv',1)
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
        self.set_attribute('ini_qt',1)
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
        self.set_attribute('ini_rv',1)
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
        self.set_attribute('ini_rt',1)
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
            logger.error('You must provide both zonal and meridional wind')
            raise ValueError('You must provide both zonal and meridional wind')

        if 'lev' not in kwargs and ulev is None and vlev is None:
            logger.error('You must provide a vertical axis either with lev or with both ulev/vlev')
            raise ValueError('You must provide a vertical axis either with lev or with both ulev/vlev')

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

###################################################################################################
#                  Setting forcing
###################################################################################################

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
        if 'time' in kwargs and kwargs['time'] is not None:
            lconstant = False
            nt, = np.array(kwargs['time']).shape
        else: # forcing is constant in time
            lconstant = True
            kwargs['time'] = [self.tstart,self.tend]
            nt = 2

        if varid in forc_vars_1D:
            # Put the expected shape of the input data
            if lconstant: 
                tmp = np.zeros((nt),dtype=np.float32)
                tmp[0] = vardata
                tmp[1] = vardata
            else:
                tmp = np.reshape(vardata,(nt,))
        else:
            # Check if lev optional argument is given
            if 'lev' not in kwargs:
                logger.error('level axis should be given for variable {0}'.format(varid))
                raise ValueError('level axis should be given for variable {0}'.format(varid))

            if len(np.array(kwargs['lev']).shape) != 1:
                logger.error('level axis should have only one dimension for variable {0}'.format(varid))
                logger.error('You may have mixed lev and height/pressure arguments')
                logger.error('lev argument is used to provide an axis within the netcdf files')

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
            logger.error('You must provide both zonal and meridional geostrophic wind')
            raise ValueError('You must provide both zonal and meridional geostrophic wind')

        if 'lev' not in kwargs and uglev is None and vglev is None:
            logger.error('You must provide a vertical axis either with lev or with both uglev/vglev')
            raise ValueError('You must provide a vertical axis either with lev or with both uglev/vglev')

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
            logger.error('You must provide either w or omega')
            raise ValueError('You must provide either w or omega')

        if w is not None and omega is not None:
            logger.error('You cannot provide both w and omega at the same time')
            raise ValueError('You cannot provide both w and omega at the same time')

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

    def add_temp_radiation_tendency(self,data,**kwargs):
        """Add a temperature radiative tendency to a Case object.
        
        Required argument:
        data -- input data as a list or a numpy array.

        See add_variable function for optional arguments.
        Note that:
        - a level axis is required (lev optional argument).
        - a levtype is required (levtype optional argument).

        If time is not provided, forcing is assumed constant in time.
        """

        self.set_attribute('radiation','tend')
        self.add_forcing_variable('tnta_rad',data,**kwargs)

    def add_theta_radiation_tendency(self,data,**kwargs):
        """Add a potential temperature radiative tendency to a Case object.
        
        Required argument:
        data -- input data as a list or a numpy array.

        See add_variable function for optional arguments.
        Note that:
        - a level axis is required (lev optional argument).
        - a levtype is required (levtype optional argument).

        If time is not provided, forcing is assumed constant in time.
        """

        self.set_attribute('radiation','tend')
        self.add_forcing_variable('tntheta_rad',data,**kwargs)

    def add_thetal_radiation_tendency(self,data,**kwargs):
        """Add a liquid-water potential temperature radiative tendency to a Case object.
        
        Required argument:
        data -- input data as a list or a numpy array.

        See add_variable function for optional arguments.
        Note that:
        - a level axis is required (lev optional argument).
        - a levtype is required (levtype optional argument).

        If time is not provided, forcing is assumed constant in time.
        """

        self.set_attribute('radiation','tend')
        self.add_forcing_variable('tnthetal_rad',data,**kwargs)

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

    def add_nudging(self,varid,data,timescale=None,z_nudging=None,p_nudging=None,lev=None,nudging_coefficient=None,lev_coef=None,**kwargs):
        """Add a nudging forcing to a Case object.
        
        Required argument:
        varid     -- id of the variable to be nudged as a string
        data      -- input data as a list or a numpy array.
        timescale -- nudging timescale in seconds (integer or float)
        z_nudging -- altitude above which nudging is applied (integer or float)
        p_nudging -- pressure altitude under which nudging is applied (integer or float)
        nudging_coefficient -- profile of nudging coefficien

        Either timescale or nudging_coefficient must be defined.
        If z_nudging and p_nudging are not provided, nudging is assumed to over the whole atmosphere

        See add_variable function for optional arguments.
        Note that:
        - a level axis is required (lev optional argument).
        - a levtype is required (levtype optional argument).

        If time is not provided, forcing is assumed constant in time.
        """

        if timescale is not None and nudging_coefficient is not None:
            logger.error('You cannot provide both a nudging timescale and a nudging coefficient profile for {0}'.format(varid))
            raise ValueError('You must provide a nudging timescale for variable {0}'.format(varid))

        if timescale is not None:
            if timescale >= 0.:
                self.set_attribute('nudging_{0}'.format(varid),float(timescale))
            else:
                logger.error('The nudging timescale for {0} is expected to be positive, but it is equal to {1}'.format(varid,timescale))
                raise ValueError('The nudging timescale for {0} is expected to be positive, but it is equal to {1}'.format(varid,timescale))

            if z_nudging is not None and p_nudging is not None:
                logger.warning('For {0}, both z_nudging and p_nudging is provided. Be sure they are consistent'.format(varid))

            if z_nudging is None and p_nudging is None:
                logger.warning('{0} will be nudged over the whole atmosphere'.format(varid))
                self.set_attribute('zh_nudging_{0}'.format(varid),0)

            if z_nudging is not None:
                self.set_attribute('zh_nudging_{0}'.format(varid),float(z_nudging))
            if p_nudging is not None:
                self.set_attribute('pa_nudging_{0}'.format(varid),float(p_nudging))



        elif nudging_coefficient is not None:
            self.set_attribute('nudging_{0}'.format(varid),-1.)
            var= 'nudging_coefficient_{0}'.format(varid)
            if lev_coef is None:
                logger.warning('No vertical levels are provided for the nudging coefficient of variable {0}'.format(varid))
                logger.warning('It is assumed it is the same as variable {0}'.format(varid))
                lev_coef = kwargs['lev']
            self.add_forcing_variable(var,nudging_coefficient,lev=lev_coef,**kwargs)

        else:
            logger.error('You must provide a nudging timescale or a nudging coefficient profile for variable {0}'.format(varid))
            raise ValueError('You must provide a nudging timescale or a nudging coefficient profile for variable {0}'.format(varid))


        var = '{0}_nud'.format(varid)
        self.add_forcing_variable(var,data,lev=lev,**kwargs)

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
            logger.error('You must provide both zonal and meridional nudging wind')
            raise ValueError('You must provide both zonal and meridional nudging wind')

        if 'lev' not in kwargs and ulev is None and vlev is None:
            logger.error('You must provide a vertical axis either with lev or with both ulev/vlev')
            raise ValueError('You must provide a vertical axis either with lev or with both ulev/vlev')

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

        self.set_attribute('radiation',"on")

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

    def add_forcing_ts(self,data,z0=None,z0h=None,z0q=None,**kwargs):
        """Add a surface temperature forcing to a Case object.
        
        This function sets a surface temperature forcing as the case surface forcing.
        In case the initial surface temperature is not defined, add it.

        Required argument:
        data -- input data as a list or a numpy array.

        If time is not provided, forcing is assumed constant in time.
        """

        self.set_attribute('surface_forcing_temp','ts')
        self.add_forcing_variable('ts_forc',data,**kwargs)

        if 'ts' not in self.var_init_list:
            if isinstance(data,float):
                self.add_init_ts(data,**kwargs)
            else:
                self.add_init_ts(data[0],**kwargs)

        if z0 is not None:
            self.set_attribute('surface_forcing_wind','z0')
            self.add_forcing_variable('z0',z0)
            if z0h is not None:
                self.add_forcing_variable('z0h',z0h)
            if z0q is not None:
                self.add_forcing_variable('z0q',z0q)

    def add_forcing_thetas(self,data,z0=None,z0h=None,z0q=None,**kwargs):
        """Add a surface potential temperature forcing to a Case object.
        
        This function sets a surface temperature forcing as the case surface forcing.
        In case the initial surface temperature is not defined, add it.

        Required argument:
        data -- input data as a list or a numpy array.

        If time is not provided, forcing is assumed constant in time.
        """

        self.set_attribute('surface_forcing_temp','thetas')
        self.add_forcing_variable('thetas_forc',data,**kwargs)

        if 'thetas' not in self.var_init_list:
            if isinstance(data,float):
                self.add_init_thetas(data,**kwargs)
            else:
                self.add_init_thetas(data[0],**kwargs)

        if z0 is not None:
            self.set_attribute('surface_forcing_wind','z0')
            self.add_forcing_variable('z0',z0)
            if z0h is not None:
                self.add_forcing_variable('z0h',z0h)
            if z0q is not None:
                self.add_forcing_variable('z0q',z0q)

    def add_surface_fluxes(self,sens=None,lat=None,time_sens=None,time_lat=None,\
                           forc_wind=None,z0=None,time_z0=None,ustar=None,time_ustar=None,**kwargs):
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
            logger.error('You must provide both sensible and latent heat fluxes')
            raise ValueError('You must provide both sensible and latent heat fluxes')

        self.set_attribute('surface_forcing_temp','surface_flux')
        self.set_attribute('surface_forcing_moisture','surface_flux')

        if forc_wind == 'z0':
            self.set_attribute("surface_forcing_wind","z0")
            if z0 is None:
                logger.error('z0 must be provided')
                raise ValueError('z0 must be provided')
            self.add_forcing_variable('z0',z0,time=time_z0)
        elif forc_wind == 'ustar':
            self.set_attribute("surface_forcing_wind","ustar")
            if ustar is None:
                logger.error('ustar must be provided')
                raise ValueError('ustar must be provided')
            self.add_forcing_variable('ustar',ustar,time=time_ustar)
        elif forc_wind is None:
            logger.error('You must specify the surface wind forcing as z0 or ustar')
            raise ValueError('You must specify the surface wind forcing as z0 or ustar')
        else:
            logger.error('Surface wind forcing unexpected: {0}. It should be either z0 or ustar'.format(forc_wind))
            raise ValueError('Surface wind forcing unexpected: {0}. It should be either z0 or ustar'.format(forc_wind))

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

###################################################################################################
#                  Case information
###################################################################################################

    def info(self):
        """Print information about a case object.

        No argument required.
        """

        for att in known_attributes:
            if att in self.attlist:
                print('{0}: {1}'.format(att,self.attributes[att]))
     
        print("######################")
        print("# Variable information")
        for var in self.var_init_list + self.var_forcing_list:
            self.variables[var].info()

###################################################################################################
#                  Reading/Writing a Case object in a netCDF file
###################################################################################################

    def write(self,fileout):
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
                logger.debug('Writing time axis for {0}'.format(var))
                self.variables[var].write(g,
                        write_time_axes=True, write_level_axes=False,
                        write_data=False, write_vertical=False)

        # Writing first only level axes
        for var in var_attributes.keys():
            if var in self.var_init_list + self.var_forcing_list:
                logger.debug('Writing level axis for {0}'.format(var))
                self.variables[var].write(g,
                        write_time_axes=False, write_level_axes=True,
                        write_data=False, write_vertical=False)

        # Writing then only vertical variables
        for var in var_attributes.keys():
            if var in self.var_init_list + self.var_forcing_list:
                logger.debug('Writing vertical variable for {0}'.format(var))
                self.variables[var].write(g, 
                        write_time_axes=False, write_level_axes=False,
                        write_data=False, write_vertical=True)

        # Finally write data
        for var in var_attributes.keys():
            if var in self.var_init_list + self.var_forcing_list:
                logger.debug('Writing data for {0}'.format(var))
                self.variables[var].write(g,
                        write_time_axes=False, write_level_axes=False,
                        write_data=True, write_vertical=False)

        # Write global attributes
        for att in known_attributes:
            if att in self.attlist:
                g.setncattr(att,self.attributes[att])

        g.close()

    def read(self,filein):
        """Read a netCDF file to store data and information into a case object.

        Required argument:
        filein -- netCDF file name as a string

        Note that the case object should be initialize first.
        """

        f = nc.Dataset(filein,'r')

        self.set_dates(f.start_date,f.end_date)

        for var in f.variables:
            if (var not in f.dimensions and var[0:6] != 'bounds' and var[0:3] not in ['zh_','pa_'])\
                    or var in ['zh_forc','pa_forc']:
                logger.debug('Reading {0}'.format(var))
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
                if logger.level == logging.DEBUG:
                    tmp.info()

        for att in known_attributes:
            try:
                self.set_attribute(att,f.getncattr(att))
            except AttributeError:
                pass
            except:
                raise

        f.close()

###################################################################################################
#                  Basic plotting
###################################################################################################

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

###################################################################################################
#                  Computing additional variables
###################################################################################################

    def compute_theta(self,pressure=None):

        if 'ta' in self.var_init_list:
            if pressure is None:
                logger.error('Pressure should be None to theta from ta')
                raise ValueError('Pressure should be None to theta from ta')
            else:
                 logger.info('Compute potential temperature from pressure and temperature')
                 theta = thermo.t2theta(p=pressure.data[0,:], temp=self.variables['ta'].data[0,:])
        elif 'thetal' in self.var_init_list:
            logger.info('Assume theta=thetal')
            theta = self.variables['thetal'].data[0,:]
        else:
            logger.error('At least ta or thetal should be given to compute theta')
            raise ValueError('At least ta or thetal should be given to compute theta')

        return theta

    def compute_thetal(self,pressure=None):

        if 'ta' in self.var_init_list:
            if pressure is None:
                logger.error('Pressure should be None to compute thetal from ta')
                raise ValueError('Pressure should be None to compute thetal from ta')
            else:
                logger.info('Compute potential temperature from pressure and temperature, assuming thetal=theta')
                thetal = thermo.t2theta(p=pressure.data[0,:], temp=self.variables['ta'].data[0,:])
        elif 'theta' in self.var_init_list:
            logger.info('Assume thetal=theta')
            thetal = self.variables['theta'].data[0,:]
        else:
            logger.error('At least ta or theta should be given to compute thetal')
            raise ValueError('At least ta or theta should be given to compute thetal')

        return thetal

    def compute_temp(self,pressure=None):

        if pressure is None:
            logger.error('Pressure should be None to compute ta from theta or thetahl')
            raise ValueError('Pressure should be None to compute ta from theta or thetahl')

        if 'theta' in self.var_init_list:
            logger.info('Compute temperature from pressure and potential temperature')
            temp = thermo.theta2t(p=pressure.data[0,:], theta=self.variables['theta'].data[0,:])
        elif 'thetal' in self.var_init_list:
            logger.info('Compute temperature from pressure and liquid potential temperature (No liquid water considered)')
            temp = thermo.theta2t(p=pressure.data[0,:], theta=self.variables['thetal'].data[0,:])
        else:
            logger.error('At least theta or thetal should be given')
            raise ValueError('At least theta or thetal should be given')

        return temp

    def compute_qv(self):

        if 'qt' in self.var_init_list:
            logger.info('Assume qv=qt')
            qv = self.variables['qt'].data[0,:]
        elif 'rv' in self.var_init_list:
            logger.info('Compute qv from rv')
            qv = thermo.rt2qt(self.variables['rv'].data[0,:])
        elif 'rt' in self.var_init_list:
            logger.info('Compute qt from rt and assume qv=qt')
            qv = thermo.rt2qt(self.variables['rt'].data[0,:])
        else:
            logger.error('Either qt, rv or rt should be defined to compute qv')
            raise ValueError('Either qt, rv or rt should be defined to compute qv')

        return qv

    def compute_qt(self):

        if 'qv' in self.var_init_list:
            logger.info('Assume qt=qv')
            qt = self.variables['qv'].data[0,:]
        elif 'rv' in self.var_init_list:
            logger.info('Compute qv from rv and assume qt=qv')
            qt = thermo.rt2qt(self.variables['rv'].data[0,:])
        elif 'rt' in self.var_init_list:
            logger.info('Compute qt from rt')
            qt = thermo.rt2qt(self.variables['rt'].data[0,:])
        else:
            logger.error('Either qv, rv or rt should be defined to compute qt')
            raise ValueError('Either qv, rv or rt should be defined to compute qt')

        return qt

    def compute_rv(self):

        if 'qv' in self.var_init_list:
            logger.info('Compute rv from qv')
            rv = thermo.qt2rt(self.variables['qv'].data[0,:])
        elif 'qt' in self.var_init_list:
            logger.info('Compute rt from qt and assume rv=rt')
            rv = thermo.qt2rt(self.variables['qt'].data[0,:])
        elif 'rt' in self.var_init_list:
            logger.info('Assume rv=rt')
            rv = selv.variables['rt'].data[0,:]
        else:
            logger.error('Either qv, qt or rt should be defined to compute rv')
            raise ValueError('Either qv, qt or rt should be defined to compute rv')

        return rv

    def compute_rt(self):

        if 'qv' in self.var_init_list:
            logger.info('Compute rt from qt, assuming qt=qv')
            rt = thermo.qt2rt(self.variables['qv'].data[0,:])
        elif 'qt' in self.var_init_list:
            logger.info('Compute rt from qt')
            rt = thermo.qt2rt(self.variables['qt'].data[0,:])
        elif 'rv' in self.var_init_list:
            logger.info('Assume rt=rv')
            rt = selv.variables['rv'].data[0,:]
        else:
            logger.error('Either qv, qt or rv should be defined to compute rt')
            raise ValueError('Either qv, qt or rv should be defined to compute rt')

        return rt

    def compute_tnta_adv(self):

        pressure = self.variables['pa_forc'].data

        if 'tntheta_adv' in self.var_forcing_list:
            logger.info('Compute tnta_adv from tntheta_adv')
            thadv = self.variables['tntheta_adv'].data
            tadv = thermo.theta2t(p=pressure,theta=thadv)
        elif 'tnthetal_adv' in self.var_forcing_list:
            logger.info('Assume tnthetal_adv=tntheta_adv and compute tnta_adv from tntheta_adv')
            thladv = self.variables['tnthetal_adv'].data
            tadv = thermo.theta2t(p=pressure,theta=thladv)
        else:
            logger.error('To compute tnta_adv, tntheta_adv or tnthetal_adv must be known')
            raise ValueError('To compute tnta_adv, tntheta_adv or tnthetal_adv must be known')

        return tadv

    def compute_tntheta_adv(self):

        pressure = self.variables['pa_forc'].data

        if 'tnta_adv' in self.var_forcing_list:
            logger.info('Compute tntheta_adv from tnta_adv')
            tadv = self.variables['tnta_adv'].data
            thadv = thermo.t2theta(p=pressure,temp=tadv)
        elif 'tnthetal_adv' in self.var_forcing_list:
            logger.info('Assume tntheta_adv=tnthetal_adv')
            thadv = self.variables['tnthetal_adv'].data
        else:
            logger.error('To compute tntheta_adv, tnta_adv or tnthetal_adv must be known')
            raise ValueError('To compute tntheta_adv, tnta_adv or tnthetal_adv must be known')

        return thadv

    def compute_tnthetal_adv(self):

        pressure = self.variables['pa_forc'].data

        if 'tnta_adv' in self.var_forcing_list:
            logger.info('Compute tnthetal_adv from tnta_adv assuming tnthetal_adv=tntheta_adv')
            tadv = self.variables['tnta_adv'].data
            thladv = thermo.t2theta(p=pressure,temp=tadv)
        elif 'tntheta_adv' in self.var_forcing_list:
            logger.info('Assume tnthetal_adv=tntheta_adv')
            thladv = self.variables['tntheta_adv'].data
        else:
            logger.error('To compute tnthetal_adv, tnta_adv or tntheta_adv must be known')
            raise ValueError('To compute tnthetal_adv, tnta_adv or tntheta_adv must be known')

        return thladv

    def compute_tnta_rad(self):

        pressure = self.variables['pa_forc'].data

        if 'tntheta_rad' in self.var_forcing_list:
            logger.info('Compute tnta_rad from tntheta_rad')
            thrad = self.variables['tntheta_rad'].data
            trad = thermo.theta2t(p=pressure,theta=thrad)
        elif 'tnthetal_rad' in self.var_forcing_list:
            logger.info('Assume tnthetal_rad=tntheta_rad and compute tnta_rad from tntheta_rad')
            thlrad = self.variables['tnthetal_rad'].data
            trad = thermo.theta2t(p=pressure,theta=thlrad)
        else:
            logger.error('To compute tnta_rad, tntheta_rad or tnthetal_rad must be known')
            raise ValueError('To compute tnta_rad, tntheta_rad or tnthetal_rad must be known')

        return trad

    def compute_tntheta_rad(self):

        pressure = self.variables['pa_forc'].data

        if 'tnta_rad' in self.var_forcing_list:
            logger.info('Compute tntheta_rad from tnta_rad')
            trad = self.variables['tnta_rad'].data
            thrad = thermo.t2theta(p=pressure,temp=trad)
        elif 'tnthetal_rad' in self.var_forcing_list:
            logger.info('Assume tntheta_rad=tnthetal_rad')
            thrad = self.variables['tnthetal_rad'].data
        else:
            logger.error('To compute tntheta_rad, tnta_rad or tnthetal_rad must be known')
            raise ValueError('To compute tntheta_rad, tnta_rad or tnthetal_rad must be known')

        return thrad

    def compute_tnthetal_rad(self):

        pressure = self.variables['pa_forc'].data

        if 'tnta_rad' in self.var_forcing_list:
            logger.info('Compute tnthetal_rad from tnta_rad assuming tnthetal_rad=tntheta_rad')
            trad = self.variables['tnta_rad'].data
            thlrad = thermo.t2theta(p=pressure,temp=trad)
        elif 'tntheta_rad' in self.var_forcing_list:
            logger.info('Assume tnthetal_rad=tntheta_rad')
            thlrad = self.variables['tntheta_rad'].data
        else:
            logger.error('To compute tnthetal_rad, tnta_rad or tntheta_rad must be known')
            raise ValueError('To compute tnthetal_rad, tnta_rad or tntheta_rad must be known')

        return thlrad

    def compute_tnqv_adv(self):

        if 'tnqt_adv' in self.var_forcing_list:
            logger.info('Assume tnqv_adv=tnqt_adv')
            qvadv = self.variables['tnqt_adv'].data
        elif 'tnrv_adv' in self.var_forcing_list:
            logger.info('Compute tnqv_adv from tnrv_adv using initial rv profile')
            rvadv = self.variables['tnrv_adv'].data
            rv = self.variables['rv'].data
            qvadv =  thermo.advrt2advqt(rt=rv,advrt=rvadv)
        elif 'tnrt_adv' in self.var_forcing_list:
            logger.info('Compute tnqv_adv from tnrv_adv assuming tnrv_adv=tnrt_adv and using initial rt profile')
            rtadv = self.variables['tnrt_adv'].data
            rt = self.variables['rt'].data
            qvadv =  thermo.advrt2advqt(rt=rt,advrt=rtadv)
        else:
            logger.error('To compute tnqv_adv, tnqt_adv, tnrv_adv or tnrt_adv must be known')
            raise ValueError('To compute tnqv_adv, tnqt_adv, tnrv_adv or tnrt_adv must be known')

        return qvadv

    def compute_tnqt_adv(self):

        if 'tnqv_adv' in self.var_forcing_list:
            logger.info('Assume tnqt_adv=tnqv_adv')
            qtadv = self.variables['tnqv_adv'].data
        elif 'tnrv_adv' in self.var_forcing_list:
            logger.info('Compute tnqt_adv from tnrt_adv assuming tnrt_adv=tnrv_adv using initial rv profile')
            rtadv = self.variables['tnrv_adv'].data
            rt = self.variables['rv'].data
            qtadv =  thermo.advrt2advqt(rt=rt,advrt=rtadv)
        elif 'tnrt_adv' in self.var_forcing_list:
            logger.info('Compute tnqt_adv from tnrt_adv using initial rt profile')
            rtadv = self.variables['tnrt_adv'].data
            rt = self.variables['rt'].data
            qtadv =  thermo.advrt2advqt(rt=rt,advrt=rtadv)
        else:
            logger.error('To compute tnqt_adv, tnqv_adv, tnrv_adv or tnrt_adv must be known')
            raise ValueError('To compute tnqt_adv, tnqv_adv, tnrv_adv or tnrt_adv must be known')

        return qtadv

    def compute_tnrv_adv(self):

        if 'tnrt_adv' in self.var_forcing_list:
            logger.info('Assume tnrv_adv=tnrt_adv')
            rvadv = self.variables['tnrt_adv'].data
        elif 'tnqv_adv' in self.var_forcing_list:
            logger.info('Compute tnrv_adv from tnqv_adv using initial qv profile')
            qvadv = self.variables['tnqv_adv'].data
            qv = self.variables['qv'].data
            rvadv =  thermo.advqt2advrt(qt=qv,advqt=qvadv)
        elif 'tnqt_adv' in self.var_forcing_list:
            logger.info('Compute tnrv_adv from tnqv_adv assuming tnqv_adv=tnqt_adv and using initial qt profile')
            qvadv = self.variables['tnqt_adv'].data
            qv = self.variables['qt'].data
            rvadv =  thermo.advqt2advrt(qt=qv,advqt=qvadv)
        else:
            logger.error('To compute tnrv_adv, tnqv_adv, tnqt_adv or tnrt_adv must be known')
            raise ValueError('To compute tnrv_adv, tnqv_adv, tnqt_adv or tnrt_adv must be known')

        return rvadv

    def compute_tnrt_adv(self):

        if 'tnrv_adv' in self.var_forcing_list:
            logger.info('Assume tnrt_adv=tnrv_adv')
            rtadv = self.variables['tnrv_adv'].data
        elif 'tnqv_adv' in self.var_forcing_list:
            logger.info('Compute tnrt_adv from tnqt_adv assuming tnqt_adv=tnqv_adv and using initial qv profile')
            qtadv = self.variables['tnqv_adv'].data
            qt = self.variables['qv'].data
            rtadv =  thermo.advqt2advrt(qt=qt,advqt=qtadv)
        elif 'tnqt_adv' in self.var_forcing_list:
            logger.info('Compute tnrt_adv from tnqt_adv using initial qt profile')
            qtadv = self.variables['tnqt_adv'].data
            qt = self.variables['qt'].data
            rtadv =  thermo.advqt2advrt(qt=qt,advqt=qtadv)
        else:
            logger.error('To compute tnrt_adv, tnqv_adv, tnqt_adv or tnrv_adv must be known')
            raise ValueError

        return rtadv

    def compute_ta_nud(self):

        pressure = self.variables['pa_forc'].data

        if 'theta_nud' in self.var_forcing_list:
            logger.info('Compute ta_nud from theta_nud')
            thnud = self.variables['theta_nud'].data
            tnud = thermo.theta2t(p=pressure,theta=thnud)
        elif 'thetal_nud' in self.var_forcing_list:
            logger.info('Compute ta_nud from theta_nud assuming theta_nud=thetal_nud')
            thnud = self.variables['thetal_nud'].data
            tnud = thermo.theta2t(p=pressure,theta=thnud)
        else:
            logger.error('To compute ta_nud, theta_nud or thetal_nud must be known')
            raise ValueError('To compute ta_nud, theta_nud or thetal_nud must be known')

        return tnud

    def compute_theta_nud(self):

        pressure = self.variables['pa_forc'].data

        if 'ta_nud' in self.var_forcing_list:
            logger.info('Compute theta_nud from ta_nud')
            tnud = self.variables['ta_nud'].data
            thnud = thermo.t2theta(p=pressure,temp=tnud)
        elif 'thetal_nud' in self.var_forcing_list:
            logger.info('Assume theta_nud=thetal_nud')
            thnud = self.variables['thetal_nud'].data
        else:
            logger.error('To compute theta_nud, ta_nud or thetal_nud must be known')
            raise ValueError('To compute theta_nud, ta_nud or thetal_nud must be known')

        return thnud

    def compute_thetal_nud(self):

        pressure = self.variables['pa_forc'].data

        if 'ta_nud' in self.var_forcing_list:
            logger.info('Compute thetal_nud from ta_nud assuming thetal_nud=theta_nud')
            tnud = self.variables['ta_nud'].data
            thlnud = thermo.t2theta(p=pressure,temp=tnud)
        elif 'theta_nud' in self.var_forcing_list:
            logger.info('Assume thetal_nud=theta_nud')
            thlnud = self.variables['theta_nud'].data
        else:
            logger.error('To compute thetal_nud, ta_nud or theta_nud must be known')
            raise ValueError('To compute thetal_nud, ta_nud or theta_nud must be known')

        return thlnud

    def compute_qv_nud(self):

        if 'qt_nud' in self.var_forcing_list:
            logger.info('Assume qv_nud=qt_nud')
            qvnud = self.variables['qt_nud'].data
        elif 'rv_nud' in self.var_forcing_list:
            logger.info('Compute qv_nud from rv_nud')
            rvnud = self.variables['rv_nud'].data
            qvnud =  thermo.rt2qt(rvnud)
        elif 'rt_nud' in self.var_forcing_list:
            logger.info('Compute qv_nud from rv_nud assuming rv_nud=rt_nud')
            rvnud = self.variables['rt_nud'].data
            qvnud =  thermo.rt2qt(rvnud)
        else:
            logger.error('To compute qv_nud, qt_nud, rv_nud or rt_nud must be known')
            raise ValueError('To compute qv_nud, qt_nud, rv_nud or rt_nud must be known')

        return qvnud

    def compute_qt_nud(self):

        if 'qv_nud' in self.var_forcing_list:
            logger.info('Assume qt_nud=qv_nud')
            qtnud = self.variables['qv_nud'].data
        elif 'rv_nud' in self.var_forcing_list:
            logger.info('Compute qt_nud from rt_nud assuming rt_nud=rv_nud')
            rtnud = self.variables['rv_nud'].data
            qtnud =  thermo.rt2qt(rtnud)
        elif 'rt_nud' in self.var_forcing_list:
            logger.info('Compute qt_nud from rt_nud')
            rtnud = self.variables['rt_nud'].data
            qtnud =  thermo.rt2qt(rtnud)
        else:
            logger.error('To compute qt_nud, qv_nud, rv_nud or rt_nud must be known')
            raise ValueError('To compute qt_nud, qv_nud, rv_nud or rt_nud must be known')

        return qtnud

    def compute_rv_nud(self):

        if 'qv_nud' in self.var_forcing_list:
            logger.info('Compute rv_nud from qv_nud')
            qvnud = self.variables['qv_nud'].data
            rvnud =  thermo.qt2rt(qvnud)
        elif 'qt_nud' in self.var_forcing_list:
            logger.info('Compute rv_nud from qv_nud assuming qv_nud=qt_nud')
            qvnud = self.variables['qt_nud'].data
            rvnud =  thermo.qt2rt(qvnud)
        elif 'rt_nud' in self.var_forcing_list:
            logger.info('Assume rv_nud=rt_nud')
            rvnud = self.variables['rt_nud'].data
        else:
            logger.error('To compute rv_nud, qv_nud, qt_nud or rt_nud must be known')
            raise ValueError('To compute rv_nud, qv_nud, qt_nud or rt_nud must be known')

        return rvnud

    def compute_rt_nud(self):

        if 'qv_nud' in self.var_forcing_list:
            logger.info('Compute rt_nud from qt_nud assuming qt_nud=qv_nud')
            qtnud = self.variables['qv_nud'].data
            rtnud =  thermo.qt2rt(qtnud)
        elif 'qt_nud' in self.var_forcing_list:
            logger.info('compute rt_nud from qt_nud')
            qtnud = self.variables['qt_nud'].data
            rtnud =  thermo.qt2rt(qtnud)
        elif 'rv_nud' in self.var_forcing_list:
            logger.info('Assume rt_nud=rv_nud')
            rtnud = self.variables['rv_nud'].data
        else:
            logger.error('To compute rt_nud, qv_nud, qt_nud or rv_nud must be known')
            raise ValueError('To compute rt_nud, qv_nud, qt_nud or rv_nud must be known')

        return rtnud

###################################################################################################
#                  Case interpolation
###################################################################################################

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
            logger.warning('No time interpolation')
            for var in self.var_init_list + self.var_forcing_list:
                VV = self.variables[var]
                dataout[var] = Variable(var, data=VV.data, name=VV.name, units=VV.units,
                        height=VV.height, pressure=VV.pressure,
                        level=VV.level, time=VV.time,
                        plotcoef=VV.plotcoef, plotunits=VV.plotunits)
                if VV.time.id != 't0':
                    dataout[var].time.id = 'time'
        else:
            timeout = Axis('time',time,name='forcing_time',units=self.tunits, calendar='gregorian')
            for var in self.var_init_list + self.var_forcing_list:
                VV = self.variables[var]
                if VV.time.id != 't0':
                    dataout[var] = VV.interpol_time(time=timeout)
                else:
                    dataout[var] = Variable(var, data=VV.data, name=VV.name, units=VV.units,
                        height=VV.height, pressure=VV.pressure,
                        level=VV.level, time=VV.time,
                        plotcoef=VV.plotcoef, plotunits=VV.plotunits)

        if lev is None:
            logger.warning('No vertical interpolation')
            for var in self.var_init_list:
                VV = dataout[var]

                if VV.level is not None:
                    if VV.level.units == 'm':
                        levout = Axis('lev',VV.level.data,name='height',units='m')
                    elif VV.level.units == 'Pa':
                        levout = Axis('lev',VV.level.data,name='air_pressure',units='Pa')
                    else:
                        logger.error('ERROR: Level type undefined for level units: {}'.format(VV.level.units))
                        raise ValueError('ERROR: Level type undefined for level units: {}'.format(VV.level.units))

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
                        logger.error('Level type undefined for level units: {0}'.format(VV.level.units))
                        raise ValueError('Level type undefined for level units: {0}'.format(VV.level.units))

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
                        dataout[var].set_coordinates('t0','zh','lat','lon')
                        dataout[var].height.set_coordinates('t0','zh','lat','lon')

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

            elif levtype == 'pressure':

                levout = Axis('lev', lev, name='air_pressure', units='Pa')
                logger.error('Pressure level type is not coded yet for interpolation')
                raise NotImplementedError('Pressure level type is not coded yet for interpolation')

            else:

                logger.error('levtype unexpected: {0}'.forma(levtype))
                raise ValueError('levtype unexpected: {0}'.forma(levtype))

        for var in self.var_init_list:
            newcase.var_init_list.append(var)
            newcase.variables[var] = dataout[var]

        for var in self.var_forcing_list:
            newcase.var_forcing_list.append(var)
            newcase.variables[var] = dataout[var]

        return newcase

###################################################################################################
#                  Adding missing initial variables
###################################################################################################

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
            logger.error('ta, theta or thetal must be given')
            raise ValueError('ta, theta or thetal must be given')

        VV = self.variables[var]
        levAxis = VV.level

        if VV.height is None and VV.pressure is None:

            logger.error('height and pressure are None for {0}'.format(var))
            raise ValueError('height and pressure are None for {0}'.format(var))

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
            
            logger.warning('Nothing to do. height and pressure already defined. Just pass')

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
                logger.error('Case unexpected: variable {0} have to be defined'.format(var))
                raise ValueError('Case unexpected: variable {0} have to be defined'.format(var))

###################################################################################################
#                  Adding missing forcing variables
###################################################################################################

    def add_missing_forcing_variables(self):

        if not(self.var_forcing_list):
            logger.warning('No forcing variable. Nothing to do')
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
                nt, = VV.time.data.shape
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

        #---- Surface temperature an potential temperature forcing
        att = 'surface_forcing_temp'
        if att in self.attlist and self.attributes[att] in ['thetas','ts']:
            if 'ts_forc' not in self.var_forcing_list and 'thetas_forc' in self.var_forcing_list:
                tmp = thermo.theta2t(p=self.variables['ps_forc'].data, theta=self.variables['thetas_forc'].data)
                self.add_variable('ts_forc', tmp, time=time)
                self.add_init_ts(tmp[0])
            elif 'ts_forc' in self.var_forcing_list and 'thetas_forc' not in self.var_forcing_list:
                tmp = thermo.t2theta(p=self.variables['ps_forc'].data, temp=self.variables['ts_forc'].data)
                self.add_variable('thetas_forc', tmp, time=time)
                self.add_init_thetas(tmp[0])
            self.attributes[att] = 'ts' # Assume this is default, even though ts and thetas are available

        #---- Height/pressure
        if height is None and pressure is None:

            logger.warning('No 3D forcing')
            #logger.error('height and pressure are None. Unexpected in add_missing_forcing_variables')
            #raise ValueError('height and pressure are None. Unexpected in add_missing_forcing_variables')

        elif pressure is None:

            logger.info('Assume pa_forc is constant over time')
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

        elif height is None:

            logger.info('Assume zh_forc is constant over time')
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
            
            logger.warning('Nothing to do. height and pressure already defined. Just pass')

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
            if att in self.attlist:
                if self.attributes[att] > 0:
                    height_loc = np.squeeze(self.variables['zh'].data)
                    pressure_loc = np.squeeze(self.variables['pa'].data)
                    # Add further description of nudging altitude/pressure
                    if 'zh_nudging_{0}'.format(var) not in self.attributes.keys():
                        zlev = thermo.plev2zlev(self.attributes['pa_nudging_{0}'.format(var)],height_loc,pressure_loc)
                        self.set_attribute('zh_nudging_{0}'.format(var),zlev)
                    if 'pa_nudging_{0}'.format(var) not in self.attributes.keys():
                        plev = thermo.zlev2plev(self.attributes['zh_nudging_{0}'.format(var)],height_loc,pressure_loc)
                        self.set_attribute('pa_nudging_{0}'.format(var),plev)
                elif self.attributes[att] == -1:
                    # the nudging profile has already been interpolated
                    pass
                elif self.attributes[att] == 0:
                    # No nudging
                    pass
                else:
                    logger.error('Case unexpected: {0}={1}'.format(att,self.attributes[att]))
                    raise ValueError('Case unexpected: {0}={1}'.format(att,self.attributes[att]))

        #---- Temperature nudging

        flag = False
        zlev = None
        plev = None
        nudging_coefficient = None
        for var in ['ta','theta','thetal']:
            att='nudging_'+var
            if att in self.attlist:
                if self.attributes[att] > 0:
                    flag = True
                    nudging_timescale = self.attributes[att]
                    if self.attributes[att] > 0: # simple nudging profile described in global attributes
                        if 'zh_{0}'.format(att) in self.attributes:
                            zlev = self.attributes['zh_{0}'.format(att)]
                        if 'pa_{0}'.format(att) in self.attributes:
                            plev = self.attributes['pa_{0}'.format(att)]
                elif self.attributes[att] == -1:
                    flag = True
                    nudging_coefficient = self.variables['nudging_coefficient_'+var]
                elif self.attributes[att] == 0:
                    # No nudging
                    pass
                else:
                    logger.error('Case unexpected: {0}={1}'.format(att,self.attributes[att]))
                    raise ValueError('Case unexpected: {0}={1}'.format(att,self.attributes[att]))

        if flag:
            if zlev is not None or plev is not None: # simple nudging profile described in global attributes
                # update nudging height/pressure for all temperature variables
                height_loc = np.squeeze(self.variables['zh'].data)
                pressure_loc = np.squeeze(self.variables['pa'].data)
                if zlev is None:
                    zlev = int(thermo.plev2zlev(plev,height_loc,pressure_loc))
                if plev is None:
                    plev = int(thermo.zlev2plev(zlev,height_loc,pressure_loc))

                for var in ['ta','theta','thetal']:
                    self.set_attribute('nudging_{0}'.format(var),nudging_timescale)
                    self.set_attribute('zh_nudging_{0}'.format(var),zlev)
                    self.set_attribute('pa_nudging_{0}'.format(var),plev)
            else:
                for var in ['ta','theta','thetal']:
                    self.add_variable('nudging_coefficient_{0}'.format(var), nudging_coefficient.data,
                                      lev=nudging_coefficient.level, time=nudging_coefficient.time,
                                      height=nudging_coefficient.height, pressure=nudging_coefficient.pressure)
                    self.set_attribute('nudging_{0}'.format(var),-1)

            # temperature nudging is active. All temperature variables are added, if needed.
            if 'ta_nud' not in self.var_forcing_list:
                tnud = self.compute_ta_nud()
                self.add_variable('ta_nud', tnud, 
                    lev=level, time=time,
                    height=height, pressure=pressure)
            if 'theta_nud' not in self.var_forcing_list:
                thnud = self.compute_theta_nud()
                self.add_variable('theta_nud', thnud, 
                    lev=level, time=time,
                    height=height, pressure=pressure)
            if 'thetal_nud' not in self.var_forcing_list:
                thlnud = self.compute_thetal_nud()
                self.add_variable('thetal_nud', thlnud, 
                    lev=level, time=time,
                    height=height, pressure=pressure)

        #---- Humidity nudging

        flag = False
        zlev = None
        plev = None
        nudging_coefficient = None
        for var in ['qv','qt','rv','rt']:
            att='nudging_'+var
            if att in self.attlist:
                if self.attributes[att] > 0:
                    flag = True
                    nudging_timescale = self.attributes[att]
                    if self.attributes[att] > 0: # simple nudging profile described in global attributes
                        if 'zh_{0}'.format(att) in self.attributes:
                            zlev = self.attributes['zh_{0}'.format(att)]
                        if 'pa_{0}'.format(att) in self.attributes:
                            plev = self.attributes['pa_{0}'.format(att)]
                elif self.attributes[att] == -1:
                    flag = True
                    nudging_coefficient = self.variables['nudging_coefficient_'+var]
                elif self.attributes[att] == 0:
                    # No nudging
                    pass
                else:
                    logger.error('Case unexpected: {0}={1}'.format(att,self.attributes[att]))
                    raise ValueError('Case unexpected: {0}={1}'.format(att,self.attributes[att]))

        if flag:
            if zlev is not None or plev is not None: # simple nudging profile described in global attributes
                # update nudging height/pressure for all humidity variables
                height_loc = np.squeeze(self.variables['zh'].data)
                pressure_loc = np.squeeze(self.variables['pa'].data)
                if zlev is None:
                    zlev = int(thermo.plev2zlev(plev,height_loc,pressure_loc))
                if plev is None:
                    plev = int(thermo.zlev2plev(zlev,height_loc,pressure_loc))

                for var in ['qv','qt','rv','rt']:
                    self.set_attribute('nudging_{0}'.format(var),nudging_timescale)
                    self.set_attribute('zh_nudging_{0}'.format(var),zlev)
                    self.set_attribute('pa_nudging_{0}'.format(var),plev)
            else:
                for var in ['qv','qt','rv','rt']:
                    self.add_variable('nudging_coefficient_{0}'.format(var), nudging_coefficient.data,
                                      lev=nudging_coefficient.level, time=nudging_coefficient.time,
                                      height=nudging_coefficient.height, pressure=nudging_coefficient.pressure)
                    self.set_attribute('nudging_{0}'.format(var),-1)

            # humidity nudging is active. All temperature variables are added, if needed.
            if 'qv_nud' not in self.var_forcing_list:
                qvnud = self.compute_qv_nud()
                self.add_variable('qv_nud', qvnud, 
                    lev=level, time=time,
                    height=height, pressure=pressure)
            if 'qt_nud' not in self.var_forcing_list:
                qtnud = self.compute_qt_nud()
                self.add_variable('qt_nud', qtnud, 
                    lev=level, time=time,
                    height=height, pressure=pressure)
            if 'rv_nud' not in self.var_forcing_list:
                rvnud = self.compute_rv_nud()
                self.add_variable('rv_nud', rvnud, 
                    lev=level, time=time,
                    height=height, pressure=pressure)
            if 'rt_nud' not in self.var_forcing_list:
                rtnud = self.compute_rt_nud()
                self.add_variable('rt_nud', rtnud, 
                    lev=level, time=time,
                    height=height, pressure=pressure)

###################################################################################################
#                  Case conversion to SCM-enabled format
###################################################################################################

    def convert2SCM(self, time=None, lev=None, levtype=None,
            usetemp=True, usetheta=True, usethetal=True):

        # Interpolation
        logger.info('#'*40)
        logger.info('#### Interpolate available variables')

        caseSCM = self.interpolate(time=time, lev=lev, levtype=levtype,
                usetemp=usetemp, usetheta=usetheta, usethetal=usethetal)

        # Add missing variables for initial state
        logger.info('#'*40)
        logger.info('#### Add missing initial variables')

        caseSCM.add_missing_init_variables()

        # Add missing forcing variables
        logger.info('#'*40)
        logger.info('#### Add missing forcing variables')

        caseSCM.add_missing_forcing_variables()

        # Final
        logger.info('#'*40)

        return caseSCM

###################################################################################################
#                  Vertically extended variables
###################################################################################################

    def extend_variable(self, varid, data=None, height=None, pressure=None, time=None, tunits=None):
        """Vertically extend the variable varid using data

        Arguments:
        varid -- string for the variable id
        data -- data to be used to extend vertically the variable varid
        height -- altitude above the reference geoide (m)
        pressure -- pressure (Pa)
        time -- time stamps
        tunits -- units to interpret time data

        Either height or pressure must be given

        if data is None, the variable is extended constantly from its highest or lowest value.
        if data is a float or an integer, the variable is extended by on level based on either (data,height) 
            or (data,pressure). height or pressure must be a float or an integer.
        if data is a 1D array, it is supposed to be a vertical profile. Either height or pressure must be given 
            and must have the same shape.
        if data is a 2D array, first dimension is supposed to be time and second to be either height or pressure. 
            height or pressure must be given and could be 1D or 2D (shape to be consistent with data)
            time and tunits must be given. time must have a shape consistent with the data first dimention. 
        """

        if height is None and pressure is None:
            logger.error('Either height or pressure must be given')
            raise ValueError('Either height or pressure must be given')

        if data is None:
            _height = height
            _pressure = pressure
            _time = None
            _tunits = None
            hmin = np.min(self.variables[varid].height.data)
            hmax = np.max(self.variables[varid].height.data)
            if np.min(height) < hmax and np.max(height) > hmin:
                logger.error("Case unexpected: cannot add levels up and down at the same time")
                raise NotImplementedError

            try:
                if np.max(height) <= hmin:
                    _data = float(self.variables[varid].data[:,0])
                elif np.min(height) >= hmax:
                    _data = float(self.variables[varid].data[:,-1])
                else:
                    logger.error("You should not be in this case")
                    raise ValueError
            except TypeError:
                _data = self.variables[varid].data[:,-1]
                _data = _data[:,np.newaxis]
                _time = self.variables[varid].time.data
                _tunits = self.variables[varid].time.units
                #print(varid, _data.shape)
                if height is not None and (isinstance(height,float) or isinstance(height,int)):
                    _height = np.zeros(_data.shape, dtype=np.float32) + height
                elif pressure is not None and (isinstance(pressure,float) or isinstance(pressure,int)):
                    _pressure = np.zeros(_data.shape, dtype=np.float32) + pressure
                else:
                    logger.error("Case unexpected in extended variable {0}".format(varid))
                    raise NotImplementedError
            except:
                raise

            self.extend_variable(varid, data=_data, height=_height, pressure=_pressure,
                    time=_time, tunits=_tunits)

        if isinstance(data,float) or isinstance(data,int):
            if height is not None:
                if isinstance(height,float) or isinstance(height,int):
                    _data = np.array([0,data])
                    _height = np.array([0,height])
                    self.extend_variable(varid, data=_data, height=_height)
                else:
                    logger.error('height and data are incompatible (data is float)')
                    raise ValueError('height and data are incompatible (data is float)')
     
            if pressure is not None:
                if isinstance(pressure,float) or isinstance(pressure,int):
                    _data = np.array([0,data])
                    _pressure = np.array([self.variables['ps'].data[0],pressure])
                    self.extend_variable(varid, data=_data, pressure=_pressure)
                else:
                    logger.error('pressure and data are incompatible (data is float)')
                    raise ValueError('pressure and data are incompatible (data is float)')

        if isinstance(data,list): data = np.array(data)
        if isinstance(height,list): height = np.array(height)
        if isinstance(pressure,list): pressure = np.array(pressure)
        if isinstance(time,list): time = np.array(time)
        
        if isinstance(data,np.ndarray):
            if len(data.shape) == 1:
                # Checking consistency
                if height is not None and height.shape != data.shape:
                    logger.error('Incompatibility between data and height shapes')
                    raise ValueError('Incompatibility between data and height shapes')
                if pressure is not None and pressure.shape != data.shape:
                    logger.error('Incompatibility between data and pressure shapes')
                    raise ValueError('Incompatibility between data and pressure shapes')

                # Extending data
                nt, _ = self.variables[varid].data.shape
                _time = self.variables[varid].time.data
                _tunits = self.variables[varid].time.units
                _data = np.tile(data,(nt,1))
                _height = height
                _pressure = pressure
                if height is not None:
                    _height = np.tile(height,(nt,1))
                if pressure is not None:
                    _pressure = np.tile(pressure,(nt,1))
                #print('CASE',_data.shape)
                _varnew = self.variables[varid].extend_vert(data=_data, height=_height, pressure=_pressure,
                                                            time=_time, tunits=_tunits)

                self.variables[varid] = _varnew

            elif len(data.shape) == 2:
                # Checking consistency
                if time is None and tunits is None:
                    logger.error('data is 2D. time and tunits must be given')
                    raise ValueError('data is 2D. time and tunits must be given')
                if data.shape[0] != time.shape[0]:
                    logger.error('Incompatibility between data and time shapes')
                    raise ValueError('Incompatibility between data and time shapes')
                if (height is not None and len(height.shape) == 2 and height.shape != data.shape)\
                        or (height is not None and len(height.shape) == 1 and height.shape[0] != data.shape[1]):
                    logger.error('Incompatibility between data and height shapes')
                    raise ValueError('Incompatibility between data and height shapes')
                if (pressure is not None and len(pressure.shape) == 2 and pressure.shape != data.shape)\
                        or (pressure is not None and len(pressure.shape) == 1 and pressure.shape[0] != data.shape[1]):
                    logger.error('Incompatibility between data and pressure shapes')
                    raise ValueError('Incompatibility between data and pressure shapes')

                # Extending data
                _varnew = self.variables[varid].extend_vert(data=data, height=height, pressure=pressure,
                                                            time=time, tunits=tunits)

                self.variables[varid] = _varnew

            else:
                logger.error('shape unexpected for data : {0}'.format(data.shape))
                raise ValueError


    def extend_init_wind(self, u=None, v=None, **kwargs):
        """Vertically extend the two wind initial components"""

        self.extend_variable('ua', data=u, **kwargs)
        self.extend_variable('va', data=v, **kwargs)

    def extend_init_temp(self, temp=None, **kwargs):
        """Vertically extend the temperarture"""

        self.extend_variable('ta', data=temp, **kwargs)

    def extend_init_theta(self, theta=None, **kwargs):
        """Vertically extend the potential temperature"""

        self.extend_variable('theta', data=theta, **kwargs)

    def extend_init_thetal(self, thetal=None, **kwargs):
        """Vertically extend the liquid-water potential temperature"""

        self.extend_variable('thetal', data=thetal, **kwargs)

    def extend_init_qv(self, qv=None, **kwargs):
        """Vertically extend the specific humidity"""

        self.extend_variable('qv', data=qv, **kwargs)

    def extend_init_qt(self, qt=None, **kwargs):
        """Vertically extend the total water"""

        self.extend_variable('qt', data=qt, **kwargs)

    def extend_init_rv(self, rv=None, **kwargs):
        """Vertically extend the water vapor mixing ratio"""

        self.extend_variable('rv', data=rv, **kwargs)

    def extend_init_rt(self, rt=None, **kwargs):
        """Vertically extend the total water mixing ratio"""

        self.extend_variable('rt', data=rt, **kwargs)

    def extend_geostrophic_wind(self, ug=None, vg=None, **kwargs):
        """Vertically extend the geostrophic wind components"""

        self.extend_variable('ug', data=ug, **kwargs)
        self.extend_variable('vg', data=vg, **kwargs)

    def extend_vertical_velocity(self, w=None, omega=None, **kwargs):
        """Vertically extend the vertical velocity, either w or omega"""

        if w is not None:
            self.extend_variable('wa', data=w, **kwargs)
        elif omega is not None:
            self.extend_variable('wap', data=omega, **kwargs)
        else:
            logger.error('At least w or omega must be given')
            raise ValueError

    def extend_temperature_advection(self, temp_adv=None, **kwargs):
        """Vertically extend the temperature large-scale advection"""

        self.extend_variable('tnt_adv', data=temp_adv, **kwargs)

    def extend_theta_advection(self, theta_adv=None, **kwargs):
        """Vertically extend the potential temperature large-scale advection"""

        self.extend_variable('tntheta_adv', data=theta_adv, **kwargs)

    def extend_thetal_advection(self, thetal_adv=None, **kwargs):
        """Vertically extend the liquid-water potential temperature large-scale advection"""

        self.extend_variable('tnthetal_adv', data=thetal_adv, **kwargs)

    def extend_qv_advection(self, qv_adv=None, **kwargs):
        """Vertically extend the specific humidity large-scale advection"""

        self.extend_variable('tnqv_adv', data=qv_adv, **kwargs)

    def extend_qt_advection(self, qt_adv=None, **kwargs):
        """Vertically extend the total water content large-scale advection"""

        self.extend_variable('tnqt_adv', data=qt_adv, **kwargs)

    def extend_rv_advection(self, rv_adv=None, **kwargs):
        """Vertically extend the water vapor mixing ratio large-scale advection"""

        self.extend_variable('tnrv_adv', data=rv_adv, **kwargs)

    def extend_rt_advection(self, rt_adv=None, **kwargs):
        """Vertically extend the total water mixing ratio large-scale advection"""

        self.extend_variable('tnrt_adv', data=rt_adv, **kwargs)
