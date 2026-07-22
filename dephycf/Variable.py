#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Variable module.

This module defines the :class:`Variable` class, a container that
associates a data array with its axes (time, level, latitude,
longitude) and, optionally, with a companion height or pressure
vertical coordinate (itself stored as a :class:`Variable`).

It also provides module-level helper functions to:

- read a :class:`Variable` from an open netCDF file (:func:`read`);
- interpolate a :class:`Variable` in time and/or in the vertical
  (:func:`interpol`).
"""

import os
import logging
import time

import netCDF4 as nc

import numpy as np
from scipy import interpolate

from .Axis import Axis
from .variables_attributes import attributes as var_attributes
from . import plotbasics

logger = logging.getLogger(__name__)


class Variable:
    """Container for a physical variable and its associated axes.

    A :class:`Variable` bundles a data array together with the axes
    it is defined on (typically ``time`` and/or ``level``), plus
    optional vertical-coordinate companions (``height`` and/or
    ``pressure``, themselves stored as :class:`Variable` instances).
    It also stores plotting-related metadata (`plotcoef`,
    `plotunits`) and provides methods to read/write netCDF data, plot
    the variable, and interpolate it in time and in the vertical.

    Attributes
    ----------
    id : :class:`str`
        Variable identifier (netCDF variable name).
    name : :class:`str` or `None`
        Standard/long name of the variable.
    units : :class:`str` or `None`
        Physical units of the variable.
    data : :class:`numpy.ndarray` or `None`
        Data array associated with the variable.
    sh : :class:`int`
        Number of dimensions of `data` (only set when `data` is not
        None).
    axes : :class:`list`
        List of :class:`Axis` (or :class:`Variable`) objects the
        variable is defined on.
    axlist : :class:`str`
        List of axis identifiers, in the same order as `axes`.
    time : :class:`Axis` or `None`
        Time axis, if any.
    level : :class:`Axis` or `None`
        Level (vertical) axis, if any.
    coord : :class:`str`
        Space-separated coordinates string (netCDF ``coordinates``
        attribute).
    plotcoef : :class:`float`
        Multiplicative coefficient applied to the data before
        plotting.
    plotunits : :class:`str` or `None`
        Units used when plotting (may differ from `units`).
    height : :class:`Variable` or `None`
        Height companion variable, used as an alternative vertical
        coordinate.
    pressure : :class:`Variable` or `None`
        Pressure companion variable, used as an alternative vertical
        coordinate.
    """

    def __init__(self, varid, data=None, name=None, units=None,
                 height=None, height_id=None, height_units=None,
                 pressure=None, pressure_id=None, pressure_units=None,
                 level=None, time=None, lat=None, lon=None,
                 axlist=None, axes=None,
                 plotcoef=1., plotunits=None):
        """Build a :class:`Variable` instance.

        Parameters
        ----------
        varid : :class:`str`
            Variable identifier (netCDF variable name).
        data : array_like, optional
            Data values. Converted to a ``numpy.float64`` array if
            given.
        name : :class:`str`, optional
            Standard/long name of the variable.
        units : :class:`str`, optional
            Physical units of the variable.
        height : :class:`Axis`, :class:`Variable` or array_like, optional
            Height vertical coordinate (or companion variable).
            Privileged over `pressure` when both are given.
        height_id : :class:`str`, optional
            Identifier to use for the height companion variable, when
            `height` is given as raw data rather than as an
            :class:`Axis`/:class:`Variable`.
        height_units : :class:`str`, optional
            Units of `height`, when given as raw data.
        pressure : :class:`Axis`, Variable or array_like, optional
            Pressure vertical coordinate (or companion variable).
        pressure_id : :class:`str`, optional
            Identifier to use for the pressure companion variable,
            when `pressure` is given as raw data.
        pressure_units : :class:`str`, optional
            Units of `pressure`, when given as raw data.
        level : :class:`Axis`, optional
            Level (vertical) axis. Ignored if `axes` is given.
        time : :class:`Axis`, optional
            Time axis. Ignored if `axes` is given.
        lat : :class:`float` or array-like, optional
            Latitude information 
        lon : :class:`float` or array-like,optional
            Longitude information 
        axlist : :class:`list`, optional
            Explicit list of axis identifiers. Used together with
            `axes` when provided.
        axes : :class:`list`, optional
            Explicit list of :class:`Axis` objects. When given,
            `time` and `level` axes are inferred from it instead of
            from the `time`/`level` parameters.
        plotcoef : :class:`float`, optional
            Default multiplicative coefficient applied before
            plotting. Overridden by a value found in
            `var_attributes`. Defaults to ``1.``.
        plotunits : :class:`str`, optional
            Default plotting units. Overridden by a value found in
            `var_attributes` when not given.
        """
        self.id = varid
        self.units = units
        self.name = name

        self.data = None
        if data is not None:
            self.data = np.array(data, dtype=np.float64)
            self.sh = len(self.data.shape)

        if axes is None:
            self.axes = []
            self.axlist = []

            self.time = time
            self.level = level

            for ax in ['time', 'level']:
                if not (self.__dict__[ax] is None):
                    self.axes.append(self.__dict__[ax])
                    self.axlist.append(self.__dict__[ax].id)
        else:
            self.axes = axes
            self.axlist = axlist

            self.time = None
            self.level = None

            for ax in axes:
                if ax.id == 't0' or ax.id[0:4] == 'time':
                    self.time = ax
                elif ax.id[0:3] == 'lev' or ax.id == 'nlev':
                    self.level = ax
                else:
                    logger.error('Axis unexpected: {0}'.format(ax.id))
                    if logger.level == logging.DEBUG:
                        ax.info()
                    raise ValueError('Axis unexpected: {0}'.format(ax.id))

        self.lat = lat
        self.lon = lon
        self.coord = " ".join(self.axlist + ['lat', 'lon'])

        try:
            self.plotcoef = var_attributes[self.id]['plotcoef']
        except KeyError:
            self.plotcoef = plotcoef
        except:
            raise

        if plotunits is None:
            try:
                self.plotunits = var_attributes[self.id]['plotunits']
            except KeyError:
                self.plotunits = self.units
            except:
                raise
        else:
            self.plotunits = plotunits

        self.height = None
        self.pressure = None

        if height is not None:  # height is privileged over pressure
            if isinstance(height, Axis) or isinstance(height, Variable):
                self.height = Variable(
                    height.id, data=height.data,
                    units=height.units, name=height.name, 
                    level=self.level, time=self.time)
                self.coord = " ".join(
                    [self.time.id, height.id, 'lat', 'lon'])
            else:
                if height_id is None:
                    height_id = 'zh_{0}'.format(self.id)
                height_name = 'height_for_{0}'.format(self.id)
                if height_units is None:
                    height_units = 'm'
                self.height = Variable(
                    height_id, data=height, 
                    units=height_units, name=height_name,
                    level=self.level, time=self.time)
                self.coord = " ".join(
                    [self.time.id, height_id, 'lat', 'lon'])

            self.height.set_coordinates(self.coord)

        elif pressure is not None:
            if isinstance(pressure, Axis) or isinstance(pressure, Variable):
                self.pressure = Variable(
                    pressure.id, data=pressure.data,
                    units=pressure.units, name=pressure.name,
                    level=self.level, time=self.time)
                self.coord = " ".join(
                    [self.time.id, pressure.id, 'lat', 'lon'])
            else:
                if pressure_id is None:
                    pressure_id = 'pa_{0}'.format(self.id)
                pressure_name = 'air_pressure_for_{0}'.format(self.id)
                if pressure_units is None:
                    pressure_units = 'Pa'
                self.pressure = Variable(
                    pressure_id, data=pressure,
                    units=pressure_units, name=pressure_name,
                    level=self.level, time=self.time)
                self.coord = " ".join(
                    [self.time.id, pressure_id, 'lat', 'lon'])

            self.pressure.set_coordinates(self.coord)

    def info(self):
        """Print a short summary of the variable to stdout.

        Displays the identifier, name, units, axes, coordinates
        string, and basic statistics (mean, min, max) of the data.
        """
        print('-' * 5, 'Variable:', self.id)
        print('-' * 10, 'Name:', self.name)
        print('-' * 10, 'Units:', self.units)
        print('-' * 10, 'Axes:', self.axlist)
        print('-' * 10, 'Coordinates:', self.coord)
        print('-' * 10, 'mean: {0:13f}; min: {1:13f}; max: {2:13f}'.format(
            np.average(self.data), np.amin(self.data), np.amax(self.data)))

    def set_coordinates(self, *coord):
        """Set the netCDF ``coordinates`` attribute string.

        Parameters
        ----------
        *coord : str
            Coordinate identifiers to join with a single space to
            build `self.coord`.
        """
        self.coord = " ".join(coord)

    def set_level(self, lev=None):
        """Replace the level axis of the variable.

        Parameters
        ----------
        lev : Axis, optional
            New level axis. If ``None``, a warning is logged and
            nothing is done.
        """
        if lev is None:
            logger.warning('level is None. Nothing to do')
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
        """Write the variable (and, optionally, its axes) to netCDF.

        Parameters
        ----------
        filein : netCDF4.Dataset
            Open (writable) netCDF file/dataset to write into.
        write_time_axes : bool, optional
            Whether to write the time axis/axes. Defaults to
            ``True``.
        write_level_axes : bool, optional
            Whether to write the level axis/axes, when the variable
            has a vertical dimension. Defaults to ``True``.
        write_data : bool, optional
            Whether to write the variable data itself. Defaults to
            ``True``.
        write_vertical : bool, optional
            Whether to also write the height/pressure companion
            variables, when the variable has a vertical dimension.
            Defaults to ``True``.
        """
        lvert = False
        for ax in self.axes:
            if ax.id[0:3] == 'lev':
                lvert = True

        if write_time_axes:
            for ax in self.axes:
                if ax.id[0:4] == 'time' or ax.id == 't0':
                    ax.write(filein)

        if write_level_axes and lvert:
            for ax in self.axes:
                if ax.id[0:3] == 'lev':
                    ax.write(filein)

        if write_data:
            if self.data is not None:
                if self.id in filein.variables:
                    logger.debug(
                        '{0} already if netCDF file. Not overwritten'.format(
                            self.id))
                else:
                    tmp = filein.createVariable(
                        self.id, np.float32, self.axlist)
                    tmp[:] = self.data
                    tmp.standard_name = self.name
                    tmp.units = self.units
                    tmp.coordinates = self.coord

        if write_vertical and lvert:
            if self.height is not None:
                if self.height.id in filein.variables:
                    logger.debug(
                        '{0} already if netCDF file. Not overwritten'.format(
                            self.height.id))
                else:
                    self.height.write(filein)

            if self.pressure is not None:
                if self.pressure.id in filein.variables:
                    logger.debug(
                        '{0} already if netCDF file. Not overwritten'.format(
                            self.pressure.id))
                else:
                    self.pressure.write(filein)

    def plot(self, rep_images=None, var2=None, label="", label2="",
             timeunits=None, levunits=None):
        """Plot the variable, optionally overlaid with a second one.

        Depending on which axes (`time`, `level`) the variable has,
        this produces a vertical profile, a time series, or a
        Hovmoller-like time/height plot, delegating the actual
        drawing to :mod:`plotbasics`.

        Parameters
        ----------
        rep_images : str, optional
            Output directory for the generated image(s).
        var2 : Variable, optional
            A second variable to overlay on the same plot for
            comparison.
        label : str, optional
            Legend label for `self` when `var2` is given.
        label2 : str, optional
            Legend label for `var2` when `var2` is given.
        timeunits : {'hours', 'days'}, optional
            Units to convert the time axis to before plotting. If
            ``None``, the native time units are used.
        levunits : {'hPa', 'Pa', 'km', 'm'}, optional
            Units to convert the level axis to before plotting. If
            ``None``, the native level units are used.
        """
        coef = self.plotcoef

        if not (self.level is None):
            if levunits is None:
                levs = self.level.data
                levunits = self.level.units
                zlabel = 'Altitude above the surface [{0}]'.format(levunits)
            elif levunits == 'hPa' and self.level.units == 'Pa':
                levs = self.level.data / 100.
                zlabel = 'Pressure [{0}]'.format(levunits)
            elif (levunits == 'hPa' and self.level.units == 'hPa') or \
                    (levunits == 'Pa' and self.level.units == 'Pa'):
                levs = self.level.data
                zlabel = 'Pressure [{0}]'.format(levunits)
            elif levunits == 'km' and self.level.units == 'm':
                levs = self.level.data / 1000.
                zlabel = 'Altitude above the surface [{0}]'.format(levunits)
            elif (levunits == 'km' and self.level.units == 'km') or \
                    (levunits == 'm' and self.level.units == 'm'):
                levs = self.level.data
                zlabel = 'Altitude above the surface [{0}]'.format(levunits)
            else:
                logger.error(
                    "Unexpected case for levunits: {0} {1}".format(
                        levunits, self.level.units))
                raise ValueError(
                    "Unexpected case for levunits: {0} {1}".format(
                        levunits, self.level.units))
            if not (var2 is None):
                if levunits is None:
                    levs2 = var2.level.data
                elif levunits == 'hPa' and var2.level.units == 'Pa':
                    levs2 = var2.level.data / 100.
                elif (levunits == 'hPa' and var2.level.units == 'hPa') or \
                        (levunits == 'Pa' and var2.level.units == 'Pa'):
                    levs2 = var2.level.data
                elif levunits == 'km' and var2.level.units == 'm':
                    levs2 = var2.level.data / 1000.
                elif (levunits == 'km' and var2.level.units == 'km') or \
                        (levunits == 'm' and var2.level.units == 'm'):
                    levs2 = var2.level.data
                else:
                    logger.error(
                        "Unexpected case for levunits (var2): "
                        "{0} {1}".format(levunits, var2.level.units))
                    logger.error(
                        "it is possible that the vertical units for DEF "
                        "and SCM files differ and the comparison cannot "
                        "be plotted ")
                    raise ValueError(
                        "Unexpected case for levunits (var2): "
                        "{0} {1}".format(levunits, var2.level.units))

        if not (self.time is None) and not (self.level is None):

            if self.time.length == 1:
                if var2 is None:
                    plotbasics.plot(
                        self.data[0, :] * coef, levs,
                        xlabel='{0} [{1}]'.format(self.id, self.plotunits),
                        ylabel=zlabel,
                        title='{0} ({1})'.format(self.name, self.plotunits),
                        rep_images=rep_images, name='{0}.png'.format(self.id),
                        yunits=levunits)
                else:
                    plotbasics.plot(
                        self.data[0, :] * coef, levs,
                        x2=var2.data[0, :] * coef,
                        y2=levs2, xlabel='{0} [{1}]'.format(
                            self.id, self.plotunits),
                        ylabel=zlabel,
                        title='{0} ({1})'.format(self.name, self.plotunits),
                        rep_images=rep_images, name='{0}.png'.format(self.id),
                        label=label, label2=label2,
                        yunits=levunits)
            else:
                if timeunits is None:
                    time = self.time.data
                    tunits = self.time.units
                else:
                    if timeunits == 'hours':
                        time = self.time.data / 3600.
                        tunits = self.time.units.replace("seconds", "hours")
                    elif timeunits == 'days':
                        time = self.time.data / 86400.
                        tunits = self.time.units.replace("seconds", "days")
                    else:
                        logger.error(
                            "timeunits unexpected for plotting: {0}".format(
                                timeunits))
                        raise NotImplementedError(
                            "timeunits unexpected for plotting: {0}".format(
                                timeunits))

                plotbasics.plot2D(
                    time, levs, self.data[:, :] * coef,
                    xlabel=tunits,
                    ylabel=zlabel,
                    title='{0} [{1}]'.format(self.id, self.plotunits),
                    rep_images=rep_images, name='{0}.png'.format(self.id),
                    yunits=levunits)

        elif not (self.time is None):

            if self.time.length > 1:

                if timeunits is None:
                    time = self.time.data
                    if not (var2 is None):
                        time2 = var2.time.data
                    tunits = self.time.units
                else:
                    if timeunits == 'hours':
                        time = self.time.data / 3600.
                        if not (var2 is None):
                            time2 = var2.time.data / 3600.
                        tunits = self.time.units.replace("seconds", "hours")
                    elif timeunits == 'days':
                        time = self.time.data / 86400.
                        if not (var2 is None):
                            time2 = var2.time.data / 86400.
                        tunits = self.time.units.replace("seconds", "days")
                    else:
                        logger.error(
                            "timeunits unexpected for plotting: {0}".format(
                                timeunits))
                        raise NotImplementedError(
                            "timeunits unexpected for plotting: {0}".format(
                                timeunits))

                if var2 is None:
                    plotbasics.plot(
                        time, self.data[:] * coef,
                        ylabel='{0} [{1}]'.format(self.id, self.plotunits),
                        xlabel=tunits,
                        title=self.name,
                        rep_images=rep_images, name='{0}.png'.format(self.id))
                else:
                    plotbasics.plot(
                        time, self.data[:] * coef,
                        x2=time2,
                        y2=var2.data[:] * coef,
                        ylabel='{0} [{1}]'.format(self.id, self.plotunits),
                        xlabel=tunits,
                        title=self.name,
                        rep_images=rep_images, name='{0}.png'.format(self.id),
                        label=label, label2=label2)
            else:
                logger.warning('no plot for variable {0}'.format(self.id))

        elif not (self.level is None):

            if var2 is None:
                plotbasics.plot(
                    self.data[:] * coef, levs,
                    xlabel='{0} [{1}]'.format(self.id, self.plotunits),
                    ylabel=zlabel,
                    title=self.name,
                    rep_images=rep_images, name='{0}.png'.format(self.id),
                    yunits=levunits)
            else:
                plotbasics.plot(
                    self.data[:] * coef, levs,
                    x2=var2.data[:] * coef,
                    y2=levs2,
                    xlabel='{0} [{1}]'.format(self.id, self.plotunits),
                    ylabel=zlabel,
                    title='{0} ({1})'.format(self.name, self.time.name),
                    rep_images=rep_images, name='{0}.png'.format(self.id),
                    label=label, label2=label2,
                    yunits=levunits)

        else:
            logger.warning('no plot for variable {0}'.format(self.id))

    def interpol_time(self, time=None):
        """Interpolate the variable onto a new time axis.

        Parameters
        ----------
        time : Axis, optional
            Target time axis. If ``None``, a warning is logged and
            `self` is returned unchanged.

        Returns
        -------
        Variable
            A new :class:`Variable` interpolated onto `time` (or
            `self`, unchanged, if no interpolation is needed/possible).

        Raises
        ------
        ValueError
            If the variable does not have a time axis.
        """
        if self.time is None:
            logger.error(
                "Time interpolation requested for variable {0} which "
                "does not have a time axis".format(self.id))
            raise ValueError(
                "Time interpolation requested for variable {0} which "
                "does not have a time axis".format(self.id))

        if time is None:
            logger.warning('time is None. Thus no time interpolation')
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

            logger.debug(
                'Variable "{0}" is an initial state variable. No need '
                'for time interpolation.'.format(self.id))
            return self

        elif l2D:  # time,level variable

            ntin, nlevin = self.data.shape
            data = np.zeros((ntout, nlevin), dtype=np.float64)
            for ilev in range(0, nlevin):
                # Pad after end date with the last value, if necessary
                ff = interpolate.interp1d(
                    self.time.data, self.data[:, ilev],
                    bounds_error=False,
                    fill_value=self.data[-1, ilev])
                data[:, ilev] = ff(time.data)

            if self.height is not None:
                height_id = self.height.id
                height_units = self.height.units
                height = np.zeros((ntout, nlevin), dtype=np.float64)
                for ilev in range(0, nlevin):
                    # Pad after end date with the last value, if necessary
                    ff = interpolate.interp1d(
                        self.height.time.data, self.height.data[:, ilev],
                        bounds_error=False,
                        fill_value=self.height.data[-1, ilev])
                    height[:, ilev] = ff(time.data)

            if self.pressure is not None:
                pressure_id = self.pressure.id
                pressure_units = self.pressure.units
                pressure = np.zeros((ntout, nlevin), dtype=np.float64)
                for ilev in range(0, nlevin):
                    # Pad after end date with the last value, if necessary
                    ff = interpolate.interp1d(
                        self.pressure.time.data,
                        self.pressure.data[:, ilev],
                        bounds_error=False,
                        fill_value=self.pressure.data[-1, ilev])
                    pressure[:, ilev] = ff(time.data)

        else:  # time only variable

            ff = interpolate.interp1d(
                self.time.data, self.data[:],
                bounds_error=False,
                fill_value=self.data[-1])  # Pad after end date with the last value, if necessary
            data = ff(time.data)

        return Variable(
            self.id, data=data, name=self.name, units=self.units,
            level=self.level, time=time,
            height=height, height_id=height_id, height_units=height_units,
            pressure=pressure, pressure_id=pressure_id,
            pressure_units=pressure_units)

    def interpol_vert(self, height=None, pressure=None, log=False):
        """Interpolate the variable onto a new vertical coordinate.

        Parameters
        ----------
        height : array_like, optional
            Target height values (1D, shared for all times, or 2D
            ``(time, level)``). Used only if the variable already has
            a height companion variable.
        pressure : array_like, optional
            Target pressure values (1D or 2D ``(time, level)``). Used
            only if the variable already has a pressure companion
            variable, and only when `height` is not usable.
        log : bool, optional
            Unused placeholder flag kept for API compatibility.
            Defaults to ``False``.

        Returns
        -------
        Variable
            `self` unchanged if there is no level axis or if neither
            `height` nor `pressure` is given; otherwise a new
            :class:`Variable` interpolated onto the target vertical
            coordinate.

        Raises
        ------
        ValueError
            If neither a usable height nor a usable pressure
            companion variable/target is available.
        """
        if self.level is None:
            logger.debug(
                "Vertical interpolation requested for variable {0}, "
                "which does not have a level axis".format(self.id))
            logger.debug("Simply return original variable")
            return self

        if height is None and pressure is None:
            logger.warning(
                "height and pressure are None. Thus no vertical "
                "interpolation")
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
                _height = np.tile(height, (ntin, 1))
            else:
                _height = height

            _, nlevout = _height.shape
            _height_id = self.height.id
            _height_units = self.height.units
            _level = Axis(
                'lev_{0}'.format(self.id), _height[0, :],
                name='height_for_{0}'.format(self.id), units='m')

            data = np.zeros((ntin, nlevout), dtype=np.float64)
            for it in range(0, ntin):
                ff = interpolate.interp1d(
                    self.height.data[it, :], self.data[it, :],
                    bounds_error=False,
                    fill_value=(self.data[it, 0], self.data[it, -1])
                )  # Pad below and above, if necessary
                data[it, :] = ff(_height[it, :])

        elif pressure is not None and self.pressure is not None:

            if len(pressure.shape) == 1:
                _pressure = np.tile(pressure, (ntin, 1))
            else:
                _pressure = pressure

            _, nlevout = _pressure.shape
            _pressure_id = self.pressure.id
            _pressure_units = self.pressure.units
            _level = Axis(
                'lev_{0}'.format(self.id), _pressure[0, :],
                name='air_pressure_for_{0}'.format(self.id), units='Pa')

            data = np.zeros((ntin, nlevout), dtype=np.float64)
            for it in range(0, ntin):
                ff = interpolate.interp1d(
                    self.pressure.data[it, :], self.data[it, :],
                    bounds_error=False,
                    fill_value=(self.data[it, -1], self.data[it, 0])
                )  # Pad below and above, if necessary
                data[it, :] = ff(_pressure[it, :])

        else:

            logger.error(
                'Case unexpected for vertical interpolation of variable '
                '{0}'.format(self.id))
            raise ValueError(
                'Case unexpected for vertical interpolation of variable '
                '{0}'.format(self.id))

        return Variable(
            self.id, data=data, name=self.name, units=self.units,
            level=_level, time=self.time,
            height=_height, height_id=_height_id, height_units=_height_units,
            pressure=_pressure, pressure_id=_pressure_id,
            pressure_units=_pressure_units)

    def extend_vert(self, height=None, pressure=None, data=None,
                     time=None, tunits=None):
        """Extend the variable in the vertical with additional levels.

        Adds levels below and/or above the current vertical range,
        using `data`/`height` (or `data`/`pressure`) as the source of
        the extra levels to insert.

        Parameters
        ----------
        height : array_like, optional
            Height values associated with `data`, used to extend the
            variable vertically. Privileged over `pressure` when the
            variable has a height companion variable. May be 1D
            (shared across time) or 2D ``(time, level)``, in which
            case `time` and `tunits` must also be given.
        pressure : array_like, optional
            Pressure values associated with `data`. Currently not
            implemented.
        data : array_like
            Data values to insert at the levels given by `height`
            (or `pressure`).
        time : array_like, optional
            Time values associated with a 2D `height`/`pressure`
            array.
        tunits : str, optional
            Units of `time`, associated with a 2D `height`/`pressure`
            array.

        Returns
        -------
        Variable
            A new :class:`Variable` with the vertical range extended.

        Raises
        ------
        ValueError
            If the variable has no level axis, if both `height` and
            `pressure` are ``None``, if `data` is ``None``, or if a
            2D `height`/`pressure` is given without `time`/`tunits`.
        NotImplementedError
            For cases not yet handled (pressure-based extension,
            non-constant number of levels to add, or simultaneous
            extension below and above).
        """
        if self.level is None:
            logger.error(
                'Vertical extension requested for variable {0}, which '
                'does not have a level axis'.format(self.id))
            raise ValueError(
                'Vertical extension requested for variable {0}, which '
                'does not have a level axis'.format(self.id))

        if height is None and pressure is None:
            logger.error("height and pressure are both None.")
            raise ValueError("height and pressure are both None.")

        if data is None:
            logger.error("data is None")
            raise ValueError("data is None")

        if height is not None and len(height.shape) == 2 and \
                time is None and tunits is None:
            logger.error(
                'As given height is a 2D array, time must be given with '
                'tunits')
            raise ValueError(
                'As given height is a 2D array, time must be given with '
                'tunits')

        # if height is not None and height.shape != data.shape:
        #     logger.error('height and input data must have the same
        #     shape: {0} vs {1}'.format(height.shape, data.shape))
        #     raise ValueError('height and input data must have the same
        #     shape: {0} vs {1}'.format(height.shape, data.shape))

        if pressure is not None and len(pressure.shape) == 2 and \
                time is None and tunits is None:
            logger.error(
                'As given pressure is a 2D array, time must be given '
                'with tunits')
            raise ValueError(
                'As given pressure is a 2D array, time must be given '
                'with tunits')

        # if pressure is not None and pressure.shape != data.shape:
        #     logger.error('pressure and input data must have the same
        #     shape')
        #     raise ValueError('pressure and input data must have the
        #     same shape')

        ntin, nlevin = self.data.shape

        _height = None
        _height_id = None
        _height_units = None

        _pressure = None
        _pressure_id = None
        _pressure_units = None

        if height is not None and self.height is not None:

            hmax = np.max(self.height.data)
            hmin = np.min(self.height.data)
            # print(hmax)
            # print(data.shape)
            # print(height.shape)

            if len(height.shape) == 2:
                __time = Axis('tmp', time, name='tmp', units=tunits)
                # __time.info()
                __level = Axis('tmp', height[0, :], name='tmp', units='m')
                vartmp = Variable(
                    'tmp', data=data, name='tmp', units='-',
                    level=__level, time=__time,
                    height=height, height_id='height', height_units='m')
                var2add = vartmp.interpol_time(time=self.time)
            elif len(height.shape) == 1:
                __height = np.tile(height, (ntin, 1))
                __level = Axis('tmp', height, name='tmp', units='m')
                __data = np.tile(data, (ntin, 1))
                var2add = Variable(
                    'tmp', data=__data, name='tmp', units='-',
                    level=__level, time=self.time,
                    height=__height, height_id='height', height_units='m')
            else:
                logger.error(
                    "Shape of given height array is unexpected: "
                    "{0}".format(height.shape))
                raise ValueError(
                    "Shape of given height array is unexpected: "
                    "{0}".format(height.shape))

            # var2add.info()
            mask_up = var2add.height.data > hmax
            mask_dn = var2add.height.data < hmin

            nlev2add_up = np.sum(mask_up, axis=1)
            nlev2add_dn = np.sum(mask_dn, axis=1)
            if np.min(nlev2add_up) != np.max(nlev2add_up):
                logger.error(
                    "Case unexpected: the number of level to add is not "
                    "constant in time: min={0} max={1}".format(
                        np.min(nlev2add_up), np.max(nlev2add_up)))
                raise NotImplementedError

            if np.min(nlev2add_dn) != np.max(nlev2add_dn):
                logger.error(
                    "Case unexpected: the number of level to add is not "
                    "constant in time: min={0} max={1}".format(
                        np.min(nlev2add_dn), np.max(nlev2add_dn)))
                raise NotImplementedError

            if np.max(nlev2add_dn) > 0 and np.max(nlev2add_up) > 0:
                logger.error(
                    "Case unexpected: cannot add levels up and down at "
                    "the same time")
                raise NotImplementedError

            nlev2add_up = int(nlev2add_up[0])
            nlev2add_dn = int(nlev2add_dn[0])

            nlev_new = nlevin + nlev2add_up + nlev2add_dn

            _data = np.zeros((ntin, nlev_new), dtype=np.float64)
            _data[:, nlev2add_dn:(nlev2add_dn + nlevin)] = self.data[:, :]
            _data[:, (nlev2add_dn + nlevin):] = var2add.data[:, mask_up[0]]
            _data[:, :nlev2add_dn] = var2add.data[:, mask_dn[0]]

            _height = np.zeros((ntin, nlev_new), dtype=np.float64)
            _height[:, nlev2add_dn:(nlev2add_dn + nlevin)] = \
                self.height.data[:, :]
            _height[:, (nlev2add_dn + nlevin):] = \
                var2add.height.data[:, mask_up[0]]
            _height[:, :nlev2add_dn] = var2add.height.data[:, mask_dn[0]]
            _height_id = self.height.id
            _height_units = self.height.units

            _level = Axis(
                'lev_{0}'.format(self.id), _height[0, :],
                name='height_for_{0}'.format(self.id),
                units=self.height.units)

        elif pressure is not None and self.pressure is not None:

            logger.error('Case not yet coded')
            raise NotImplementedError

        else:

            logger.error(
                'Case unexpected for vertical extension of variable '
                '{0}'.format(self.id))
            raise ValueError(
                'Case unexpected for vertical extension of variable '
                '{0}'.format(self.id))

        return Variable(
            self.id, data=_data, name=self.name, units=self.units,
            level=_level, time=self.time,
            height=_height, height_id=_height_id, height_units=_height_units,
            pressure=_pressure, pressure_id=_pressure_id,
            pressure_units=_pressure_units,
            plotcoef=self.plotcoef, plotunits=self.plotunits)


def read(name, filein):
    """Read a :class:`Variable` from an open netCDF file.

    Reads the variable `name` from `filein`, together with its
    dimensions (as :class:`Axis` objects) and, if present, its
    height or pressure vertical-coordinate companion variable
    (identified from the netCDF ``coordinates`` attribute).

    Parameters
    ----------
    name : str
        Name of the netCDF variable to read.
    filein : netCDF4.Dataset
        Open netCDF file/dataset to read from.

    Returns
    -------
    Variable
        The variable read from the file, with its axes and, if
        applicable, its height/pressure companion variable attached.
    """
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

        tmp_ax = Axis(ax, filein[ax][:], **kwargs)
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

    # else:
    #     print 'ERROR: vertical variable for {0} is unexpected'.format(name)
    #     raise ValueError

    try:
        varout = Variable(
            name, data=tmp[:], name=tmp.standard_name, units=tmp.units,
            height=height, height_id=height_id, height_units=height_units,
            pressure=pressure, pressure_id=pressure_id,
            pressure_units=pressure_units,
            axes=axes, axlist=axlist)
    except AttributeError:
        varout = Variable(
            name, data=tmp[:], units=tmp.units,
            height=height, height_id=height_id, height_units=height_units,
            pressure=pressure, pressure_id=pressure_id,
            pressure_units=pressure_units,
            axes=axes, axlist=axlist)
    except:
        raise

    return varout


def interpol(var, levout=None, timeout=None, log=False):
    """Interpolate a variable onto new level and/or time axes.

    Parameters
    ----------
    var : Variable
        Variable to interpolate.
    levout : Axis, optional
        Target level axis. If ``None``, the original level axis (and
        data) is kept.
    timeout : Axis, optional
        Target time axis. If ``None`` (or if the variable has a
        single time step), the original time axis (and data) is
        kept.
    log : bool, optional
        If ``True``, interpolate in log-space along the vertical
        (only used for the time-and-level branch and the level-only
        branch). Defaults to ``False``.

    Returns
    -------
    Variable
        A new :class:`Variable` interpolated onto `levout` and/or
        `timeout`, as applicable.

    Raises
    ------
    ValueError
        If the input variable has neither a time nor a level axis.

    Notes
    -----
    TODO: Extrapolation to be revisited.
    """
    # TODO: Extrapolation to be revisited

    if not (var.time is None) and not (var.level is None):

        ntin, nlevin = var.time.length, var.level.length
        if not (levout is None):
            nlevout = levout.length
            tmp = np.zeros((ntin, nlevout, 1, 1), dtype=np.float64)
            for it in range(0, ntin):
                if log:
                    ff = interpolate.interp1d(
                        var.level.data, np.log(var.data[it, :, 0, 0]),
                        bounds_error=False,
                        fill_value=np.log(var.data[it, -1, 0, 0]))
                    tmp[it, :, 0, 0] = np.exp(ff(levout.data))
                else:
                    ff = interpolate.interp1d(
                        var.level.data, var.data[it, :, 0, 0],
                        bounds_error=False,
                        fill_value=var.data[it, -1, 0, 0])
                    tmp[it, :, 0, 0] = ff(levout.data)
        else:
            tmp = var.data
            nlevout = nlevin
            levout = var.data.level

        if not (timeout is None) and ntin > 1:
            ntout = timeout.length
            tmp2 = np.zeros((ntout, nlevout, 1, 1), dtype=np.float64)
            for ilev in range(0, nlevout):
                ff = interpolate.interp1d(
                    var.time.data, tmp[:, ilev, 0, 0],
                    bounds_error=False, fill_value=tmp[-1, ilev, 0, 0])
                tmp2[:, ilev, 0, 0] = ff(timeout.data)
        else:
            tmp2 = tmp
            ntout = ntin
            timeout = var.time

        varout = Variable(
            var.id, name=var.name, units=var.units, data=tmp2,
            time=timeout, level=levout, lat=var.lat, lon=var.lon)

    elif not (var.time is None):

        ntin = var.time.length
        if not (timeout is None) and ntin > 1:
            ntout = timeout.length
            tmp = np.zeros((ntout, 1, 1), dtype=np.float64)
            ff = interpolate.interp1d(
                var.time.data, var.data[:, 0, 0],
                bounds_error=False, fill_value=var.data[-1, 0, 0])
            tmp[:, 0, 0] = ff(timeout.data)
        else:
            tmp = np.log(var.data)
            tmp = var.data
            ntout = ntin
            timeout = var.time

        varout = Variable(
            var.id, name=var.name, units=var.units, data=tmp,
            time=timeout, lat=var.lat, lon=var.lon)

    elif not (var.level is None):

        nlevin = var.level.length
        if not (levout is None):
            nlevout = levout.length
            tmp = np.zeros((nlevout, 1, 1), dtype=np.float64)
            if log:
                ff = interpolate.interp1d(
                    var.level.data, np.log(var.data[:, 0, 0]),
                    bounds_error=False, fill_value=var.data[-1, 0, 0])
                tmp[:, 0, 0] = np.exp(ff(levout.data))
            else:
                ff = interpolate.interp1d(
                    var.level.data, var.data[:, 0, 0],
                    bounds_error=False, fill_value=var.data[-1, 0, 0])
                tmp[:, 0, 0] = ff(levout.data)
        else:
            tmp = var.data
            nlevout = nlevin
            levout = data.level

        varout = Variable(
            var.id, name=var.name, units=var.units, data=tmp,
            level=levout, lat=var.lat, lon=var.lon)

    else:
        logger.error('Weird, time and level are None for {0}'.format(var))
        raise ValueError(
            'Weird, time and level are None for {0}'.format(var))

    return varout
