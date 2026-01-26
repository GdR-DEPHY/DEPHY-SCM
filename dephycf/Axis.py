#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Axis handling module.

This module defines the :class:`Axis` class, representing a one-dimensional
numeric axis typically used in scientific datasets (e.g., NetCDF files).
"""

import numpy as np

__all__ = ["Axis"]


class Axis:
    """
    Represents a numeric data axis, typically used in scientific datasets
    (e.g., NetCDF files).

    Parameters
    ----------
    axis_id : :class:`str`
        Short name of the axis identifier (e.g., 'time', 'lat', 'lon').
    data : array-like of :class:`float`
        Sequence of numeric values for the axis. Will be converted to
        :class:`numpy.ndarray` with dtype ``float64``.
    name : :class:`str`, optional
        Human-readable name for the axis.
    units : :class:`str`, optional
        Units associated with the axis values.
    \**kwargs : :class:`dict`, optional
        Additional arbitrary metadata. Each key/value pair is added as
        an attribute of the object and written to NetCDF as a variable
        attribute. For example::

            Axis("height", [0, 10, 20], name="altitude", units="m", source="radar")

        will create an Axis with a custom attribute ``source="radar"``.

    Raises
    ------
    ValueError
        If `data` is empty, None, or cannot be converted to a numeric NumPy array.
    """

    def __init__(self, axis_id, data, name=None, units=None, **kwargs):
        """Initialize an Axis object."""
        if data is None:
            raise ValueError("`data` cannot be None.")

        self.id = axis_id
        self.data = np.array(data, dtype=np.float64)

        if self.data.size == 0:
            raise ValueError("`data` cannot be empty.")
        if self.data.ndim != 1:
            raise ValueError("`data` must be one-dimensional.")

        self.length = self.data.shape[0]
        self.units = units
        self.name = name

        # Store extra attributes from keyword arguments
        for key, value in kwargs.items():
            setattr(self, key, value)

    def __repr__(self):
        """Return a concise string representation of the axis."""
        return (
            f"Axis(id='{self.id}', length={self.length}, "
            f"name='{self.name}', units='{self.units}')"
        )

    def info(self, data=False):
        """
        Print metadata information about the axis.

        Parameters
        ----------
        data : :class:`bool`, optional
            If True, also prints the full array of axis values.
            Default is False.

        Returns
        -------
        None
            Prints information to stdout.            
        """
        print("-" * 10, "axis id:", self.id)
        print("-" * 10, "name:", self.name)
        print("-" * 10, "units:", self.units)
        print("-" * 10, "length:", self.length)

        # Display additional attributes (if any)
        for key, value in self.__dict__.items():
            if key not in ["id", "name", "units", "length", "data"]:
                print("-" * 10, f"{key}: {value}")

        # Optionally display axis data
        if data:
            print("-" * 10, "data:", self.data)

    def write(self, file_obj):
        """
        Write this axis as a dimension and variable to a NetCDF file.

        Parameters
        ----------
        file_obj : :class:`netCDF4.Dataset`
            A NetCDF file opened in write or append mode.

        Raises
        ------
        AttributeError
            If `file_obj` does not support the required NetCDF API methods.

        Returns
        -------
        None
            Writes data and metadata into the provided NetCDF file object.
        """
        if self.id not in file_obj.dimensions:
            # Create dimension
            file_obj.createDimension(self.id, self.length)

            # Create variable for this dimension
            axis_var = file_obj.createVariable(self.id, np.float64, (self.id,))
            axis_var[:] = self.data

            # Add standard NetCDF metadata
            if self.name is not None:
                axis_var.standard_name = self.name
            if self.units is not None:
                axis_var.units = self.units

            # Add custom metadata attributes
            for key, value in self.__dict__.items():
                if key not in ["id", "units", "name", "data", "length"]:
                    axis_var.setncattr(key, value)
