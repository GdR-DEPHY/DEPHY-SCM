.. _variable-section:

**********************************
Variable (:mod:`dephycf.Variable`)
**********************************

.. _variable:

This module provides the :class:`Variable` class, a container that
associates a data array with its axes (time, level, latitude,
longitude) and, optionally, with a companion vertical-coordinate
variable (height and/or pressure), itself stored as an instance of
:class:`Variable`.

It also provides two module-level functions:

- :func:`read`, to read a :class:`Variable` from an already open
  netCDF file;
- :func:`interpol`, to interpolate a :class:`Variable` in time and/or
  along the vertical.

Example
"""""""

.. code-block:: python

   import netCDF4 as nc
   from dephycf import Variable

   # Read a variable from a netCDF file
   with nc.Dataset('example.nc') as filein:
       var = Variable.read('ta', filein)
       var.info()

   # Time interpolation
   var_interp = var.interpol_time(time=new_time_axis)

Notes
"""""

.. note::
   The :func:`interpol` function contains a TODO regarding the
   revision of the extrapolation strategy (see the source code).

.. seealso::
   :mod:`Axis`, :mod:`plotbasics`, :mod:`variables_attributes`

API Reference
"""""""""""""

.. automodule:: dephycf.Variable
   :members:
   :undoc-members:
   :show-inheritance:
