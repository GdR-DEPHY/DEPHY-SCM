Attributes catalog
==================

This page documents the catalog of physical variables defined in
:mod:`dephycf.variables_attributes`.

Overview
--------

The variable metadata are stored in the ordered dictionary:

- :data:`dephycf.variables_attributes.attributes`

Each variable is identified by a short name used as its netCDF identifier
and associated with physical metadata such as units and CF standard names.

Metadata fields
---------------

Each entry in ``variables_attributes`` may contain:

================= =========================================
Key               Meaning
================= =========================================
``name``           CF standard name or descriptive name
``units``          Physical units
``plotcoef``       Scaling factor applied before plotting
``plotunits``      Units used for plotting
================= =========================================

Variables
---------

.. autofunction:: dephycf.variables_attributes.iter_attributes

Raw definition
--------------

For completeness, the full Python definition of the catalog is provided below.

.. literalinclude:: ../../dephycf/variables_attributes.py
   :language: python
   :start-after: attributes = OrderedDict([
   :end-before: ])
