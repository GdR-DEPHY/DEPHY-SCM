.. _tutorial-section:

Tutorial 
========

From a Case Description to an SCM-Ready Driver
----------------------------------------------

This tutorial walks through a complete, real-world example: the
``ARMCU/REF`` case (ARM-Cumulus). It illustrates the two-stage
workflow used throughout the DEPHY-SCM case library:

1. **"DEF" driver** — encode the case *as originally described* in
   the literature (`Brown et al. (2002) <https://doi.org/10.1256/003590002320373210>`_,
   `Lenderink et al. (2004) <https://doi.org/10.1256/qj.03.122>`_ with :class:`~dephycf.Case.Case`,
   and write it out to a first netCDF file.
2. **"SCM" driver** — read that file back, complete and extend it, and
   interpolate it onto regular or finer time and vertical grids, thereby producing the final driver file used to
   actually run single-column model simulations.

The full scripts are extensively reproduced in available in the :ref:`ARMCU/REF <ARMCU_REF-section>` case section
— this page reproduces and annotates them.

.. note::
   The ``:lines:`` ranges used by the ``literalinclude`` directives
   below refer to the scripts as they exist today. If those scripts
   are edited, the line numbers here will need to be updated
   accordingly (or replaced with ``:start-after:``/``:end-before:``
   markers).

Stage 1 — Building the "DEF" driver
--------------------------------------

Setup
~~~~~~~

The script starts with the usual imports, plus two flags controlling
whether to plot the case and print debugging information:

.. literalinclude:: ../ARMCU/REF/driver_DEF.py
   :language: python
   :lines: 15-27

Creating the case and its metadata
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A :class:`~dephycf.Case.Case` is created with its identifier 
(note the ``CASE/SUBCASE`` convention), location,
start/end dates and surface type, then documented with a few global
attributes (:meth:`~dephycf.Case.Case.set_title`,
:meth:`~dephycf.Case.Case.set_reference`,
:meth:`~dephycf.Case.Case.set_author`,
:meth:`~dephycf.Case.Case.set_script`):

.. literalinclude:: ../ARMCU/REF/driver_DEF.py
   :language: python
   :lines: 33-44

Initial state
~~~~~~~~~~~~~~~

The initial surface pressure is added first
(:meth:`~dephycf.Case.Case.add_init_ps`), followed by the initial
vertical profiles of potential temperature, total water mixing ratio
and wind, given here on the same altitude grid
(:meth:`~dephycf.Case.Case.add_init_theta`,
:meth:`~dephycf.Case.Case.add_init_rt`,
:meth:`~dephycf.Case.Case.add_init_wind`):

.. literalinclude:: ../ARMCU/REF/driver_DEF.py
   :language: python
   :lines: 50-70

An initial profile of turbulent kinetic energy is defined and added
with :meth:`~dephycf.Case.Case.add_init_tke`:

.. literalinclude:: ../ARMCU/REF/driver_DEF.py
   :language: python
   :lines: 72-83

Forcing
~~~~~~~~~

ARMCU implements a geopstrophic wind forcing, using a time-constant geostrophic wind, equal to the initial wind profile. It is
prescribed with :meth:`~dephycf.Case.Case.add_geostrophic_wind`:

.. literalinclude:: ../ARMCU/REF/driver_DEF.py
   :language: python
   :lines: 89-90

The large-scale advection of potential temperature (together with
its radiative tendency, via ``include_rad=True``) and of total water
are built from a time/height table and added with
:meth:`~dephycf.Case.Case.add_theta_advection` and
:meth:`~dephycf.Case.Case.add_rt_advection`. Note the conversion to
SI units (K s-1, kg kg-1 s-1) expected by the DEPHY format:

.. literalinclude:: ../ARMCU/REF/driver_DEF.py
   :language: python
   :lines: 92-130

Finally, the surface sensible and latent heat fluxes are prescribed
directly (rather than a surface temperature), together with a
roughness length, via :meth:`~dephycf.Case.Case.add_surface_fluxes`.
Note the time indices should start from 0 as the time axis unit counts from the starting date of the case

.. literalinclude:: ../ARMCU/REF/driver_DEF.py
   :language: python
   :lines: 132-146

Writing and plotting
~~~~~~~~~~~~~~~~~~~~~~~

The case is written to a netCDF file with
the :meth:`~dephycf.Case.Case.write` method, optionally inspected with the
:meth:`~dephycf.Case.Case.info` method, and optionally plotted with the
:meth:`~dephycf.Case.Case.plot` method:

.. literalinclude:: ../ARMCU/REF/driver_DEF.py
   :language: python
   :lines: 152-162

This produces ``ARMCU_REF_DEF_driver.nc`` which implements the DEPHY standard for 
the originally described case (Note ``DEF`` which emphasizes this point).
This file is not necessarily SCM-ready: different time and vertical grids for 
the initial state and forcing levels, variables not necessarily useful 
as state variables commonly used in climate and NWP models... 
The SCM driver scripts enhanced the originally-defined case.


Stage 2 — Converting to an SCM-ready driver
-----------------------------------------------

Setup and reading the "DEF" case back
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The second script starts a *new* :class:`~dephycf.Case.Case` and
populates it by reading back the netCDF file written in stage 1,
with :meth:`~dephycf.Case.Case.read`:

.. literalinclude:: ../ARMCU/REF/driver_SCM.py
   :language: python
   :lines: 15-41

Extending and completing the profiles
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Before interpolation, the case is completed with information beyond
what the original description provided:

- the wind and total water profiles are vertically extended with
  :meth:`~dephycf.Case.Case.extend_init_wind` and
  :meth:`~dephycf.Case.Case.extend_init_rt` assuming constant value
  (vertical grids of NWP and climate models often extend far up);
- the temperature profile is extended above 12 km using an
  independent ERA5 reanalysis profile, via
  :meth:`~dephycf.Case.Case.extend_init_theta` (note the use of the
  :mod:`~dephycf.thermo` module to convert temperature to potential
  temperature);
- a surface skin temperature time series, from independent
  observations, is added with
  :meth:`~dephycf.Case.Case.add_surface_skin_temp`.

.. literalinclude:: ../ARMCU/REF/driver_SCM.py
   :language: python
   :lines: 48-73

Interpolating and converting
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A target vertical grid (regular, 10 m resolution) and a target time
grid (regular, 30-minute time step, covering the full simulation
period) are defined, then the whole case is interpolated onto them
and completed in a single call to
:meth:`~dephycf.Case.Case.convert2SCM`. Note this function also compute 
another sets of state variables, assuming a few hypotheses.
A couple of metadata attributes are updated on the resulting, new case:

.. literalinclude:: ../ARMCU/REF/driver_SCM.py
   :language: python
   :lines: 75-92

Writing and plotting
~~~~~~~~~~~~~~~~~~~~~~~

The SCM-ready case is written out with
the :meth:`~dephycf.Case.Case.write` method, plotted with
the :meth:`~dephycf.Case.Case.plot` method, and optionally compared to the
original "DEF" case with the :meth:`~dephycf.Case.Case.plot_compare` method —
a useful sanity check that the interpolation/extension steps did not
distort the original forcing:

.. literalinclude:: ../ARMCU/REF/driver_SCM.py
   :language: python
   :lines: 98-109

This produces ``ARMCU_REF_SCM_driver.nc``: a case implementation on a regular
time/height grid, with all the variables required to drive a
single-column model, ready to be used directly.

Summary
---------

.. list-table::
   :header-rows: 1

   * - Stage
     - Input
     - Key methods
     - Output
   * - DEF
     - Transcribed case description
     - ``add_init_*``, ``add_*_advection``, ``add_surface_fluxes``,
       ``write``
     - ``CASE_SUBCASE_DEF_driver.nc``
   * - SCM
     - The DEF driver file (+ optional auxiliary data)
     - ``read``, ``extend_init_*``, ``add_*``, ``convert2SCM``, ``write``
     - ``CASE_SUBCASE_SCM_driver.nc``

.. seealso::
   :ref:`Case section <case-section>` for the full :class:`~dephycf.Case.Case` API
   reference, and :ref:`Variable section <variable-section>` for the underlying
   :class:`~dephycf.Variable.Variable` objects.

