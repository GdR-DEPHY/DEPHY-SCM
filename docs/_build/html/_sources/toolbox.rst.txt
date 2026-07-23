``dephycf`` toolbox
-----------------------------------

``dephycf`` is a pure-Python package built around a small set of
Python objects:

:class:`~dephycf.Axis.Axis`
    A labelled coordinate axis (time, level...), with its data,
    units and metadata.

:class:`~dephycf.Variable.Variable`
    A data array associated with its axes and, where relevant, a
    height or pressure companion variable used as an alternative
    vertical coordinate. Provides netCDF read/write, plotting, and
    time/vertical interpolation.

:class:`~dephycf.Case.Case`
    The main entry point: a full case (initial state + forcing),
    with global metadata (dates, location, surface type...). It
    exposes convenience methods to add each standard DEPHY variable
    (initial profiles, forcings, nudging, radiative tendencies...),
    to derive missing-but-required variables from whatever is
    available (:meth:`~dephycf.Case.Case.add_missing_init_variables`,
    :meth:`~dephycf.Case.Case.add_missing_forcing_variables`), to
    interpolate a case onto a new time/vertical grid
    (:meth:`~dephycf.Case.Case.interpolate`), and to convert a case
    end-to-end into an SCM-ready driver file
    (:meth:`~dephycf.Case.Case.convert2SCM`,
    :meth:`~dephycf.Case.Case.write`).

Supporting modules cover unit conversions and thermodynamic
relations (``thermo``), physical constants (``constants``), known
variable/attribute names (``variables_attributes``, ``attributes``),
and basic plotting helpers (``plotbasics``).

A typical workflow looks like:

.. code-block:: python

   from dephycf.Case import Case

   case = Case('MYCASE/REF', lat=..., lon=...,
               startDate=..., endDate=...)

   case.add_init_ps(...)
   case.add_init_temp(..., lev=..., levtype='pressure')
   case.add_init_wind(u=..., v=..., ulev=..., vlev=...)
   # ... add forcing variables ...

   case.add_missing_init_variables()
   case.add_missing_forcing_variables()

   scm_case = case.convert2SCM(lev=target_levels, levtype='pressure')
   scm_case.write('MYCASE_REF_SCM_driver.nc')
