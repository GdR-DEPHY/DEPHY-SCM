Welcome to DEPHY-SCM!
=====================

`DEPHY-SCM <https://github.com/GdR-DEPHY/DEPHY-SCM>`_ is a 
community initiative of the French `GdR DEPHY
<https://web.lmd.jussieu.fr/DEPHY/dephy.html>`_ research group,
providing:

- a **standard netCDF format** for describing the initial state and
  forcing of Single-Column Model (SCM) and Large-Eddy Simulation
  (LES) case studies, so that the same input file can be used to 
  drive many different models in a consistent way;
- a **library of ready-to-use SCM cases** in that format of
  well-known boundary-layer and convection cases such as
  ``GABLS1``, ``GABLS4``, ``ARMCU``, ``BOMEX``, ``RICO``, ``DYNAMO``,
  ``SANDU``, and many others;
- a **Python-based toolbox**, the ``dephycf`` package, used both to build
  and maintain that case library and to help modeling groups
  produce their own DEPHY-compliant driver files from their own data.

The goal is to lower the cost of running the same case across
different models (LES, SCMs embedded in various climate/NWP models) and
thereby to ease the use of 1D/LES simulations for model development and 
tuning and facilitate multi-model intercomparison studies. The initiative
builds also on earlier community efforts such as GCSS/GASS.

The DEPHY format standard itself — the precise netCDF conventions
(dimensions, required/optional variables, naming, units, global
attributes) — is described in the project's *Common Format*
specification document, versioned independently of the code (see the
PDF at the root of the repository). This toolbox implements and
helps use that standard; it does not redefine it.

.. toctree::
   :maxdepth: 2
   :caption: Contents

   install
   toolbox
   example
   cases/index
   api/index
