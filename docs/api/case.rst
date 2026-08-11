.. _case-section:

**************************
Case (:mod:`dephycf.Case`)
**************************

``Case`` Module
==================

.. currentmodule:: dephycf.Case

.. automodule:: dephycf.Case
   :no-members:
   :no-undoc-members:

Usage Example
"""""""""""""

.. code-block:: python

   from dephycf.Case import Case

   # Create a case
   case = Case('EXAMPLE/control',
               lat=45.0, lon=5.0,
               startDate='2020-01-01 00:00:00',
               endDate='2020-01-01 06:00:00')
   case.set_title('An example DEPHY case')

   # Add initial-state variables
   case.add_init_ps(101325.)
   case.add_init_temp(temp_profile, lev=level_data, levtype='altitude')
   case.add_init_wind(u=u_profile, v=v_profile,
                       ulev=level_data, vlev=level_data)

   # Complete the case and convert it to a fully SCM-ready case
   scm_case = case.convert2SCM(lev=target_levels, levtype='altitude')

   # Write it out to a netCDF file
   scm_case.write('EXAMPLE_control_SCM_driver.nc')

.. seealso::
   :mod:`Variable`, :mod:`Axis`, :mod:`thermo`, :mod:`constants`,
   :mod:`attributes`, :mod:`variables_attributes`

API Reference
"""""""""""""

.. autoclass:: Case
   :no-members:

Construction and global attributes
""""""""""""""""""""""""""""""""""

.. automethod:: Case.__init__
.. automethod:: Case.set_dates
.. automethod:: Case.set_latlon
.. automethod:: Case.set_attribute
.. automethod:: Case.set_title
.. automethod:: Case.set_comment
.. automethod:: Case.set_reference
.. automethod:: Case.set_author
.. automethod:: Case.set_modifications
.. automethod:: Case.set_script
.. automethod:: Case.set_case_type

Adding and removing variables
"""""""""""""""""""""""""""""

.. automethod:: Case.remove_variable
.. automethod:: Case.add_variable
.. automethod:: Case.add_init_variable
.. automethod:: Case.add_forcing_variable
.. automethod:: Case.add_latitude
.. automethod:: Case.add_longitude
.. automethod:: Case.add_orography

Convenience wrappers for common initial-state variables
"""""""""""""""""""""""""""""""""""""""""""""""""""""""

.. automethod:: Case.add_init_ps
.. automethod:: Case.add_init_ts
.. automethod:: Case.add_init_thetas
.. automethod:: Case.add_init_height
.. automethod:: Case.add_init_pressure
.. automethod:: Case.add_init_temp
.. automethod:: Case.add_init_theta
.. automethod:: Case.add_init_thetal
.. automethod:: Case.add_init_qv
.. automethod:: Case.add_init_qt
.. automethod:: Case.add_init_rv
.. automethod:: Case.add_init_rt
.. automethod:: Case.add_init_hur
.. automethod:: Case.add_init_wind
.. automethod:: Case.add_init_tke

Convenience wrappers for common forcing variables
"""""""""""""""""""""""""""""""""""""""""""""""""

.. automethod:: Case.add_z0
.. automethod:: Case.add_z0h
.. automethod:: Case.add_z0q
.. automethod:: Case.add_surface_pressure_forcing
.. automethod:: Case.add_pressure_forcing
.. automethod:: Case.add_height_forcing
.. automethod:: Case.add_geostrophic_wind
.. automethod:: Case.add_vertical_velocity
.. automethod:: Case.add_temp_advection
.. automethod:: Case.add_theta_advection
.. automethod:: Case.add_thetal_advection
.. automethod:: Case.add_temp_radiation_tendency
.. automethod:: Case.add_theta_radiation_tendency
.. automethod:: Case.add_thetal_radiation_tendency
.. automethod:: Case.add_qv_advection
.. automethod:: Case.add_wind_advection
.. automethod:: Case.add_qt_advection
.. automethod:: Case.add_rv_advection
.. automethod:: Case.add_rt_advection
.. automethod:: Case.add_nudging
.. automethod:: Case.add_wind_nudging
.. automethod:: Case.add_temp_nudging
.. automethod:: Case.add_theta_nudging
.. automethod:: Case.add_thetal_nudging
.. automethod:: Case.add_qv_nudging
.. automethod:: Case.add_qt_nudging
.. automethod:: Case.add_rv_nudging
.. automethod:: Case.add_rt_nudging
.. automethod:: Case.add_ozone
.. automethod:: Case.activate_radiation
.. automethod:: Case.deactivate_radiation
.. automethod:: Case.add_surface_temp
.. automethod:: Case.add_surface_skin_temp
.. automethod:: Case.add_rad_ts
.. automethod:: Case.add_forcing_ts
.. automethod:: Case.add_forcing_thetas
.. automethod:: Case.add_surface_fluxes
.. automethod:: Case.set_betaevap
.. automethod:: Case.deactivate_surface_evaporation

Case information, plotting, netCDF read/write
"""""""""""""""""""""""""""""""""""""""""""""

.. automethod:: Case.info
.. automethod:: Case.write
.. automethod:: Case.read
.. automethod:: Case.plot
.. automethod:: Case.plot_compare

Computing derived variables
"""""""""""""""""""""""""""

.. automethod:: Case.compute_theta
.. automethod:: Case.compute_thetal
.. automethod:: Case.compute_temp
.. automethod:: Case.compute_qv
.. automethod:: Case.compute_qt
.. automethod:: Case.compute_rv
.. automethod:: Case.compute_rt
.. automethod:: Case.compute_hur
.. automethod:: Case.compute_tnta_adv
.. automethod:: Case.compute_tntheta_adv
.. automethod:: Case.compute_tnthetal_adv
.. automethod:: Case.compute_tnta_rad
.. automethod:: Case.compute_tntheta_rad
.. automethod:: Case.compute_tnthetal_rad
.. automethod:: Case.compute_tnqv_adv
.. automethod:: Case.compute_tnqt_adv
.. automethod:: Case.compute_tnrv_adv
.. automethod:: Case.compute_tnrt_adv
.. automethod:: Case.compute_ta_nud
.. automethod:: Case.compute_theta_nud
.. automethod:: Case.compute_thetal_nud
.. automethod:: Case.compute_qv_nud
.. automethod:: Case.compute_qt_nud
.. automethod:: Case.compute_rv_nud
.. automethod:: Case.compute_rt_nud

Interpolation and SCM conversion
""""""""""""""""""""""""""""""""

.. automethod:: Case.interpolate

.. note::
   The :meth:`interpolate` method contains placeholder
   ``usetemp``/``usetheta``/``usethetal`` parameters that are
   currently unused (see the source code).

.. automethod:: Case.add_missing_init_variables
.. automethod:: Case.add_missing_forcing_variables
.. automethod:: Case.convert2SCM

Vertically extending variables
""""""""""""""""""""""""""""""

.. automethod:: Case.extend_variable
.. automethod:: Case.extend_init_pressure
.. automethod:: Case.extend_init_wind
.. automethod:: Case.extend_init_temp
.. automethod:: Case.extend_init_theta
.. automethod:: Case.extend_init_thetal
.. automethod:: Case.extend_init_qv
.. automethod:: Case.extend_init_qt
.. automethod:: Case.extend_init_rv
.. automethod:: Case.extend_init_rt
.. automethod:: Case.extend_init_hur
.. automethod:: Case.extend_geostrophic_wind
.. automethod:: Case.extend_vertical_velocity
.. automethod:: Case.extend_temperature_advection
.. automethod:: Case.extend_theta_advection
.. automethod:: Case.extend_thetal_advection
.. automethod:: Case.extend_qv_advection
.. automethod:: Case.extend_qt_advection
.. automethod:: Case.extend_rv_advection
.. automethod:: Case.extend_wind_advection
.. automethod:: Case.extend_rt_advection
