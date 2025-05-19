from collections import OrderedDict

attributes = OrderedDict([
        ###########################
        # Initial variables
        ###########################
        ('ps',       {'name': 'surface_air_pressure',                       'units': 'Pa', 'plotcoef': 0.01,  'plotunits': 'hPa'}),
        ('zh',       {'name': 'height_above_surface',                       'units': 'm',  'plotcoef': 0.001, 'plotunits': 'km'}),
        ('pa',       {'name': 'air_pressure',                               'units': 'Pa', 'plotcoef': 0.01,  'plotunits': 'hPa'}),
        ('ua',       {'name': 'eastward_wind',                              'units': 'm s-1'}),
        ('va',       {'name': 'northward_wind',                             'units': 'm s-1'}),
        ('ta',       {'name': 'air_temperature',                            'units': 'K'}),
        ('theta',    {'name': 'air_potential_temperature',                  'units': 'K'}),
        ('thetav',   {'name': 'air_virtual_potential_temperature',          'units': 'K'}),
        ('thetal',   {'name': 'air_liquid_potential_temperature',           'units': 'K'}),
        ('hur',      {'name': 'relative_humidity',                          'units': '1', 'plotcoef':  100., 'plotunits': '%'}),
        ('qv',       {'name': 'specific_humidity',                          'units': '1', 'plotcoef': 1000., 'plotunits': 'g kg-1'}),
        ('qt',       {'name': 'mass_fraction_of_water_in_air',              'units': '1', 'plotcoef': 1000., 'plotunits': 'g kg-1'}),
        ('ql',       {'name': 'mass_fraction_of_cloud_liquid_water_in_air', 'units': '1', 'plotcoef': 1000., 'plotunits': 'g kg-1'}),
        ('qi',       {'name': 'mass_fraction_of_cloud_ice_water_in_air',    'units': '1', 'plotcoef': 1000., 'plotunits': 'g kg-1'}),        
        ('rv',       {'name': 'humidity_mixing_ratio',                      'units': '1', 'plotcoef': 1000., 'plotunits': 'g kg-1'}),
        ('rt',       {'name': 'water_mixing_ratio',                         'units': '1', 'plotcoef': 1000., 'plotunits': 'g kg-1'}),
        ('rl',       {'name': 'cloud_liquid_water_mixing_ratio',            'units': '1', 'plotcoef': 1000., 'plotunits': 'g kg-1'}),
        ('ri',       {'name': 'cloud_ice_water_mixing_ratio',               'units': '1', 'plotcoef': 1000., 'plotunits': 'g kg-1'}),        
        ('tke',      {'name': 'specific_turbulent_kinetic_energy',          'units': 'm2 s-2'}),
        ('hur',      {'name': 'relative_humidity',                          'units': '%'}),
        # Surface
        ('ts',    {'name': 'surface_temperature',                  'units': 'K'}),
        ('thetas',{'name': 'surface_potential_temperature',        'units': 'K'}),
        ('mrsos', {'name': 'mass_content_of_water_in_soil_layer',  'units': 'kg m-2'}),
        ###########################
        # Forcing
        ###########################
        # General
        ('lat',     {'name': 'latitude',                     'units': 'degrees_north'}),
        ('lon',     {'name': 'longitude',                    'units': 'degrees_east'}),
        ('ps_forc', {'name': 'forcing_surface_air_pressure', 'units': 'Pa', 'plotcoef': 0.01,  'plotunits': 'hPa'}),
        ('zh_forc', {'name': 'forcing_height_above_surface', 'units': 'm',  'plotcoef': 0.001, 'plotunits': 'km'}),
        ('pa_forc', {'name': 'forcing_air_pressure',         'units': 'Pa', 'plotcoef': 0.01,  'plotunits': 'hPa'}),
        # Geostrophic forcings
        ('ug', {'name': 'geostrophic_eastward_wind',  'units': 'm s-1'}),
        ('vg', {'name': 'geostrophic_northward_wind', 'units': 'm s-1'}),
        # Advective forcings
        ('tnta_adv',     {'name': 'tendency_of_air_temperature_due_to_advection',                  'units': 'K s-1', 'plotcoef': 86400.,       'plotunits': 'K day-1'}),
        ('tntheta_adv',  {'name': 'tendency_of_air_potential_temperature_due_to_advection',        'units': 'K s-1', 'plotcoef': 86400.,       'plotunits': 'K day-1'}),
        ('tnthetal_adv', {'name': 'tendency_of_air_liquid_potential_temperature_due_to_advection', 'units': 'K s-1', 'plotcoef': 86400.,       'plotunits': 'K day-1'}),
        ('tnqv_adv',     {'name': 'tendency_of_specific_humidity_due_to_advection',                'units': 's-1',   'plotcoef': 86400.*1000., 'plotunits': 'g kg-1 day-1'}),
        ('tnqt_adv',     {'name': 'tendency_of_mass_fraction_of_water_in_air_due_To_advection',    'units': 's-1',   'plotcoef': 86400.*1000., 'plotunits': 'g kg-1 day-1'}),
        ('tnrv_adv',     {'name': 'tendency_of_humidity_mixing_ratio_due_to_advection',            'units': 's-1',   'plotcoef': 86400.*1000., 'plotunits': 'g kg-1 day-1'}),
        ('tnrt_adv',     {'name': 'tendency_of_water_mixing_ratio_due_to_advection',               'units': 's-1',   'plotcoef': 86400.*1000., 'plotunits': 'g kg-1 day-1'}),
        ('tnua_adv',     {'name': 'tendency_of_eastward_wind_due_to_advection',                    'units': 'm s-2', 'plotcoef': 86400.,       'plotunits': 'K day-1'}),
        ('tnva_adv',     {'name': 'tendency_of_northward_wind_due_to_advection',                   'units': 'm s-2', 'plotcoef': 86400.,       'plotunits': 'K day-1'}),
        # Radiative tendency
        ('tnta_rad',     {'name': 'tendency_of_air_temperature_due_to_radiative_heating',                  'units': 'K s-1', 'plotcoef': 86400., 'plotunits': 'K day-1'}),
        ('tntheta_rad',  {'name': 'tendency_of_air_potential_temperature_due_to_radiative_heating',        'units': 'K s-1', 'plotcoef': 86400., 'plotunits': 'K day-1'}),
        ('tnthetal_rad', {'name': 'tendency_of_air_liquid_potential_temperature_due_to_radiative_heating', 'units': 'K s-1', 'plotcoef': 86400., 'plotunits': 'K day-1'}),
        # Vertical velocity forcing
        ('wa',  {'name': 'upward_air_velocity',                 'units': 'm s-1', 'plotcoef': 100.,  'plotunits': 'cm s-1'}),
        ('wap', {'name': 'lagrangian_tendency_of_air_pressure', 'units': 'Pa s-1','plotcoef': 864.,  'plotunits': 'hPa day-1'}),
        # Nudging profiles  
        ('ua_nud',                     {'name': 'nudging_eastward_wind',                                    'units': 'm s-1'}),
        ('va_nud',                     {'name': 'nudging_northward_wind',                                   'units': 'm s-1'}),
        ('ta_nud',                     {'name': 'nudging_air_temperature',                                  'units': 'K'}),
        ('theta_nud',                  {'name': 'nudging_air_potential_temperature',                        'units': 'K'}),
        ('thetal_nud',                 {'name': 'nudging_air_liquid_potential_temperature',                 'units': 'K'}),
        ('qv_nud',                     {'name': 'nudging_specific_humidity',                                'units': '1', 'plotcoef': 1000., 'plotunits': 'g kg-1'}),
        ('qt_nud',                     {'name': 'nudging_mass_fraction_of_water_in_air',                    'units': '1', 'plotcoef': 1000., 'plotunits': 'g kg-1'}),
        ('rv_nud',                     {'name': 'nudging_humidity_mixing_ratio',                            'units': '1', 'plotcoef': 1000., 'plotunits': 'g kg-1'}),
        ('rt_nud',                     {'name': 'nudging_water_mixing_ratio',                               'units': '1', 'plotcoef': 1000., 'plotunits': 'g kg-1'}),
        ('nudging_coefficient_ua',     {'name': 'nudging_coefficient_for_eastward_wind',                    'units': 's-1'}),
        ('nudging_coefficient_va',     {'name': 'nudging_coefficient_for_northward_wind',                   'units': 's-1'}),
        ('nudging_coefficient_ta',     {'name': 'nudging_coefficient_for_air_temperature',                  'units': 's-1'}),
        ('nudging_coefficient_theta',  {'name': 'nudging_coefficient_for_air_potential_temperature',        'units': 's-1'}),
        ('nudging_coefficient_thetal', {'name': 'nudging_coefficient_for_air_liquid_potential_temperature', 'units': 's-1'}),
        ('nudging_coefficient_qv',     {'name': 'nudging_coefficient_for_specific_humidity',                'units': 's-1'}),    
        ('nudging_coefficient_qt',     {'name': 'nudging_coefficient_for_mass_fraction_of_water_in_air',    'units': 's-1'}),    
        ('nudging_coefficient_rv',     {'name': 'nudging_coefficient_for_humidity_mixing_ratio',            'units': 's-1'}),    
        ('nudging_coefficient_rt',     {'name': 'nudging_coefficient_for_water_mixing_ratio',               'units': 's-1'}),              
        # Surface forcing
        ('orog',       {'name': 'surface_altitude',                             'units': 'm'}),
        ('hfss',       {'name': 'surface_upward_sensible_heat_flux',            'units': 'W m-2'}),
        ('hfls',       {'name': 'surface_upward_latent_heat_flux',              'units': 'W m-2'}),
        ('wpthetap_s', {'name': 'surface_upward_potentail_temperature_flux',    'units': 'K m s-1'}),
        ('wpqvp_s',    {'name': 'surface_upward_specific_humidity_flux',        'units': 'm s-1'}),
        ('wpqtp_s',    {'name': 'surface_upward_water_mass_fraction_flux',      'units': 'm s-1'}),
        ('wprvp_s',    {'name': 'surface_upward_humidity_mixing_ratio_flux',    'units': 'm s-1'}),
        ('wprtp_s',    {'name': 'surface_upward_water_mixing_ratio_flux',       'units': 'm s-1'}),
        ('ts_forc',    {'name': 'forcing_surface_temperature',                  'units': 'K'}),
        ('thetas_forc',{'name': 'forcing_surface_potential_temperature',        'units': 'K'}),
        ('tskin',      {'name': 'surface_skin_temperature',                     'units': 'K'}),
        ('ustar',      {'name': 'surface_friction_velocity',                    'units': 'm s-1'}),
        ('z0',         {'name': 'surface_roughness_length_for_momentum_in_air', 'units': 'm'}),
        ('z0h',        {'name': 'surface_roughness_length_for_heat_in_air',     'units': 'm'}),
        ('z0q',        {'name': 'surface_roughness_length_for_humidity_in_air', 'units': 'm'}),
        ('beta',       {'name': 'soil_water_stress_factor',                     'units': '-'}),
        ('mrsos_forc', {'name': 'forcing_mass_content_of_water_in_soil_layer',  'units': 'kg m-2'}),
        # Atmospheric composition
        ('o3', {'name': 'mole_fraction_of_ozone_in_air', 'units': '1'}),
        # Radiation
        ('alb',  {'name': 'surface_albedo',              'units': '1'}),
        ('emis', {'name': 'surface_longwave_emissivity', 'units': '1'}),
        ('sza',  {'name': 'solar_zenith_angle',          'units': 'degree'}),
        ('i0',   {'name': 'solar_irradiance',            'units': 'W m-2'}),
        ])
