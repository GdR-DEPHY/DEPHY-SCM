from collections import OrderedDict

attributes = OrderedDict([
        ###########################
        # Initial variables
        ###########################
        ('ps',       {'name': 'Surface pressure',                   'units': 'Pa',      'plotcoef': 0.01,  'plotunits': 'hPa'}),
        ('height',   {'name': 'Height above ground',                'units': 'm'}),
        ('pressure', {'name': 'Pressure',                           'units': 'Pa',      'plotcoef': 0.01,  'plotunits': 'hPa'}),
        ('u',        {'name': 'Zonal wind',                         'units': 'm s-1'}),
        ('v',        {'name': 'Meridional wind',                    'units': 'm s-1'}),
        ('temp',     {'name': 'Temperature',                        'units': 'K'}),
        ('theta',    {'name': 'Potential temperature',              'units': 'K'}),
        ('thetav',   {'name': 'Virtual potential temperature',      'units': 'K'}),
        ('thetal',   {'name': 'Liquid potential temperature',       'units': 'K'}),
        ('thetas',   {'name': 'Conservative potential temperature', 'units': 'K'}),
        ('qv',       {'name': 'Specific humidity',                  'units': 'kg kg-1', 'plotcoef': 1000., 'plotunits': 'g kg-1'}),
        ('qt',       {'name': 'Total water content',                'units': 'kg kg-1', 'plotcoef': 1000., 'plotunits': 'g kg-1'}),
        ('rv',       {'name': 'Water vapor mixing ratio',           'units': 'kg kg-1', 'plotcoef': 1000., 'plotunits': 'g kg-1'}),
        ('rt',       {'name': 'Total Water mixing ratio',           'units': 'kg kg-1', 'plotcoef': 1000., 'plotunits': 'g kg-1'}),
        ('ql',       {'name': 'Liquid water content',               'units': 'kg kg-1', 'plotcoef': 1000., 'plotunits': 'g kg-1'}),
        ('qi',       {'name': 'Ice water content',                  'units': 'kg kg-1', 'plotcoef': 1000., 'plotunits': 'g kg-1'}),
        ('tke',      {'name': 'Turbulent kinetic energy',           'units': 'm2 s-2'}),
        ('rh',       {'name': 'Relative humidity',                  'units': '%'}),
        ###########################
        # Forcing
        ###########################
        ('ps_forc',       {'name': 'Surface pressure for forcing', 'units': 'Pa', 'plotcoef': 0.01,  'plotunits': 'hPa'}),
        ('height_forc',   {'name': 'Height for forcing',  'units': 'm'}),
        ('pressure_forc', {'name': 'Pressure for forcing','units': 'Pa', 'plotcoef': 0.01,  'plotunits': 'hPa'}),
        # Geostrophic forcings
        ('ug', {'name': 'Geostrophic zonal wind',      'units': 'm s-1'}),
        ('vg', {'name': 'Geostrophic meridional wind', 'units': 'm s-1'}),
        # Advective forcings
        ('tadv',   {'name': 'Temperature large-scale advection',              'units': 'K s-1',       'plotcoef': 86400.,       'plotunits': 'K day-1'}),
        ('thadv',  {'name': 'Potential large-scale advection',                'units': 'K s-1',       'plotcoef': 86400.,       'plotunits': 'K day-1'}),
        ('thladv', {'name': 'Liquid potential large-scale advection',         'units': 'K s-1',       'plotcoef': 86400.,       'plotunits': 'K day-1'}),
        ('qvadv',  {'name': 'Specific humidity large-scale advection',        'units': 'kg kg-1 s-1', 'plotcoef': 86400.*1000., 'plotunits': 'g kg-1 day-1'}),
        ('qtadv',  {'name': 'Total water content large-scale advection',      'units': 'kg kg-1 s-1', 'plotcoef': 86400.*1000., 'plotunits': 'g kg-1 day-1'}),
        ('rvadv',  {'name': 'Water vapor mixing ratio large-scale advection', 'units': 'kg kg-1 s-1', 'plotcoef': 86400.*1000., 'plotunits': 'g kg-1 day-1'}),
        ('rtadv',  {'name': 'Total water mixing ratio large-scale advection', 'units': 'kg kg-1 s-1', 'plotcoef': 86400.*1000., 'plotunits': 'g kg-1 day-1'}),
        ('uadv',   {'name': 'Zonal wind large-scale advection',               'units': 'm s-2',       'plotcoef': 86400.,       'plotunits': 'K day-1'}),
        ('vadv',   {'name': 'Meridional wind large-scale advection',          'units': 'm s-2',       'plotcoef': 86400.,       'plotunits': 'K day-1'}),
        # Radiative tendency
        ('trad',  {'name': 'Radiative temperature tendency',            'units': 'K s-1',       'plotcoef': 86400.,       'plotunits': 'K day-1'}),
        ('thrad',  {'name': 'Radiative potential temperature tendency', 'units': 'K s-1',       'plotcoef': 86400.,       'plotunits': 'K day-1'}),
        # Vertical velocity forcing
        ('w',     {'name': 'Vertical velocity',          'units': 'm s-1', 'plotcoef': 100.,  'plotunits': 'cm s-1'}),
        ('omega', {'name': 'Vertical pressure velocity', 'units': 'Pa s-1'}),
        # Nudging profiles
        ('temp_nudging', {'name': 'Temperature profile for nudging',            'units': 'K'}),
        ('theta_nudging', {'name': 'Potential temperature profile for nudging', 'units': 'K'}),
        ('qv_nudging', {'name': 'Specific humidity profile for nudging',        'units': 'kg kg-1'}),
        ('qt_nudging', {'name': 'Total water content profile for nudging',      'units': 'kg kg-1'}),
        ('rv_nudging', {'name': 'Water vapor mixing ratio profile for nudging', 'units': 'kg kg-1'}),
        ('rt_nudging', {'name': 'Total water mixing ratio profile for nudging', 'units': 'kg kg-1'}),
        # Surface forcing
        ('sfc_sens_flx', {'name': 'Surface sensible heat flux',               'units': 'W m-2'}),
        ('sfc_lat_flx',  {'name': 'Surface latent heat flux',                 'units': 'W m-2'}),
        ('wpthetap',     {'name': 'Surface flux of potential temperature',    'units': 'K m s-1'}),
        ('wpqvp',        {'name': 'Surface flux of specific humidity',        'units': 'kg kg-1 m s-1'}),
        ('wpqtp',        {'name': 'Surface flux of total water',              'units': 'kg kg-1 m s-1'}),
        ('wprvp',        {'name': 'Surface flux of water vapor mixing ratio', 'units': 'kg kg-1 m s-1'}),
        ('wprtp',        {'name': 'Surface flux of total water mixing ratio', 'units': 'kg kg-1 m s-1'}),
        ('ts',           {'name': 'Surface temperature',                      'units': 'K'}),
        ('ustar',        {'name': 'Surface friction velocity',                'units': 'm s-1'}),
        ])
