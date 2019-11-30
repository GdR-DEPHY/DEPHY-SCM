from collections import OrderedDict

attributes = OrderedDict([
        ###########################
        # Initial variables
        ###########################
        ('ps',       {'name': 'Surface pressure',             'units': 'Pa',      'plotcoef': 0.01,  'plotunits': 'hPa'}),
        ('height',   {'name': 'Height above ground',          'units': 'm'}),
        ('pressure', {'name': 'Pressure',                     'units': 'Pa',      'plotcoef': 0.01,  'plotunits': 'hPa'}),
        ('u',        {'name': 'Zonal wind',                   'units': 'm s-1'}),
        ('v',        {'name': 'Meridional wind',              'units': 'm s-1'}),
        ('temp',     {'name': 'Temperature',                  'units': 'K'}),
        ('theta',    {'name': 'Potential temperature',        'units': 'K'}),
        ('thetal',   {'name': 'Liquid potential temperature', 'units': 'K'}),
        ('qv',       {'name': 'Specific humidity',            'units': 'kg kg-1', 'plotcoef': 1000., 'plotunits': 'g kg-1'}),
        ('qt',       {'name': 'Total water content',          'units': 'kg kg-1', 'plotcoef': 1000., 'plotunits': 'g kg-1'}),
        ('rv',       {'name': 'Water vapor mixing ratio',     'units': 'kg kg-1', 'plotcoef': 1000., 'plotunits': 'g kg-1'}),
        ('rt',       {'name': 'Total Water mixing ratio',     'units': 'kg kg-1', 'plotcoef': 1000., 'plotunits': 'g kg-1'}),
        ('ql',       {'name': 'Liquid water content',         'units': 'kg kg-1', 'plotcoef': 1000., 'plotunits': 'g kg-1'}),
        ('qi',       {'name': 'Ice water content',            'units': 'kg kg-1', 'plotcoef': 1000., 'plotunits': 'g kg-1'}),
        ('tke',      {'name': 'Turbulent kinetic energy',     'units': 'm2 s-2'}),
        ###########################
        # Forcing
        ###########################
        ('ps_forc',       {'name': 'Forcing surface pressure', 'units': 'Pa', 'plotcoef': 0.01,  'plotunits': 'hPa'}),
        ('height_forc',   {'name': 'Forcing height',           'units': 'm'}),
        ('pressure_forc', {'name': 'Forcing pressure',         'units': 'Pa', 'plotcoef': 0.01,  'plotunits': 'hPa'}),
        # Geostrophic forcings
        ('ug', {'name': 'Geostrophic zonal wind',      'units': 'm s-1'}),
        ('vg', {'name': 'Geostrophic meridional wind', 'units': 'm s-1'}),
        # Advective forcings
        ('tadv',  {'name': 'Large-scale advection of temperature',              'units': 'K s-1',       'plotcoef': 86400.,       'plotunits': 'K day-1'}),
        ('thadv', {'name': 'Large-scale advection of potential temperature',    'units': 'K s-1',       'plotcoef': 86400.,       'plotunits': 'K day-1'}),
        ('qvadv', {'name': 'Large-scale advection of specific humidity',        'units': 'kg kg-1 s-1', 'plotcoef': 86400.*1000., 'plotunits': 'g kg-1 day-1'}),
        ('qtadv', {'name': 'Large-scale advection of total water content',      'units': 'kg kg-1 s-1', 'plotcoef': 86400.*1000., 'plotunits': 'g kg-1 day-1'}),
        ('rvadv', {'name': 'Large-scale advection of water vapor mixing ratio', 'units': 'kg kg-1 s-1', 'plotcoef': 86400.*1000., 'plotunits': 'g kg-1 day-1'}),
        ('rtadv', {'name': 'Large-scale advection of totat water mixing ratio', 'units': 'kg kg-1 s-1', 'plotcoef': 86400.*1000., 'plotunits': 'g kg-1 day-1'}),
        # Vertical velocity forcing
        ('w',     {'name': 'Vertical velocity',          'units': 'm s-1', 'plotcoef': 100.,  'plotunits': 'cm s-1'}),
        ('omega', {'name': 'Pressure vertical velocity', 'units': 'Pa s-1'}),
        # Surface forcing
        ('sfc_sens_flx', {'name': 'Surface sensible heat flux', 'units': 'W m-2'}),
        ('sfc_lat_flx',  {'name': 'Surface latent heat flux',   'units': 'W m-2'}),
        ('ts',           {'name': 'Surface temperature',        'units': 'K'}),
        ])
