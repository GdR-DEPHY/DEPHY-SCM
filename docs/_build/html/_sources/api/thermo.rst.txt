*********************************************************************
Thermodyanmical functions and other utilities (:mod:`dephycf.thermo`)
*********************************************************************

Example
"""""""

You can import and use the :mod:`dephycf.thermo` module directly:

    >>> from dephycf import thermo
    
You may start by initializing a vertical coordinate

    >>> z = np.arange(0.0, 12001.0, 500.0)  # m
    
and then using it to further initialize a temperature vertical profile starting à 300 K at the surface and decreasing with a lapse rate of 6.5 K km-1

    >>> temp = thermo._init_temp(z)
    >>> print(temp)
    [300.   296.75 293.5  290.25 287.   283.75 280.5  277.25 274.   270.75
     267.5  264.25 261.   257.75 254.5  251.25 248.   244.75 241.5  238.25
     235.   231.75 228.5  225.25 222.  ]

You may also initialize a vertical profile of specific humidity:

    >>> qv = thermo._init_qv(z)
    >>> print(qv)
    [1.00000000e-02 7.78800783e-03 6.06530660e-03 4.72366553e-03
     3.67879441e-03 2.86504797e-03 2.23130160e-03 1.73773943e-03
     1.35335283e-03 1.05399225e-03 8.20849986e-04 6.39278612e-04
     4.97870684e-04 3.87742078e-04 3.01973834e-04 2.35177459e-04
     1.83156389e-04 1.42642339e-04 1.11089965e-04 8.65169520e-05
     6.73794700e-05 5.24751840e-05 4.08677144e-05 3.18278080e-05
     2.47875218e-05]

Introducing a surface pressure, you can then compute the pressure levels:

    >>> ps = 100000.0 # Pa
    >>> pressure = thermo.z2p(z, ps, ta=temp, qv=qv)
    >>> print(pressure)
    [100000.          94464.9259531   89174.17250373  84121.18329611
     79299.02674886  74700.54043306  70318.43955459  66145.39805999
     62174.10895315  58397.32887947  54807.9108365   51398.82794401
     48163.19049423  45094.2579596   42185.44722181  39430.33797296
     36822.67600293  34356.37490862  32025.51662654  29824.35108947
     27747.29523191  25788.93151266  23944.00607998  22207.42667342
     20574.26033205]

It can be used to compute the potential temperature profile:

    >>> theta = thermo.t2theta(pressure, temp)
    >>> print(theta)
    [300.         301.61731794 303.26724712 304.94943714 306.66379723
     308.41045089 310.18970005 312.00199686 313.84792172 315.72816636
     317.64352085 319.59486402 321.58315636 323.60943516 325.67481133
     327.78046767 329.92765834 332.1177093  334.35201971 336.63206395
     338.95939443 341.335645   343.76253495 346.24187358 348.77556536]

Here follows the thermodynmical functions and utilities defined in :mod:`dephycf.thermo`.

API Reference
"""""""""""""

.. automodule:: dephycf.thermo
   :members:
   :undoc-members:
   :show-inheritance:
   :member-order: bysource
