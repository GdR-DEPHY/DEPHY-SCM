**************************
Axis (:mod:`dephycf.Axis`)
**************************

Example
"""""""

    >>> import numpy as np
    >>> from netCDF4 import Dataset
    >>> from axis import Axis
    >>> time_axis = Axis("time", np.arange(0, 24),
    ...     name="Time", units="hours since 2000-01-01 00:00:00")
    >>> time_axis.info()
    ---------- axis id: time
    ---------- name: Time
    ---------- units: hours since 2000-01-01 00:00:00
    ---------- length: 24
    >>> with Dataset("example.nc", "w") as nc:
    ...     axis.write(nc)

API Reference
"""""""""""""

.. automodule:: dephycf.Axis
   :members:
   :undoc-members:
   :show-inheritance:


