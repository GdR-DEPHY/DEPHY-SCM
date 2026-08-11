Installation
============


Requirements
""""""""""""

``dephycf`` is based on Python 3. Python 3.9 has been tested positively.
It needs the following packages:

- `netCDF4 <https://unidata.github.io/netcdf4-python/>`_
- `numpy <https://numpy.org/>`_
- `scipy <https://scipy.org/>`_
- `xarray <https://docs.xarray.dev/en/stable/>`_ (only for MAGIC cases)
- `matplotlib <https://matplotlib.org/>`_

Optional packages:

- `metpy <https://unidata.github.io/MetPy/latest/index.html>`_ (to compute or inverse relative humidity)

Install ``dephycf``
"""""""""""""""""""

A simple clone of the Github repository should work:

Note that the full repository has currently a size of ??, as it comes with already prepared netCDF files for each case. Future versions will remove these netCDF files and make them computed at installation.

Using ``dephycf``
"""""""""""""""""

The tools is provided as a Python module named dephycf. To use it just update your ``PYTHONPATH``, e.g., sourcing a ``setenv`` file, where you have defined where dephycf is installed: 

.. code-block:: bash
   
   cat setenv
   >> export PYTHONPATH=${DEPHYCF}:$PYTHONPATH
   source setenv. 

You may need to update it. A pip install is in progress.
