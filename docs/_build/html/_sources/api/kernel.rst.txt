****************************************
Kernel smoothing (:mod:`dephycf.kernel`)
****************************************

Overview
""""""""

This module implements Gaussian kernel smoothing.

The Gaussian kernel is defined as

.. math::

   w(x)=
   \exp\left(
   -\frac{(x-\mu)^2}
          {2\sigma^2}
   \right)

where :math:`\mu` is the kernel center and
:math:`\sigma` the kernel standard deviation.

Example
"""""""

.. code-block:: python

    import numpy as np
    from dephycf.kernel import gaussian_kernel_mean

    x = np.linspace(0, 100, 200)
    y = np.sin(x)

    ys = gaussian_kernel_mean(y, x, sigma=5)

API Reference
"""""""""""""

.. automodule:: dephycf.kernel
   :members:
   :undoc-members:
   :show-inheritance:
