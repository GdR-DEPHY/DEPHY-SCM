**********************************************
Plotting Utilities (:mod:`dephycf.plotbasics`)
**********************************************

This page documents the plotting utilities provided in
``dephycf.plotbasics``.

These functions are wrappers around Matplotlib that simplify
the creation of line and contour plots for DEPHY-CF case analysis.

Examples
--------

1D Line Plot
~~~~~~~~~~~~

.. code-block:: python

   import numpy as np
   from dephycf import plotting

   x = np.linspace(0, 10, 100)
   y = np.sin(x)

   plotting.plot(x, y, name="line_plot.png", xlabel="X", ylabel="sin(x)")

This creates and saves a simple line plot ``line_plot.png``, which is displayed below.

.. plot::

   import numpy as np
   from dephycf import plotbasics

   x = np.linspace(0, 10, 100)
   y = np.sin(x)

   # Call plotting function
   plotbasics.plot(x, y, name="line_plot.png", xlabel="X", ylabel="sin(x)")

   # Reload the saved figure for doc display
   import matplotlib.pyplot as plt
   img = plt.imread("line_plot.png")
   plt.imshow(img)
   plt.axis("off")

2D Contour Plot
~~~~~~~~~~~~~~~

.. code-block:: python

   import numpy as np
   from dephycf import plotting

   x = np.linspace(0, 2*np.pi, 50)
   y = np.linspace(0, 1, 20)
   X, Y = np.meshgrid(x, y, indexing="ij")
   Z = np.sin(X) * np.cos(2*np.pi*Y)

   plotting.plot2D(x, y, Z, name="contour_plot.png",
                   xlabel="X", ylabel="Y")

This creates and saves a filled contour plot ``contour_plot.png``, which is displayed below.

.. plot::

   import numpy as np
   from dephycf import plotbasics

   x = np.linspace(0, 2*np.pi, 50)
   y = np.linspace(0, 1, 20)
   X, Y = np.meshgrid(x, y, indexing="ij")
   Z = np.sin(X) * np.cos(2*np.pi*Y)

   # Call plotting function
   plotbasics.plot2D(x, y, Z, name="contour_plot.png",
                   xlabel="X", ylabel="Y")

   # Reload the saved figure for doc display
   import matplotlib.pyplot as plt
   img = plt.imread("contour_plot.png")
   plt.imshow(img)
   plt.axis("off")

Module API
----------

.. automodule:: dephycf.plotbasics
   :members:
   :undoc-members:
   :show-inheritance:
