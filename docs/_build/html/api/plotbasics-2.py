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