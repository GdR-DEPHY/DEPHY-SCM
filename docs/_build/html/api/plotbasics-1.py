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