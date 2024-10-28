import numpy as np
from scipy.interpolate import make_interp_spline
import matplotlib.pyplot as plt


# Dataset
x = np.array(['N', 'O', 'P', 'Q', 'R', 'S', 'T', 'U', 'V'])

y = np.array([0.419931, 0.522658, 0.533691, 0.560025, 0.587951, 0.604969, 0.611993, 0.599856, 0.583438])
error= np.array([0.75216, 0.864439, 0.857716, 0.871517, 0.88493, 0.8944, 0.902694])

# Plotting the Graph
plt.scatter(x, y)
#plt.errorbar(x, y, error, fmt='None', capsize = 3)
plt.grid()
plt.xlabel("Set of coils", fontsize=15)
plt.ylabel("Loss of efficiency (%)", fontsize=15)
plt.title('Failure rectangular coils')
plt.show()
