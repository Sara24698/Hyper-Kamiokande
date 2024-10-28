import numpy as np
from scipy.interpolate import make_interp_spline
import matplotlib.pyplot as plt

# Dataset
x = np.array([1, 1.6, 2, 2.5, 3, 3.5, 4])
y = np.array([49.63, 18.81, 16.01, 14.28, 12.63, 11.35, 9.74])

x2 = np.array([1, 1.6, 2, 2.5, 3, 3.5, 4])
y2 = np.array([28.07, 17.49, 16.53, 14.93, 13.63, 12.23, 10.61])

x3 = np.array([1, 1.6, 2, 2.5, 3, 3.5, 4])
y3 = np.array([19.54, 17.81, 16.86, 15.72, 14.4, 13.04, 11.66])

x4 = np.array([1, 1.6, 2, 2.5, 3, 3.5, 4])
y4 = np.array([65.73, 25.99, 18.69, 13.85, 11.8, 10.34, 8.69])

X_Y_Spline = make_interp_spline(x, y)
X_Y_Spline2 = make_interp_spline(x2, y2)
X_Y_Spline3 = make_interp_spline(x3, y3)
X_Y_Spline4 = make_interp_spline(x4, y4)


# Returns evenly spaced numbers
# over a specified interval.
X_ = np.linspace(x.min(), x.max(), 500)
Y_ = X_Y_Spline(X_)

X2_ = np.linspace(x2.min(), x2.max(), 500)
Y2_ = X_Y_Spline2(X2_)

X3_ = np.linspace(x3.min(), x3.max(), 500)
Y3_ = X_Y_Spline3(X3_)

X4_ = np.linspace(x4.min(), x4.max(), 500)
Y4_ = X_Y_Spline4(X4_)

# Plotting the Graph
plt.plot(X_, Y_, label='3 m')
plt.plot(X2_, Y2_, label='2 m')
plt.plot(X3_, Y3_, label='1 m')
plt.plot(X4_, Y4_, label='4 m')
plt.grid()
#plt.ylim(0,80)
#plt.title('% PMTs with excess as a function of distance to coils ')
plt.xlabel("Distance from PMTs to coils (m)", fontsize=15)
plt.ylabel("Proportion of PMTs with excess (%)", fontsize=15)
plt.legend(title='Distance between coils (m)')
plt.show()

