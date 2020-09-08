import numpy as np
from math import pi


rho = 3.0
X = 2.0

phi = np.linspace(0, 2*pi, 12)
x = [X for p in phi]
y = rho * np.cos(phi)
z = rho * np.sin(phi)

np.savetxt('Points.dat', np.c_[x, y, z], fmt = '%f %f %f')

