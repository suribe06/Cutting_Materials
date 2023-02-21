import numpy as np
from scipy import interpolate
import matplotlib.pyplot as plt

# Load data from .dat file
data_CrN = np.loadtxt('curva_morse_CrN.dat')
x_CrN = data_CrN[:, 0]  # distance
y_CrN = data_CrN[:, 1] * 13.6056980658938 # energy
# Create an interpolation function for interaction potential between Cr and N atoms
f = interpolate.interp1d(x_CrN, y_CrN, kind='cubic', fill_value="extrapolate")

# Load data from .dat file
data_CrCr = np.loadtxt('curva_morse_Cr.dat')
x_CrCr = data_CrCr[:, 0]  # distance
y_CrCr = data_CrCr[:, 1] * 13.6056980658938 # energy
# Create an interpolation function for interaction potential between Cr atoms
g = interpolate.interp1d(x_CrCr, y_CrCr, kind='cubic', fill_value="extrapolate")

x_min = np.min([np.min(x_CrCr), np.min(x_CrN)])
x_max = np.min([np.max(x_CrCr), np.max(x_CrN)])

x_values = np.linspace(x_min, x_max, 500)
energy_CrN = f(x_values)
energy_CrCr = g(x_values)

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 5))
ax1.plot(x_values, energy_CrN, color='m')
ax1.set_xlabel('Distance (Å)')
ax1.set_ylabel('Interaction Energy (eV)')
ax1.set_title('Interaction potential between Cr-N atoms')
ax2.plot(x_values, energy_CrCr, color='b')
ax2.set_xlabel('Distance (Å)')
ax2.set_ylabel('Interaction Energy (eV)')
ax2.set_title('Interaction potential between Cr-Cr atoms')
plt.tight_layout()
plt.show()
