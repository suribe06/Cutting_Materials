import numpy as np
from scipy import interpolate
import matplotlib.pyplot as plt

eV_factor = 13.6056980658938
max_xvalue = 11

def energyFunctionCrN():
    # Load data from .dat file
    data_CrN = np.loadtxt('curva_morse_CrN.dat')
    x_CrN = data_CrN[:, 0]  # distance
    y_CrN = data_CrN[:, 1] * eV_factor # energy
    # Create an interpolation function for interaction potential between Cr and N atoms
    f = interpolate.interp1d(x_CrN, y_CrN, kind='cubic', fill_value="extrapolate")
    def new_f(x):
        return f(x) - f(max_xvalue)
    return new_f

def energyFunctionCrCr():
    # Load data from .dat file
    data_CrCr = np.loadtxt('curva_morse_Cr.dat')
    x_CrCr = data_CrCr[:, 0]  # distance
    y_CrCr = data_CrCr[:, 1] * eV_factor # energy
    # Create an interpolation function for interaction potential between Cr atoms
    g = interpolate.interp1d(x_CrCr, y_CrCr, kind='cubic', fill_value="extrapolate")
    def new_g(x):
        return g(x) - g(max_xvalue)
    return new_g

def energyFunctionNN():
    # Load data from .dat file
    data_NN = np.loadtxt('bond-lenght.dat')
    x_NN = data_NN[:, 0]  # distance
    y_NN = data_NN[:, 1] * eV_factor # energy
    # Create an interpolation function for interaction potential between N atoms
    h = interpolate.interp1d(x_NN, y_NN, kind='cubic', fill_value="extrapolate")
    def new_h(x):
        return h(x) - h(max_xvalue)
    return new_h

def testMethod():
    x_values = np.linspace(1, max_xvalue, 500)
    f = energyFunctionCrN()
    g = energyFunctionCrCr()
    h = energyFunctionNN()
    energy_CrN = f(x_values)
    energy_CrCr = g(x_values)
    energy_NN = h(x_values)

    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(12, 6))
    ax1.plot(x_values, energy_CrN, color='m')
    ax1.set_xlabel('Distance (Å)')
    ax1.set_ylabel('Interaction Energy (eV)')
    ax1.set_title('Cr-N atoms')
    ax2.plot(x_values, energy_CrCr, color='b')
    ax2.set_xlabel('Distance (Å)')
    ax2.set_ylabel('Interaction Energy (eV)')
    ax2.set_title('Cr-Cr atoms')
    ax3.plot(x_values, energy_NN, color='r')
    ax3.set_xlabel('Distance (Å)')
    ax3.set_ylabel('Interaction Energy (eV)')
    ax3.set_title('N-N atoms')
    plt.tight_layout()
    plt.show()

#testMethod()