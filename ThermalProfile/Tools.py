"""
This file contains some code that is used multiple times, to avoid duplicate code.
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import scipy.linalg
import scipy.optimize

# constants
R = 470e3  # [m], park 2016
surfaceTemperature = 230
boundaryRadius = R - 70e3
emissivity = 0.91
StefanBoltzmannConstant = 5.670374419e-8

crustConductivity = 2.22
crustDensity = 1570
crustHeatCapacity = 2108
k_crust = crustConductivity / crustDensity / crustHeatCapacity

mantleConductivity = 2.7
mantleDensity = 2520
mantleHeatCapacity = 1e3
k_mantle = mantleConductivity / mantleDensity / mantleHeatCapacity

heatGeneration = 2162 * (2162 - 917) / (2600 - 917) * (430e-9 * 0.0619e-3 * 0.5 ** (4.5e9/1250e6)
                                                       + 130e-9 * 0.0204e-3 * 0.5 ** (4.5e9/14000e6)
                                                       + 17.5e-9 * 0.401e-3 * 0.5 ** (4.5e9/704e6)
                                                       + 52.4e-9 * 0.104e-3 * 0.5 ** (4.5e9/4470e6))
print(heatGeneration)
heatGeneration = 5e-15

def radiationThermalProfileCoefficients():
    # calculates coefficients for the thermal profile given the constraint of outgoing thermal radiation
    # T(r) = - (Hr^2 / 6k) - (c1 / kr) + c2
    c1crust = k_crust * R ** 2 * (heatGeneration * R / (3 * k_crust) + emissivity * StefanBoltzmannConstant
                              * surfaceTemperature ** 4)
    c2crust = surfaceTemperature + heatGeneration * R ** 2 / (6 * k_crust) + c1crust / (k_crust * R)

    boundaryTemperature = (- heatGeneration * boundaryRadius ** 2 / (6 * k_crust) - c1crust / (k_crust * boundaryRadius)
                           + c2crust)

    c2mantle = boundaryTemperature + heatGeneration * boundaryRadius ** 2 / (6 * k_mantle)

    return c1crust, c2crust, c2mantle

def continuityThermalProfileCoefficients():
    # calculates coefficients for the thermal profile give the constraint of continuity of the gradients between layers
    # T(r) = - (Hr^2 / 6k) - (c1 / kr) + c2
    c1crust = heatGeneration * boundaryRadius ** 3 / 3 * (1 - k_crust / k_mantle)
    c2crust = surfaceTemperature + heatGeneration * R ** 2 / (6 * k_crust) + c1crust / (k_crust * R)

    boundaryTemperature = (- heatGeneration * boundaryRadius ** 2 / (6 * k_crust) - c1crust / (k_crust * boundaryRadius)
                           + c2crust)

    c2mantle = boundaryTemperature + heatGeneration * boundaryRadius ** 2 / (6 * k_mantle)

    return c1crust, c2crust, c2mantle
