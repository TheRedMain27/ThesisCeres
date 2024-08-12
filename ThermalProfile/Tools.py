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
boundaryRadius = R - 70e3
emissivity = 0.91
StefanBoltzmannConstant = 5.670374419e-8
surfaceTemperature = (1.37e3 / 2.77 ** 2 * (1 - 0.09) / (4 * emissivity * StefanBoltzmannConstant)) ** (1/4)

serpentiniteConductivity = 2.7  # (Osako et al., 2010)
serpentiniteHeatCapacity = 1e3  # (Osako et al., 2010)

iceConductivity = 3.48  # Engineering Toolbox
iceHeatCapacity = 1389  # Engineering Toolbox

crustDensity = 1570
crustRockFraction = (crustDensity - 917) / (2600 - 917)

crustConductivity = crustRockFraction * serpentiniteConductivity + (1 - crustRockFraction) * iceConductivity
crustHeatCapacity = crustRockFraction * serpentiniteHeatCapacity + (1 - crustRockFraction) * iceHeatCapacity
k_crust = crustConductivity / crustDensity / crustHeatCapacity

mantleConductivity = serpentiniteConductivity
mantleDensity = 2520
mantleHeatCapacity = serpentiniteHeatCapacity
k_mantle = mantleConductivity / mantleDensity / mantleHeatCapacity

# heatGeneration = (26.4e-6 * 2.96e-8 * 0.5 ** (4.5e9/14e9)
#                   + 94.6e-6 * 8.08e-9 * 0.5 ** (4.5e9/4.47e9)
#                   + 29.2e-6 * 5.53e-4 * 0.5 ** (4.5e9/1.25e9)
#                   + 569e-6 * 8.08e-9 * 0.5 ** (4.5e9/704e6))
# print(heatGeneration)

mantleHeatGeneration = 1.2e-14
crustHeatGeneration = crustRockFraction * mantleHeatGeneration

# def radiationThermalProfileCoefficients():
#     # calculates coefficients for the thermal profile given the constraint of outgoing thermal radiation
#     # T(r) = - (Hr^2 / 6k) - (c1 / kr) + c2
#     c1crust = k_crust * R ** 2 * (heatGeneration * R / (3 * k_crust) + emissivity * StefanBoltzmannConstant
#                               * surfaceTemperature ** 4)
#     c2crust = surfaceTemperature + heatGeneration * R ** 2 / (6 * k_crust) + c1crust / (k_crust * R)
#
#     boundaryTemperature = (- heatGeneration * boundaryRadius ** 2 / (6 * k_crust) - c1crust / (k_crust * boundaryRadius)
#                            + c2crust)
#
#     c2mantle = boundaryTemperature + heatGeneration * boundaryRadius ** 2 / (6 * k_mantle)
#
#     return c1crust, c2crust, c2mantle

def continuityThermalProfileCoefficients():
    # calculates coefficients for the thermal profile give the constraint of continuity of the gradients between layers
    # T(r) = - (Hr^2 / 6k) - (c1 / kr) + c2
    c1crust = boundaryRadius ** 3 / 3 * (crustHeatGeneration - k_crust / k_mantle * mantleHeatGeneration)
    c2crust = surfaceTemperature + crustHeatGeneration * R ** 2 / (6 * k_crust) + c1crust / (k_crust * R)

    boundaryTemperature = (- crustHeatGeneration * boundaryRadius ** 2 / (6 * k_crust)
                           - c1crust / (k_crust * boundaryRadius) + c2crust)

    c2mantle = boundaryTemperature + mantleHeatGeneration * boundaryRadius ** 2 / (6 * k_mantle)

    return c1crust, c2crust, c2mantle
