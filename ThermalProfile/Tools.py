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
albedo = 0.09 # (Li et al., 2006)
emissivity = 1 - albedo
StefanBoltzmannConstant = 5.670374419e-8
solarConstant = 1.37e3 # (Imke de Pater & Jack J. Lissauer, Fundamental Planetary Science)
semiMajorAxis = 2.77 # AU (Imke de Pater & Jack J. Lissauer, Fundamental Planetary Science)
surfaceTemperature = ((solarConstant / semiMajorAxis ** 2 * (1 - albedo) / (4 * emissivity * StefanBoltzmannConstant))
                      ** (1/4))

serpentiniteConductivity = 2.7  # (Osako et al., 2010)
serpentiniteHeatCapacity = 1e3  # (Osako et al., 2010)
serpentiniteDensity = 2600

iceConductivity = 3.48  # Engineering Toolbox
iceHeatCapacity = 1389  # Engineering Toolbox
iceDensity = 917

crustDensity = 1570 # From crust modelling
crustRockFraction = (crustDensity - iceDensity) / (serpentiniteDensity - iceDensity)

crustConductivity = crustRockFraction * serpentiniteConductivity + (1 - crustRockFraction) * iceConductivity
crustHeatCapacity = crustRockFraction * serpentiniteHeatCapacity + (1 - crustRockFraction) * iceHeatCapacity
k_crust = crustConductivity / crustDensity / crustHeatCapacity

mantleConductivity = serpentiniteConductivity
mantleDensity = 2520 # From crust modelling
mantleHeatCapacity = serpentiniteHeatCapacity
k_mantle = mantleConductivity / mantleDensity / mantleHeatCapacity

freezingDepressionConstant = 1.86 # Aylward, Gordon; Findlay, Tristan (2002), SI Chemical Data 5th ed.
NaClSolutionParticles = 2
NaClMolarMass = 58.443

# heatGeneration = serpentiniteDensity * (64e-9 * 1.69e-12 * 0.5 ** (4.5e9/14e9)
#                   + 16e-9 * 1.5e-12 * 0.5 ** (4.5e9/4.47e9)
#                   + 160e-6 * 5.56e-12 * 0.5 ** (4.5e9/1.25e9)
#                   + 16e-9 * 6.46e-14 * 0.5 ** (4.5e9/704e6)) # Treimann, 1986; Grott Breuer, 2008
# print(heatGeneration)

def continuityThermalProfileCoefficients(mantleHeatGeneration, crustHeatGeneration):
    # calculates coefficients for the thermal profile give the constraint of continuity of the gradients between layers
    # T(r) = - (Hr^2 / 6k) - (c1 / kr) + c2
    c1crust = boundaryRadius ** 3 / 3 * (crustHeatGeneration - k_crust / k_mantle * mantleHeatGeneration)
    c2crust = surfaceTemperature + crustHeatGeneration * R ** 2 / (6 * k_crust) + c1crust / (k_crust * R)

    boundaryTemperature = (- crustHeatGeneration * boundaryRadius ** 2 / (6 * k_crust)
                           - c1crust / (k_crust * boundaryRadius) + c2crust)

    c2mantle = boundaryTemperature + mantleHeatGeneration * boundaryRadius ** 2 / (6 * k_mantle)

    return c1crust, c2crust, c2mantle
