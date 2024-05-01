"""
This file contains some code that is used multiple times, to avoid duplicate code.
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import scipy.linalg
import scipy.optimize

# constants
G = 6.67430e-11
R = 470e3  # [m], park 2016
M = 62.62736e9 / G  # [kg], konopliv 2018
I = 0.375 * M * R ** 2  #[kg*m2] mao and mckinnon, 2018

iceDensity = 917  # [kg/m3]
olivineDensity = 3320  # [kg/m3]
waterDensity = 1000  # [kg/m3]
dissolvedSalt = 360  # [kg/m3]
brineDensity = waterDensity + dissolvedSalt  # [kg/m3]

waterMolarMass = 2 * 1.008 + 15.999  # [g/mol]
olivineMolarMass = (55.845 + 24.305) + 28.085 + 4 * 15.999  # [g/mol], same amount of Fe and Mg in the olivine
saltMolarMass = 22.99 + 35.45  # [g/mol], NaCl

serpentiniteDensity = 2600
serpentiniteMolarMass = 3 * (55.845 + 24.305) / 2 + 2 * 28.085 + 9 * 15.999 + 4 * 1.008

ironDensity = 7874
ironMolarMass = 55.845

def profileCalculator(rlist, profile, dr):
    # calculates mass [kg], gravity [m/s2] and pressure [Pa] for a density profile
    # rlist is the numpy array of radial values at which the profile is evaluated
    # profile is a numpy array with the densities from r=0 to r=R
    # dr is the step size of the profile

    # check if the inputs are self-consistent
    if rlist.shape[0] != round(R/dr) + 1:
        raise ValueError("rlist does not have the expected length based on given dr \n"+\
                         "Expected length (R="+str(R)+", dr="+str(dr)+"): " + str(round(R/dr)+1)+\
                         "\nActual length: " + str(rlist.shape[0]))

    # initialize list of masses
    Mlist = np.zeros(rlist.shape[0])

    # integrate mass
    for i, r in enumerate(rlist[:-1]):
        dMdr = 4 * np.pi * profile[i] * r ** 2
        Mlist[i + 1] = Mlist[i] + dMdr * dr

    # calculate gravitational acceleration
    glist = (G * Mlist) / rlist ** 2
    # correct for division by 0
    glist[0] = 0

    # initialize list of pressures
    plist = np.zeros((rlist.shape[0]))

    # integrate pressure
    for i, r in enumerate(np.flip(rlist)[:-1]):
        dpdr = G * 4 / 3 * np.pi * np.flip(profile)[i] ** 2 * r
        plist[i + 1] = plist[i] + dpdr * dr

    # flip list of pressures right way up
    plist = np.flip(plist)

    return Mlist, glist, plist