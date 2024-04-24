import math
import numpy as np
import scipy.linalg
import matplotlib
import matplotlib.pyplot as plt
from Tools import *

def determineCoefficients():
    # solves linear density profile constrained by mass and moment of inertia
    # returns coefficients [c0, c1], where rho = c0 - c1 * r
    A = np.array([[1/3, -R/4], [1/5, -R/6]])
    b = np.array([M / (4 * math.pi * R ** 3), I / (4 * math.pi * R ** 5)])
    return scipy.linalg.solve(A, b)

if __name__ == "__main__":
    coefficients = determineCoefficients()
    c0 = coefficients[0]
    c1 = coefficients[1]

    dr = 1
    rlist = np.arange(0, R + dr, dr)
    profile = c0 - c1 * rlist

    MVector, gVector, pVector = profileCalculator(rlist, profile, dr)

    matplotlib.rcParams.update({'font.size': 18})
    fig, axs = plt.subplots(1, 3, figsize=(15, 6))
    fig.supylabel("Radius [km]")
    axs[0].plot(profile, rlist/1e3)
    axs[0].set_xlabel("Density [kg/m$^{3}$]")
    axs[1].plot(MVector/1e20, rlist/1e3)
    axs[1].set_xlabel("Integrated Mass [$10^{20}$ kg]")
    axs[2].plot(pVector/1e6, rlist/1e3)
    axs[2].set_xlabel("Pressure [MPa]")
    plt.tight_layout()
    plt.savefig(r"Images/LinearModel.pdf")
    plt.show()