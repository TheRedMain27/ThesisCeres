import math
import numpy as np
import scipy.linalg
import matplotlib
import matplotlib.pyplot as plt
from Tools import *

def determineLinearCoefficients():
    # solves linear density profile constrained by mass and moment of inertia
    # returns coefficients [c0, c1], where rho = c0 - c1 * r
    A = np.array([[1/3, -R/4],
                  [1/5, -R/6]])
    b = np.array([M / (4 * math.pi * R ** 3),
                  I / (4 * math.pi * R ** 5)])
    return scipy.linalg.solve(A, b)
def determineQuadraticCoefficients(c0):
    # solves quadratic density profile constrained by mass and moment of inertia
    # returns coefficients [c1, c2], given c0, where rho = c0 + c1 * r + c2 * r^2
    A = np.array([[R/4, (R**2)/5],
                  [R/6, (R**2)/7]])
    b = np.array([M / (4 * math.pi * R ** 3) - c0 / 3,
                  I / (4 * math.pi * R ** 5) - c0 / 5])
    return scipy.linalg.solve(A, b)

def determineCubicCoefficients(c0):
    # solves cubic density profile constrained by mass and moment of inertia
    # returns coefficients [c2, c3], given c0, where rho = c0 + c2 * r^2 + c3 * r^3
    # assumption: slope of density at r = 0 is zero --> c1 = 0
    A = np.array([[(R**2)/5, (R**3)/6],
                  [(R**2)/7, (R**3)/8]])
    b = np.array([M / (4 * math.pi * R ** 3) - c0 / 3,
                  I / (4 * math.pi * R ** 5) - c0 / 5])
    return scipy.linalg.solve(A, b)

if __name__ == "__main__":
    dr = 1
    rlist = np.arange(0, R + dr, dr)

    matplotlib.rcParams.update({'font.size': 16})
    fig, axs = plt.subplots(1, 3, figsize=(15, 6))
    fig.supylabel("Radius [km]")
    fig.supxlabel("Density [kg/m$^{3}$]")

    linearCoefficients = determineLinearCoefficients()
    c0 = linearCoefficients[0]
    c1 = linearCoefficients[1]
    profile = c0 - c1 * rlist
    axs[0].plot(profile, rlist / 1e3)
    axs[0].set_title("Linear")

    c0list = np.arange(1000, 5000, 1000)
    for c0 in c0list:
        coefficients = determineQuadraticCoefficients(c0)
        c1 = coefficients[0]
        c2 = coefficients[1]

        profile = c0 + c1 * rlist + c2 * rlist ** 2

        axs[1].plot(profile, rlist/1e3, label = r"$\rho_{core}$ = " + str(c0))
        axs[1].legend()
        axs[1].set_title("Quadratic")

    c0list = np.arange(1000, 5000, 1000)
    for c0 in c0list:
        coefficients = determineCubicCoefficients(c0)
        c2 = coefficients[0]
        c3 = coefficients[1]

        profile = c0 + c2 * rlist ** 2 + c3 * rlist ** 3

        axs[2].plot(profile, rlist / 1e3, label=r"$\rho_{core}$ = " + str(c0))
        axs[2].legend()
        axs[2].set_title("Cubic")

    plt.tight_layout()
    plt.savefig(r"Images/PolynomialModels.pdf")
    plt.show()
