import math
import numpy as np
import scipy.linalg
import matplotlib
import matplotlib.pyplot as plt

G = 6.67430e-11
R = 470e3  # [m], park 2016
M = 62.62736e9 / G  # [kg], konopliv 2018
I = 0.37 * M * R ** 2  #[kg*m2] park, 2016

def determine2LayerDensities(innerRadius):
    # solves 2 layer densities constrained by mass and moment of inertia
    # innerRadius in meters
    A = np.array([[float(innerRadius ** 3), float(R ** 3 - innerRadius ** 3)],
                   [float(innerRadius ** 5), float(R ** 5 - innerRadius ** 5)]])
    b = np.array([M / (4/3 * math.pi), I / (8/15 * math.pi)])
    rho = scipy.linalg.solve(A, b)
    return rho[0], rho[1]

def calculate2LayerProfiles(innerRadius, innerDensity, outerDensity, dr = 1):
    # set up values of r at which variables will be calculated
    rlist = np.arange(0, R, dr)

    # set up rlists for mass calculation
    # mantle gets an additional point to calculate the first point of the crust
    rlistMantle = np.arange(0, innerRadius + dr, dr)
    rlistCrust = np.arange(innerRadius, R, dr)

    # initialize empty lists to hold calculated values
    MlistMantle = np.zeros((rlistMantle.shape[0]))
    MlistCrust = np.zeros((rlistCrust.shape[0]))

    # use Euler integrator to calculate mantle masses
    for i, r in enumerate(rlistMantle[:-1]):
        dMdrMantle = 4 * np.pi * innerDensity * r ** 2
        MlistMantle[i + 1] = MlistMantle[i] + dMdrMantle * dr

    # set first point of crust using mantle
    MlistCrust[0] = MlistMantle[-1]
    # delete last point of mantle (it will be included in crust)
    MlistMantle = MlistMantle[:-1]

    # use Euler integrator to calculate crust masses
    for i, r in enumerate(rlistCrust[:-1]):
        dMdrCrust = 4 * np.pi * outerDensity * r ** 2
        MlistCrust[i + 1] = MlistCrust[i] + dMdrCrust * dr

    # tie together mass lists and scale them
    Mlist = np.concatenate((MlistMantle, MlistCrust))
    Mlist1020 = Mlist / 1e20

    # calculate gravitational acceleration
    g = (G * Mlist) / rlist ** 2
    # correct for division by zero
    g[0] = 0

    # set up rlists for pressure calculation
    # crust gets an additional point to calculate the first point of the mantle
    rlistMantle = np.arange(0, innerRadius, dr)
    rlistCrust = np.arange(innerRadius - dr, R, dr)

    # initialize empty lists to hold calculated values
    plistMantle = np.zeros((rlistMantle.shape[0]))
    plistCrust = np.zeros((rlistCrust.shape[0]))

    # use Euler integrator to calculate crust pressures
    for i, r in enumerate(np.flip(rlistCrust)[:-1]):
        dpdrCrust = G * 4 / 3 * np.pi * outerDensity ** 2 * r
        plistCrust[i + 1] = plistCrust[i] + dpdrCrust * dr

    # set first point of mantle using crust
    plistMantle[0] = plistCrust[-1]
    # delete last point of crust (it will be included in mantle)
    plistCrust = plistCrust[:-1]

    # use Euler integrator to calculate mantle pressures
    for i, r in enumerate(np.flip(rlistMantle)[:-1]):
        dpdrMantle = G * 4 / 3 * np.pi * innerDensity ** 2 * r
        plistMantle[i + 1] = plistMantle[i] + dpdrMantle * dr

    # tie together pressure lists and scale them
    plist = np.concatenate((plistCrust, plistMantle))
    plistMPa = np.flip(plist) / 1e6

    # get radius in km
    rlistkm = rlist / 1e3

    return rlistkm, Mlist1020, g, plistMPa


if __name__ == "__main__":
    mantleRadius = np.arange(200e3, 469e3, 10e3)
    innerDensity = []
    outerDensity = []
    feasibleRadius = []
    Pmax = []
    for radius in mantleRadius:
        rhoInner, rhoOuter = determine2LayerDensities(radius)
        if rhoInner < 4000 and rhoOuter > 1000:
            rVector, MVector, gVector, PVector = calculate2LayerProfiles(radius, rhoInner, rhoOuter)
            print("Mantle radius " + str(radius) + " [m]")
            print("innerDensity = " + str(rhoInner) + " [kg/m3]")
            print("outerDensity = " + str(rhoOuter) + " [kg/m3]")
            print("P = " + str(PVector[0]) + " [MPa]")
            feasibleRadius.append(radius)
            innerDensity.append(rhoInner)
            outerDensity.append(rhoOuter)
            Pmax.append(PVector[0])

    matplotlib.rcParams.update({'font.size': 18})
    fig, axs = plt.subplots(1, 2, figsize=(10, 6))
    fig.supylabel("Radius [km]")
    axs[0].plot(innerDensity, np.array(feasibleRadius) / 1e3, label = "inner")
    axs[0].plot(outerDensity, np.array(feasibleRadius) / 1e3, label = "outer")
    axs[0].set_xlabel("Density [kg/m$^{3}$]")
    axs[0].legend()
    axs[1].plot(Pmax, np.array(feasibleRadius) / 1e3)
    axs[1].set_xlabel("Maximum Pressure [m/s$^{2}$]")
    plt.tight_layout()
    plt.savefig(r"Images/TwoLayerModels.pdf")
    plt.show()

