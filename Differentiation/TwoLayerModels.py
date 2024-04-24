import math
import numpy as np
import scipy.linalg
import matplotlib
import matplotlib.pyplot as plt
from Tools import *

def determine2LayerDensities(innerRadius):
    # solves 2 layer densities constrained by mass and moment of inertia
    # innerRadius in meters
    A = np.array([[float(innerRadius ** 3), float(R ** 3 - innerRadius ** 3)],
                   [float(innerRadius ** 5), float(R ** 5 - innerRadius ** 5)]])
    b = np.array([M / (4/3 * math.pi), I / (8/15 * math.pi)])
    rho = scipy.linalg.solve(A, b)
    return rho[0], rho[1]

if __name__ == "__main__":
    dr = 1
    rlist = np.arange(0, R + dr, dr)

    mantleRadius = np.arange(200e3, 469e3, 10e3)
    innerDensity = []
    outerDensity = []
    feasibleRadius = []
    Pmax = []

    for radius in mantleRadius:
        rhoInner, rhoOuter = determine2LayerDensities(radius)
        if rhoInner < 4000 and rhoOuter > 1000:
            profile = np.ones(rlist.shape[0])
            profile[:round(radius/dr)+1] *= rhoInner
            profile[round(radius/dr)+1:] *= rhoOuter
            MVector, gVector, PVector = profileCalculator(rlist, profile, dr)

            print("Mantle radius " + str(radius) + " [m]")
            print("innerDensity = " + str(rhoInner) + " [kg/m3]")
            print("outerDensity = " + str(rhoOuter) + " [kg/m3]")
            print("P = " + str(PVector[0]/1e6) + " [MPa]")

            feasibleRadius.append(radius)
            innerDensity.append(rhoInner)
            outerDensity.append(rhoOuter)
            Pmax.append(PVector[0]/1e6)

    matplotlib.rcParams.update({'font.size': 18})
    fig, axs = plt.subplots(1, 2, figsize=(10, 6))
    fig.supylabel("Layer Boundary Radius [km]")
    axs[0].plot(innerDensity, np.array(feasibleRadius) / 1e3, label = r"$\rho_{inner}$")
    axs[0].plot(outerDensity, np.array(feasibleRadius) / 1e3, label = r"$\rho_{outer}$")
    axs[0].set_xlabel("Density [kg/m$^{3}$]")
    axs[0].legend()
    axs[1].plot(Pmax, np.array(feasibleRadius) / 1e3)
    axs[1].set_xlabel("Maximum Pressure [m/s$^{2}$]")
    plt.tight_layout()
    plt.savefig(r"Images/TwoLayerModels.pdf")
    plt.show()

