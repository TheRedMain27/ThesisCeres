import scipy.optimize
import matplotlib
import matplotlib.pyplot as plt
from Tools import *

def innerRadiusDifference(outerDensity, innerDensity):
    innerRadius1 = (3 * M / (4 * np.pi) - outerDensity * R ** 3) / (innerDensity - outerDensity)
    innerRadius1 = np.sign(innerRadius1) * (np.abs(innerRadius1)) ** (1 / 3)
    innerRadius2 = (15 * I / (8 * np.pi) - outerDensity * R ** 5) / (innerDensity - outerDensity)
    innerRadius2 = np.sign(innerRadius2) * (np.abs(innerRadius2)) ** (1 / 5)
    return innerRadius1 - innerRadius2

def determine2LayerModel(innerDensity, outerDensityGuess):
    # solves 2 layer model constrained by mass and moment of inertia
    # innerDensity in kg/m^3 (same as g/L)
    outerDensity = scipy.optimize.fsolve(innerRadiusDifference, outerDensityGuess, innerDensity)
    innerRadius = (3 * M / (4 * np.pi) - outerDensity[0] * R ** 3) / (innerDensity - outerDensity[0])
    innerRadius = np.sign(innerRadius) * (np.abs(innerRadius)) ** (1 / 3)
    return outerDensity[0], innerRadius


if __name__ == "__main__":
    outerDensityGuess = 2000
    innerDensityRange = [2400, 2450, 2525, 2650]  # from Ruesch et al., 2019
    labels = [r"$\rho_{mantle}$ = " + str(density) for density in innerDensityRange]

    dr = 1
    rlist = np.arange(0, R + dr, dr)

    matplotlib.rcParams.update({'font.size': 18})
    fig, axs = plt.subplots(1, 2, figsize=(10, 6))
    fig.supylabel("Radius [km]")
    axs[0].set_xlabel("Density [kg/m$^{3}$]")
    axs[1].set_xlabel("Pressure [MPa]")

    for i, innerDensity in enumerate(innerDensityRange):
        outerDensity, innerRadius = determine2LayerModel(innerDensity, outerDensityGuess)

        print("Inner Density = " + str(innerDensity) + " [kg/m3]")
        print("Outer Density = " + str(outerDensity) + " [kg/m3]")
        print("Crustal Thickness = " + str((R - innerRadius) / 1e3) + " [km]")

        profile = np.ones(rlist.shape[0])
        profile[:round(innerRadius / dr) + 1] *= innerDensity
        profile[round(innerRadius / dr) + 1:] *= outerDensity
        MVector, gVector, PVector = profileCalculator(rlist, profile, dr)

        axs[0].plot(profile, rlist / 1e3, label = labels[i])
        axs[1].plot(PVector / 1e6, rlist / 1e3, label = labels[i])

    axs[0].legend()
    axs[1].legend()
    plt.tight_layout()
    plt.savefig(r"Images/SelectedModels.pdf")
    plt.show()