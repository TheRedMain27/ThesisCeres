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

    # 2200<rho<2650 from Ruesch et al., 2019; rho>2400 from TwoLayerModels
    innerDensityRange = np.arange(2400, 2650 + 10, 10)

    dr = 1
    rlist = np.arange(0, R + dr, dr)

    matplotlib.rcParams.update({'font.size': 18})
    fig, axs = plt.subplots(1, 2, figsize=(10, 6))
    fig.supxlabel("Inner Density [kg/m$^{3}$]")
    axs[0].set_ylabel("Outer Density [kg/m$^{3}$]")
    axs[1].set_ylabel("Crustal Thickness [km]")

    selectedValue = 2500
    outerDensity, innerRadius = determine2LayerModel(selectedValue, outerDensityGuess)

    axs[0].axvline(selectedValue, color="black", alpha=0.5, linestyle="--")
    axs[0].axhline(outerDensity, color="black", alpha=0.5, linestyle="--")

    axs[1].axvline(selectedValue, color="black", alpha=0.5, linestyle="--")
    axs[1].axhline((R - innerRadius)/1e3, color="black", alpha=0.5, linestyle="--")

    baseI = I
    for IModification in [-0.001, 0, 0.001]:
        I = baseI + IModification * M * R ** 2
        Inorm = I / (M * R ** 2)

        outerDensityList = []
        crustalThicknessList = []

        for innerDensity in innerDensityRange:
            outerDensity, innerRadius = determine2LayerModel(innerDensity, outerDensityGuess)
            outerDensityList.append(outerDensity)
            crustalThicknessList.append((R - innerRadius) / 1e3)

        axs[0].plot(innerDensityRange, outerDensityList, label=r"$\frac{I}{MR^{2}}$="+str(Inorm))
        axs[1].plot(innerDensityRange, crustalThicknessList, label=r"$\frac{I}{MR^{2}}$="+str(Inorm))

    I = baseI

    axs[0].legend()
    axs[1].legend()
    plt.tight_layout()
    plt.savefig(r"Images/SelectedModelRange.pdf")
    plt.show()