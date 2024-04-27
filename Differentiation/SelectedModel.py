import scipy.optimize
import matplotlib
import matplotlib.pyplot as plt
from Tools import *

def innerRadiusDifference(innerDensity, outerDensity):
    innerRadius1 = (3 * M / (4 * np.pi) - outerDensity * R ** 3) / (innerDensity - outerDensity)
    innerRadius1 = np.sign(innerRadius1) * (np.abs(innerRadius1)) ** (1 / 3)
    innerRadius2 = (15 * I / (8 * np.pi) - outerDensity * R ** 5) / (innerDensity - outerDensity)
    innerRadius2 = np.sign(innerRadius2) * (np.abs(innerRadius2)) ** (1 / 5)
    return innerRadius1 - innerRadius2

def determine2LayerModel(outerDensity, innerDensityGuess):
    # solves 2 layer model constrained by mass and moment of inertia
    # innerDensity in kg/m^3 (same as g/L)
    innerDensity = scipy.optimize.fsolve(innerRadiusDifference, innerDensityGuess, outerDensity)
    innerRadius = (3 * M / (4 * np.pi) - outerDensity * R ** 3) / (innerDensity[0] - outerDensity)
    innerRadius = np.sign(innerRadius) * (np.abs(innerRadius)) ** (1 / 3)
    return innerDensity[0], innerRadius


if __name__ == "__main__":
    innerDensityGuess = 3000

    # rho > 1660 from Bland et al., 2016 assuming NaCl salt from De Sanctis et al., 2016
    # rho < 2000 from TwoLayerModels
    outerDensityRange = np.arange(1650, 1970 + 10, 10)

    dr = 1
    rlist = np.arange(0, R + dr, dr)

    matplotlib.rcParams.update({'font.size': 18})
    fig, axs = plt.subplots(1, 2, figsize=(10, 6))
    fig.supxlabel("Outer Density [kg/m$^{3}$]")
    axs[0].set_ylabel("Inner Density [kg/m$^{3}$]")
    axs[1].set_ylabel("Crustal Thickness [km]")

    selectedValue = 1750
    innerDensity, innerRadius = determine2LayerModel(selectedValue, innerDensityGuess)

    axs[0].axvline(selectedValue, color="black", alpha=0.5, linestyle="--")
    axs[0].axhline(innerDensity, color="black", alpha=0.5, linestyle="--")
    secondXaxis0 = axs[0].twiny()
    secondXaxis0.set_xticks([selectedValue])
    secondXaxis0.set_xticklabels([str(selectedValue)], fontsize=12)
    secondYaxis0 = axs[0].twinx()
    secondYaxis0.set_yticks([innerDensity])
    secondYaxis0.set_yticklabels([str(round(innerDensity))], fontsize=12)

    axs[1].axvline(selectedValue, color="black", alpha=0.5, linestyle="--")
    axs[1].axhline((R - innerRadius)/1e3, color="black", alpha=0.5, linestyle="--")
    secondXaxis1 = axs[1].twiny()
    secondXaxis1.set_xticks([selectedValue])
    secondXaxis1.set_xticklabels([str(selectedValue)], fontsize=12)
    secondYaxis1 = axs[1].twinx()
    secondYaxis1.set_yticks([(R - innerRadius) / 1e3])
    secondYaxis1.set_yticklabels([str(round((R - innerRadius) / 1e3))], fontsize=12)

    baseI = I
    for IModification in [-0.001, 0, 0.001]:
        I = baseI + IModification * M * R ** 2
        Inorm = I / (M * R ** 2)

        innerDensityList = []
        crustalThicknessList = []

        for outerDensity in outerDensityRange:
            innerDensity, innerRadius = determine2LayerModel(outerDensity, innerDensityGuess)
            innerDensityList.append(innerDensity)
            crustalThicknessList.append((R - innerRadius) / 1e3)

        axs[0].plot(outerDensityRange, innerDensityList, label=r"$\frac{I}{MR^{2}}$="+str(Inorm))
        axs[1].plot(outerDensityRange, crustalThicknessList, label=r"$\frac{I}{MR^{2}}$="+str(Inorm))

    I = baseI

    secondXaxis0.set_xlim(axs[0].get_xlim())
    secondYaxis0.set_ylim(axs[0].get_ylim())
    secondXaxis1.set_xlim(axs[1].get_xlim())
    secondYaxis1.set_ylim(axs[1].get_ylim())

    axs[0].legend()
    axs[1].legend()
    plt.tight_layout()
    plt.savefig(r"Images/SelectedModelRange.pdf")
    plt.savefig(r"Images/SelectedModelRange.png")
    plt.show()