import sys
sys.path.append("../Differentiation")
from Tools import *
from Differentiation.Tools import *

def thermalProfile(rlist):
    c1crust, c2crust, c2mantle = continuityThermalProfileCoefficients()

    crustList = rlist[rlist >= boundaryRadius]
    mantleList = rlist[rlist < boundaryRadius]

    crustTemperatures = (- crustHeatGeneration * crustList ** 2 / (6 * k_crust) - c1crust / (k_crust * crustList)
                         + c2crust)
    mantleTemperatures = (- mantleHeatGeneration * mantleList ** 2 / (6 * k_mantle) + c2mantle)

    temperatures = np.concatenate((mantleTemperatures, crustTemperatures))

    return temperatures

if __name__ == "__main__":
    dr = 5e3
    rlist = np.arange(0, 471e3, dr)
    temperatures = thermalProfile(rlist)
    print(temperatures)

    densityProfile = np.array([mantleDensity] * len(rlist[rlist < boundaryRadius]) +
                              [crustDensity] * len(rlist[rlist >= boundaryRadius]))

    pressures = profileCalculator(rlist, densityProfile, dr)[2]

    phaseDiagramPressures = np.array([1, 2, 5, 10, 15] + list(np.arange(20, 110, 10)) + list(np.arange(120, 220, 20)))
    phaseDiagramTemperatures = np.array([-0.064, -0.14, -0.37, -0.75, -1.14, -1.54, -2.36, -3.21, -4.09, -5, -5.94,
                                         -6.91, -7.91, -8.94, -11.09, -13.35, -15.73, -18.22, -20.83]) + 273.15
    phaseDiagramTemperatures -= 1.86 * 2 * 360 / 58.443 / 917

    matplotlib.rcParams.update({'font.size': 18})
    plt.figure(figsize=(6, 6))
    plt.xlabel("Temperature [K]")
    plt.ylabel("Radius [km]")
    plt.plot(temperatures, rlist / 1e3)
    plt.tight_layout()
    plt.savefig(r"Images/TemperatureProfile.pdf")
    plt.savefig(r"Images/TemperatureProfile.png")
    plt.show()

    matplotlib.rcParams.update({'font.size': 18})
    plt.figure(figsize=(6, 6))
    plt.xlabel("Pressure [MPa]")
    plt.ylabel("Radius [km]")
    plt.plot(pressures / 1e6, rlist / 1e3)
    plt.tight_layout()
    plt.show()

    matplotlib.rcParams.update({'font.size': 18})
    plt.figure(figsize=(6, 6))
    plt.xlabel("Temperature [K]")
    plt.ylabel("Pressure [MPa]")
    plt.plot(temperatures, pressures / 1e6)
    plt.plot(phaseDiagramTemperatures, phaseDiagramPressures)
    plt.tight_layout()
    plt.savefig(r"Images/PTdiagram.pdf")
    plt.savefig(r"Images/PTdiagram.png")
    plt.show()