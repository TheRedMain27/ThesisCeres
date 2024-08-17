import sys
sys.path.append("../Differentiation")
from Tools import *
from Differentiation.Tools import *

def thermalProfile(rlist, mantleHeatGeneration):
    crustHeatGeneration = crustRockFraction * mantleHeatGeneration

    c1crust, c2crust, c2mantle = continuityThermalProfileCoefficients(mantleHeatGeneration, crustHeatGeneration)

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
    temperatures1 = thermalProfile(rlist, 0.8e-14)
    temperatures2 = thermalProfile(rlist, 1.05e-14)
    temperatures3 = thermalProfile(rlist, 1.3e-14)

    densityProfile = np.array([mantleDensity] * len(rlist[rlist < boundaryRadius]) +
                              [crustDensity] * len(rlist[rlist >= boundaryRadius]))

    pressures = profileCalculator(rlist, densityProfile, dr)[2]

    phaseDiagramPressures = np.array([1, 2, 5, 10, 15] + list(np.arange(20, 110, 10)) + list(np.arange(120, 180, 20)))
    phaseDiagramTemperatures = np.array([-0.064, -0.14, -0.37, -0.75, -1.14, -1.54, -2.36, -3.21, -4.09, -5, -5.94,
                                         -6.91, -7.91, -8.94, -11.09, -13.35, -15.73]) + 273.15
                                        # Engineering Toolbox
    phaseDiagramTemperatures1 = (phaseDiagramTemperatures - freezingDepressionConstant * NaClSolutionParticles * 180e3
                                 / NaClMolarMass / iceDensity)
    phaseDiagramTemperatures2 = (phaseDiagramTemperatures - freezingDepressionConstant * NaClSolutionParticles * 360e3
                                 / NaClMolarMass / iceDensity)

    matplotlib.rcParams.update({'font.size': 16})
    fig, axs = plt.subplots(1, 2, figsize=(12, 6))
    fig.supxlabel("Temperature [K]")
    axs[0].set_ylabel("Radius [km]")
    axs[0].plot(temperatures1, rlist / 1e3)
    axs[0].plot(temperatures2, rlist / 1e3)
    axs[0].plot(temperatures3, rlist / 1e3)

    axs[1].set_ylabel("Pressure [MPa]")
    axs[1].plot(temperatures1, pressures / 1e6, label = r"H=0.8$\cdot 10^{-14}$")
    axs[1].plot(temperatures2, pressures / 1e6, label = r"H=1.05$\cdot 10^{-14}$")
    axs[1].plot(temperatures3, pressures / 1e6, label = r"H=1.3$\cdot 10^{-14}$")
    axs[1].plot(phaseDiagramTemperatures, phaseDiagramPressures, label = "0 g/L NaCl")
    axs[1].plot(phaseDiagramTemperatures1, phaseDiagramPressures, label = "180 g/L NaCl")
    axs[1].plot(phaseDiagramTemperatures2, phaseDiagramPressures, label = "360 g/L NaCl")
    axs[1].axhline(22.5, color = "grey", linestyle = '--')
    plt.gca().invert_yaxis()

    handles, labels = plt.gca().get_legend_handles_labels()
    fig.legend(handles, labels, loc='upper right')
    plt.tight_layout()
    plt.savefig(r"Images/Temperature.pdf")
    plt.savefig(r"Images/Temperature.png")
    plt.show()