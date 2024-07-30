from Tools import *

def thermalProfile(rlist, constraint):
    if constraint == "radiation":
        c1crust, c2crust, c2mantle = radiationThermalProfileCoefficients()
    elif constraint == "continuity":
        c1crust, c2crust, c2mantle = continuityThermalProfileCoefficients()
    else:
        return

    crustList = rlist[rlist >= boundaryRadius]
    mantleList = rlist[rlist < boundaryRadius]

    crustTemperatures = (- heatGeneration * crustList ** 2 / (6 * k_crust) - c1crust / (k_crust * crustList)
                         + c2crust)
    mantleTemperatures = (- heatGeneration * mantleList ** 2 / (6 * k_mantle) + c2mantle)

    temperatures = np.concatenate((mantleTemperatures, crustTemperatures))

    return temperatures

if __name__ == "__main__":
    rlist = np.arange(0, 471e3, 5e3)
    temperatures = thermalProfile(rlist, "radiation")
    print(temperatures)
    temperatures = thermalProfile(rlist, "continuity")
    print(temperatures)
