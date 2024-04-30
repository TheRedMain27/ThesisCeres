from Tools import *
import TwoLayerModels

def testProfileCalculator():
    density = 3 * M / (4 * math.pi * R ** 3)
    dr = 1
    rlist = np.arange(0, R + dr, dr)
    profile = np.ones(rlist.shape[0]) * density
    Mlist, glist, plist = profileCalculator(rlist, profile, dr)

    rlistkm = rlist / 1e3
    Mlist1020 = Mlist / 1e20
    plistMPa = plist / 1e6

    matplotlib.rcParams.update({'font.size': 18})
    fig, axs = plt.subplots(1, 3, figsize=(15, 6))
    fig.supylabel("Radius [km]")
    axs[0].plot(Mlist1020, rlistkm)
    axs[0].set_xlabel("Integrated Mass [$10^{20}$ kg]")
    axs[1].plot(glist, rlistkm)
    axs[1].set_xlabel("Gravity [m/s$^{2}$]")
    axs[2].plot(plistMPa, rlistkm)
    axs[2].set_xlabel("Pressure [MPa]")
    plt.tight_layout()
    plt.show()

def test2LayerDensity():
    testRadius = [1, 100000, 300000, 450000, 469999]
    for i, innerRadius in enumerate(testRadius):
        rhoInner, rhoOuter = TwoLayerModels.determine2LayerDensities(innerRadius)
        Mdiff = (TwoLayerModels.M - 4 / 3 * math.pi *
                 (rhoOuter * TwoLayerModels.R ** 3 + (rhoInner - rhoOuter) * innerRadius ** 3)) / TwoLayerModels.M
        Idiff = (TwoLayerModels.I - 8 / 15 * math.pi *
                 (rhoOuter * TwoLayerModels.R ** 5 + (rhoInner - rhoOuter) * innerRadius ** 5)) / TwoLayerModels.I
        print("test " + str(i + 1) + " : " + str((Mdiff < 1e-10 and Idiff < 1e-10)))

if __name__ == "__main__":
    testProfileCalculator()