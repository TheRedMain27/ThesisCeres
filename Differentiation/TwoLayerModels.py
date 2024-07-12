from Tools import *

def determine2LayerDensities(innerRadius):
    # solves 2 layer densities constrained by mass and moment of inertia
    # innerRadius in meters
    A = np.array([[float(innerRadius ** 3), float(R ** 3 - innerRadius ** 3)],
                   [float(innerRadius ** 5), float(R ** 5 - innerRadius ** 5)]])
    b = np.array([M / (4/3 * np.pi), I / (8/15 * np.pi)])
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
        if rhoInner < 4000 and rhoOuter > 900:
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
    fig.supylabel(r"$r_{bound}$ [km]")

    ax2 = axs[0].twiny()
    axs[0].axvline(innerDensity[-1], color="black", alpha=0.5, linestyle="--")
    axs[0].axvline(outerDensity[0], color="black", alpha=0.5, linestyle="--")
    ax2.set_xticks([outerDensity[0], innerDensity[-1]])
    ax2.set_xticklabels([str(round(outerDensity[0])), str(round(innerDensity[-1]))], fontsize=12)

    axs[0].plot(innerDensity, np.array(feasibleRadius) / 1e3, label=r"$\rho_{m}$")
    axs[0].plot(outerDensity, np.array(feasibleRadius) / 1e3, label=r"$\rho_{c}$")
    ax2.set_xlim(axs[0].get_xlim())
    axs[0].set_xlabel("Density [kg/m$^{3}$]")
    axs[0].legend()
    axs[1].plot(Pmax, np.array(feasibleRadius) / 1e3)
    axs[1].set_xlabel("Maximum Pressure [MPa]")

    plt.tight_layout()
    plt.savefig(r"Images/TwoLayerModels.pdf")
    plt.savefig(r"Images/TwoLayerModels.png")
    plt.show()