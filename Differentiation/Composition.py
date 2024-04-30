from Tools import *
from TwoLayerModels import determine2LayerDensities

def determineCompositions(outerDensity, innerDensity):
    # Determines compositions of inner and outer layers from their densities
    # Assumptions:
    # 1. crust is made up of water ice and serpentinite
    # 2. mantle is made up of brine with 360 g/L dissolved NaCl and serpentinite
    # 3. serpentinite density is 3320 kg/m3

    crustserpentiniteFraction = (outerDensity - iceDensity) / (serpentiniteDensity - iceDensity)
    mantleserpentiniteFraction = (innerDensity - brineDensity) / (serpentiniteDensity - brineDensity)

    return crustserpentiniteFraction, mantleserpentiniteFraction

def atomCounter(mantleRadius, mantleDensity, crustDensity):
    # determines atomic proportions based on densities and crustal thickness
    # returns array with atomic proportions for [O, Si, Fe+Mg(50/50), H, Na+Cl (50/50)]

    mantleVolume = (4/3) * np.pi * mantleRadius ** 3
    mantleMass = mantleVolume * mantleDensity
    crustVolume = (4/3) * np.pi * R ** 3 - mantleVolume
    crustMass = crustVolume * crustDensity

    crustFraction, mantleFraction = determineCompositions(crustDensity, mantleDensity)

    serpentiniteMass = crustFraction * crustMass + mantleFraction * mantleMass
    iceMass = (1 - crustFraction) * crustMass
    brineMass = (1 - mantleFraction) * mantleMass
    waterMass = (waterDensity / brineDensity) * brineMass
    saltMass = (dissolvedSalt / brineDensity) * brineMass

    serpentiniteMolecules = serpentiniteMass / serpentiniteMolarMass  # [kmol]
    waterMolecules = (iceMass + waterMass) / waterMolarMass
    saltMolecules = saltMass / saltMolarMass

    # dependent on which silicate is used!
    oxygenAtoms = 9 * serpentiniteMolecules + waterMolecules
    siliconAtoms = 2 * serpentiniteMolecules
    metalAtoms = 3 * serpentiniteMolecules
    hydrogenAtoms = 4 * serpentiniteMolecules + 2 * waterMolecules
    saltAtoms = 2 * saltMolecules

    totalAtoms = oxygenAtoms + siliconAtoms + metalAtoms + hydrogenAtoms + saltAtoms

    return np.array([oxygenAtoms, siliconAtoms, metalAtoms, hydrogenAtoms, saltAtoms]) / totalAtoms

if __name__ == "__main__":
    dr = 1
    rlist = np.arange(0, R + dr, dr)

    innerRadius = np.arange(240e3, 460e3, 10e3)

    feasibleRadius = []
    oxygenFraction = []
    siliconFraction = []
    metalFraction = []
    hydrogenFraction = []
    saltFraction = []

    for radius in innerRadius:
        rhoInner, rhoOuter = determine2LayerDensities(radius)
        if rhoInner <= serpentiniteDensity and rhoOuter >= iceDensity:
            feasibleRadius.append(radius)

            proportions = atomCounter(radius, rhoInner, rhoOuter)

            oxygenFraction.append(proportions[0])
            siliconFraction.append(proportions[1])
            metalFraction.append(proportions[2])
            hydrogenFraction.append(proportions[3])
            saltFraction.append(proportions[4])

    cumulativeOxygenFraction = np.array(oxygenFraction)
    cumulativeSiliconFraction = (cumulativeOxygenFraction + np.array(siliconFraction))
    cumulativeMetalFraction = (cumulativeSiliconFraction + np.array(metalFraction))
    cumulativeHydrogenFraction = (cumulativeMetalFraction + np.array(hydrogenFraction))
    cumulativeSaltFraction = (cumulativeHydrogenFraction + np.array(saltFraction))

    matplotlib.rcParams.update({'font.size': 18})
    plt.figure(figsize=(6, 6))
    plt.ylabel("Layer Boundary Radius [km]")
    plt.xlabel("Atomic Proportion [%]")

    plt.plot(cumulativeOxygenFraction * 100, np.array(feasibleRadius) / 1e3, label="O")
    plt.plot(cumulativeSiliconFraction * 100, np.array(feasibleRadius) / 1e3, label="Si")
    plt.plot(cumulativeMetalFraction * 100, np.array(feasibleRadius) / 1e3, label="Fe + Mg")
    plt.plot(cumulativeHydrogenFraction * 100, np.array(feasibleRadius) / 1e3, label="H")
    plt.plot(cumulativeSaltFraction * 100, np.array(feasibleRadius) / 1e3, label="Na + Cl")

    plt.xlim([0, 100])
    plt.legend()
    plt.tight_layout()
    plt.show()