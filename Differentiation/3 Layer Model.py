from Tools import *

brineFraction = 0.2
rockFraction = 0

coreDensity = olivineDensity
mantleDensity = ((1 - brineFraction) * serpentiniteDensity + brineFraction * brineDensity)
crustDensity = (1 - rockFraction) * iceDensity + rockFraction * serpentiniteDensity
def mantleRadiusDifference(coreRadius):
    mantleRadius1 = ((3 * M / (4 * np.pi) - crustDensity * R ** 3 -
                      (coreDensity - mantleDensity - crustDensity) * coreRadius ** 3)
                     / (mantleDensity - crustDensity))
    mantleRadius1 = np.sign(mantleRadius1) * (np.abs(mantleRadius1)) ** (1 / 3)
    mantleRadius2 = ((15 * I / (8 * np.pi) - crustDensity * R ** 5 -
                      (coreDensity - mantleDensity - crustDensity) * coreRadius ** 5)
                     / (mantleDensity - crustDensity))
    mantleRadius2 = np.sign(mantleRadius2) * (np.abs(mantleRadius2)) ** (1 / 5)
    return mantleRadius1 - mantleRadius2

coreRadiusGuess = 300e3

# for brineFraction in range(0, 0.15, 0.01):
#     for rockFraction in range(0, 0.6, 0.05):
#         coreDensity = olivineDensity
#         mantleDensity = ((1 - brineFraction) * serpentiniteDensity + brineFraction * brineDensity)
#         crustDensity = (1 - rockFraction) * iceDensity + rockFraction * serpentiniteDensity

coreRadius = scipy.optimize.fsolve(mantleRadiusDifference, coreRadiusGuess)
mantleRadius = ((3 * M / (4 * np.pi) - crustDensity * R ** 3 -
                     (coreDensity - mantleDensity - crustDensity) * coreRadius ** 3)
                     / (mantleDensity - crustDensity))
mantleRadius = np.sign(mantleRadius) * (np.abs(mantleRadius)) ** (1 / 3)

print("Core radius: " + str(round(coreRadius[0] / 1e3)) + " [km]")
print("Mantle radius: " + str(round(mantleRadius[0] / 1e3)) + " [km]")

coreVolume = (4/3) * np.pi * coreRadius ** 3
coreMass = coreVolume * coreDensity
mantleVolume = (4/3) * np.pi * mantleRadius ** 3
crustVolume = (4/3) * np.pi * R ** 3 - mantleVolume

olivineMass = coreMass
serpentiniteMass = ((1 - brineFraction) * mantleVolume * serpentiniteDensity
                    + rockFraction * crustVolume * serpentiniteDensity)
brineMass = brineFraction * mantleVolume * brineDensity
waterMass = (waterDensity / brineDensity) * brineMass
saltMass = (dissolvedSalt / brineDensity) * brineMass
iceMass = (1 - rockFraction) * crustVolume * iceDensity

olivineMolecules = olivineMass / olivineMolarMass # [kmol]
serpentiniteMolecules = serpentiniteMass / serpentiniteMolarMass
waterMolecules = (iceMass + waterMass) / waterMolarMass
saltMolecules = saltMass / saltMolarMass

oxygenAtoms = 4 * olivineMolecules + 9 * serpentiniteMolecules  # + 9 * serpentiniteMolecules + waterMolecules
siliconAtoms = 2 * serpentiniteMolecules + olivineMolecules
metalAtoms = 3 * serpentiniteMolecules + 2 * olivineMolecules
hydrogenAtoms = 2 * waterMolecules + 4 * serpentiniteMolecules
saltAtoms = 2 * saltMolecules

totalAtoms = (oxygenAtoms + siliconAtoms + metalAtoms) / 0.95  # + hydrogenAtoms + saltAtoms

print("Oxygen : " + str(round(oxygenAtoms[0] / totalAtoms[0] * 100, 3)) + "%")
print("Silicon : " + str(round(siliconAtoms[0] / totalAtoms[0] * 100, 3)) + "%")
print("Metals : " + str(round(metalAtoms[0] / totalAtoms[0] * 100, 3)) + "%")
# print("Hydrogen : " + str(round(hydrogenAtoms[0] / totalAtoms[0] * 100, 3)) + "%")
# print("Salt Atoms: " + str(round(saltAtoms[0] / totalAtoms[0] * 100, 3)) + "%")

print("Silicon-Metals Fraction: " + str(round(siliconAtoms[0] / metalAtoms[0], 3)))