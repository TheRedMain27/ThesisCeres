from Tools import *

brineFraction = 0.15
rockFraction = 0.15

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

coreRadius = scipy.optimize.fsolve(mantleRadiusDifference, coreRadiusGuess)
mantleRadius = ((3 * M / (4 * np.pi) - crustDensity * R ** 3 -
                     (coreDensity - mantleDensity - crustDensity) * coreRadius ** 3)
                     / (mantleDensity - crustDensity))
mantleRadius = np.sign(mantleRadius) * (np.abs(mantleRadius)) ** (1 / 3)

print(coreRadius / 1e3, mantleRadius / 1e3)

# for brineFraction in np.arange(0, 1.01, 0.01):
#     for rockFraction in np.arange(0, 1.01, 0.01):
#         coreDensity = olivineDensity
#         mantleDensity = ((1 - brineFraction) * serpentiniteDensity + brineFraction * brineDensity)
#         crustDensity = (1 - rockFraction) * iceDensity + rockFraction * serpentiniteDensity
#
#         coreRadius = scipy.optimize.fsolve(mantleRadiusDifference, coreRadiusGuess)
#         mantleRadius = ((3 * M / (4 * np.pi) - crustDensity * R ** 3 -
#                              (coreDensity - mantleDensity - crustDensity) * coreRadius ** 3)
#                              / (mantleDensity - crustDensity))
#         mantleRadius = np.sign(mantleRadius) * (np.abs(mantleRadius)) ** (1 / 3)
#
#         coreVolume = (4/3) * np.pi * coreRadius ** 3
#         coreMass = coreVolume * coreDensity
#         mantleVolume = (4/3) * np.pi * mantleRadius ** 3
#         crustVolume = (4/3) * np.pi * R ** 3 - mantleVolume
#
#         olivineMass = coreMass
#         serpentiniteMass = ((1 - brineFraction) * mantleVolume * serpentiniteDensity
#                             + rockFraction * crustVolume * serpentiniteDensity)
#         brineMass = brineFraction * mantleVolume * brineDensity
#         waterMass = (waterDensity / brineDensity) * brineMass
#         saltMass = (dissolvedSalt / brineDensity) * brineMass
#         iceMass = (1 - rockFraction) * crustVolume * iceDensity
#
#         olivineMolecules = olivineMass / olivineMolarMass # [kmol]
#         serpentiniteMolecules = serpentiniteMass / serpentiniteMolarMass
#         waterMolecules = (iceMass + waterMass) / waterMolarMass
#         saltMolecules = saltMass / saltMolarMass
#
#         oxygenAtoms = 4 * olivineMolecules + 9 * serpentiniteMolecules  # + 9 * serpentiniteMolecules + waterMolecules
#         siliconAtoms = 2 * serpentiniteMolecules + olivineMolecules
#         ironAtoms = 2 * olivineMolecules
#         magnesiumAtoms = 3 * serpentiniteMolecules
#         hydrogenAtoms = 2 * waterMolecules + 4 * serpentiniteMolecules
#         saltAtoms = 2 * saltMolecules
#
#         totalAtoms = (oxygenAtoms + siliconAtoms + ironAtoms + magnesiumAtoms) / 0.95  # + hydrogenAtoms + saltAtoms
#
#         metalFraction = ironAtoms[0] / magnesiumAtoms[0]
#         siliconFraction = siliconAtoms[0] / (ironAtoms[0] + magnesiumAtoms[0])
#
#         if 0 <= coreRadius <= mantleRadius and coreRadius <= mantleRadius <= R \
#             and 0.8 < metalFraction < 1.2 and 0.2 < siliconFraction < 0.5:
#             print("Brine Fraction: " + str(brineFraction))
#             print("Rock Fraction: " + str(rockFraction))
#             print("Core radius: " + str(round(coreRadius[0] / 1e3)) + " [km]")
#             print("Mantle radius: " + str(round(mantleRadius[0] / 1e3)) + " [km]")
#             print("Iron-Magnesium Fraction: " + str(round(ironAtoms[0] / magnesiumAtoms[0], 3)))
#             print("Silicon-Metals Fraction: " + str(round(siliconAtoms[0] / (ironAtoms[0] + magnesiumAtoms[0]), 3)))


print("Oxygen : " + str(round(oxygenAtoms[0] / totalAtoms[0] * 100, 3)) + "%")
print("Silicon : " + str(round(siliconAtoms[0] / totalAtoms[0] * 100, 3)) + "%")
print("Iron : " + str(round(ironAtoms[0] / totalAtoms[0] * 100, 3)) + "%")
print("Magnesium : " + str(round(magnesiumAtoms[0] / totalAtoms[0] * 100, 3)) + "%")
# print("Hydrogen : " + str(round(hydrogenAtoms[0] / totalAtoms[0] * 100, 3)) + "%")
# print("Salt Atoms: " + str(round(saltAtoms[0] / totalAtoms[0] * 100, 3)) + "%")

print("Iron-Magnesium Fraction: " + str(round(ironAtoms[0] / magnesiumAtoms[0], 3)))
print("Silicon-Metals Fraction: " + str(round(siliconAtoms[0] / (ironAtoms[0] + magnesiumAtoms[0]), 3)))