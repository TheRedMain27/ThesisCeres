from Tools import *

brineFraction = 0
rockFraction = 0

coreDensity = ironDensity
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

coreRadiusGuess = 50e3

coreRadius = scipy.optimize.fsolve(mantleRadiusDifference, coreRadiusGuess)
mantleRadius = ((3 * M / (4 * np.pi) - crustDensity * R ** 3 -
                     (coreDensity - mantleDensity - crustDensity) * coreRadius ** 3)
                     / (mantleDensity - crustDensity))
mantleRadius = np.sign(mantleRadius) * (np.abs(mantleRadius)) ** (1 / 3)

print(coreRadius / 1e3, mantleRadius / 1e3)

coreVolume = (4/3) * np.pi * mantleRadius ** 3
coreMass = coreVolume * coreDensity
mantleVolume = (4/3) * np.pi * mantleRadius ** 3
crustVolume = (4/3) * np.pi * R ** 3 - mantleVolume
crustMass = crustVolume * crustDensity

ironMass = coreMass
serpentiniteMass = (1 - brineFraction) * mantleVolume * serpentiniteDensity
brineMass = brineFraction * mantleVolume * brineDensity
waterMass = (waterDensity / brineDensity) * brineMass
saltMass = (dissolvedSalt / brineDensity) * brineMass
iceMass = crustMass

ironMolecules = ironMass / ironMolarMass # [kmol]
serpentiniteMolecules = serpentiniteMass / serpentiniteMolarMass
waterMolecules = (iceMass + waterMass) / waterMolarMass
saltMolecules = saltMass / saltMolarMass

ironAtoms = ironMolecules
oxygenAtoms = 9 * serpentiniteMolecules + waterMolecules
siliconAtoms = 2 * serpentiniteMolecules
magnesiumAtoms = 3 * serpentiniteMolecules
hydrogenAtoms = 2 * waterMolecules
saltAtoms = 2 * saltMolecules

totalAtoms = ironAtoms + oxygenAtoms + siliconAtoms + magnesiumAtoms + hydrogenAtoms + saltAtoms

print(np.array([ironAtoms, oxygenAtoms, siliconAtoms, magnesiumAtoms, hydrogenAtoms, saltAtoms]) / totalAtoms * 100)