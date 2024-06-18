from Tools import *

mantleDensity = 2000; # (King, 2018)
crustDensity = 910; # from optimization
referenceDepth = 11.7e3; # from optimization

mantleRadius = R - referenceDepth

def coreDensityDifference(coreRadius):
    coreDensity1 = (3 * M / (4 * np.pi) - crustDensity * (R ** 3 - mantleRadius ** 3) -
                    mantleDensity * (mantleRadius ** 3 - coreRadius ** 3)) / coreRadius ** 3
    coreDensity2 = (15 * I / (8 * np.pi) - crustDensity * (R ** 5 - mantleRadius ** 5) -
                     mantleDensity * (mantleRadius ** 5 - coreRadius ** 5 )) / coreRadius ** 5
    return coreDensity1 - coreDensity2

coreRadiusGuess = 100e3

coreRadius = scipy.optimize.fsolve(coreDensityDifference, coreRadiusGuess)
coreDensity = (3 * M / (4 * np.pi) - crustDensity * (R ** 3 - mantleRadius ** 3) -
                    mantleDensity * (mantleRadius ** 3 - coreRadius ** 3)) / coreRadius ** 3

print(coreRadius / 1e3, coreDensity)