import math
import TwoLayerModels

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
    test2LayerDensity()