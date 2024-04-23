import math
import MaxPressure

def test2LayerDensity():
    testRadius = [1, 100000, 300000, 450000, 469999]
    for i, innerRadius in enumerate(testRadius):
        rhoInner, rhoOuter = MaxPressure.determine2LayerDensities(innerRadius)
        Mdiff = (MaxPressure.M - 4 / 3 * math.pi *
                 (rhoOuter * MaxPressure.R ** 3 + (rhoInner - rhoOuter) * innerRadius ** 3)) / MaxPressure.M
        Idiff = (MaxPressure.I - 8 / 15 * math.pi *
                 (rhoOuter * MaxPressure.R ** 5 + (rhoInner - rhoOuter) * innerRadius ** 5)) / MaxPressure.I
        print("test " + str(i + 1) + " : " + str((Mdiff < 1e-10 and Idiff < 1e-10)))

if __name__ == "__main__":
    test2LayerDensity()