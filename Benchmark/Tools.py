"""
This file contains some code that is used multiple times, to avoid duplicate code.
"""
import numpy as np
import matplotlib.pyplot as plt

# constants
topoDataType = np.dtype("int16")
topoDataType = topoDataType.newbyteorder(">")

topoFileName = "Data/Earth2014.BED2014.1min.geod.bin"  # (Hirt and Rexer, 2015)
maskFileName = "Data/MSK2014_landtypes.1min.geod.bin"  # (Hirt and Rexer, 2015)

meanContinental = 33.6  # [km], (Reguzzoni and Sampietro, 2015)
meanOceanic = 12.9  # [km], (Reguzzoni and Sampietro, 2015)

compensationDepth = 19.2e3  # [m]

def readTopo():
    # reads topography data from file and returns np array
    data = np.fromfile(topoFileName, dtype=topoDataType, count=-1)
    return np.flip(np.reshape(data, (10800, 21600)), 0)

def getLandMask():
    # reads the land mask
    data = np.fromfile(maskFileName, dtype=topoDataType, count=-1)
    data = np.flip(np.reshape(data, (10800, 21600)), 0)
    data[data != 2] = 1
    data[data == 2] = 0
    return data

def calculateAverageThickness():
    mask = getLandMask()
    mask[mask == 0] = meanOceanic
    mask[mask == 1] = meanContinental
    meanThickness = np.mean(mask)
    print("The mean moho thickness is: " + str(meanThickness) + " [km].")

if __name__ == "__main__":
    calculateAverageThickness()