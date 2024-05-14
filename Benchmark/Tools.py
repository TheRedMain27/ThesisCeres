"""
This file contains some code that is used multiple times, to avoid duplicate code.
"""
import numpy as np
import matplotlib.pyplot as plt

# constants
G = 6.67430e-11
R = 6371e3  # [m]
M = 3.986004418e14 / G  # [kg]
g0 = 9.80665  # [m/s2]

topoDataType = np.dtype("int16")
topoDataType = topoDataType.newbyteorder(">")

topoFileName = "Data/Earth2014.BED2014.5min.geod.bin"  # (Hirt and Rexer, 2015)
maskFileName = "Data/MSK2014_landtypes.1min.geod.bin"  # (Hirt and Rexer, 2015)

meanContinentalThickness = 33.6  # [km], (Reguzzoni and Sampietro, 2015)
meanOceanicThickness = 12.9  # [km], (Reguzzoni and Sampietro, 2015)

meanContinentalDensity = 2.65e3  # [kg/m3], (Tewari, Prasad and Kumar, 2018)
meanOceanicDensity = 2.9e3  # [kg/m3], (Tewari, Prasad and Kumar, 2018)

compensationDepth = 60e3  # [m]
crustDensity = 2.81e3  # [kg/m3]
mantleDensity = 4.5e3  # [kg/m3]

mohoFileName = "Data/MohoCrustalThicknessSzwillusEtAl2019"  # (Szwillus et al., 2019)

def readTopo():
    # reads topography data from file and returns np array
    data = np.fromfile(topoFileName, dtype=topoDataType, count=-1)
    return np.flip(np.reshape(data, (2160, 4320)), 0)

def decreaseTopoResolution(topo, resolution):
    # decreases the resolution by taking the mean value of 2:1 rectangles of such size that
    # the desired resolution is achieved
    if (resolution[0] > topo.shape[0] or resolution[1] > topo.shape[1] or resolution[1] / resolution[0] != 2 or
            topo.shape[0] % resolution[0] > 1e6 or topo.shape[1] % resolution[1] > 1e6):
        raise ValueError("Original resolution (2160, 4320) should be divisible by new resolution")
    return (topo.reshape(resolution[0], int(topo.shape[0] / resolution[0]), resolution[1],
                         int(topo.shape[1] / resolution[1])).mean(axis=1).mean(axis=2))

def getLandMask():
    # reads the land mask
    data = np.fromfile(maskFileName, dtype=topoDataType, count=-1)
    data = np.flip(np.reshape(data, (10800, 21600)), 0)
    data[data != 2] = 1
    data[data == 2] = 0
    return data

def calculateAverageThickness():
    mask = getLandMask()
    mask[mask == 0] = meanOceanicThickness
    mask[mask == 1] = meanContinentalThickness
    meanThickness = np.mean(mask)
    print("The mean moho thickness is: " + str(meanThickness) + " [km].")

def calculateAverageDensity():
    mask = getLandMask()
    mask[mask == 0] = meanOceanicDensity
    mask[mask == 1] = meanContinentalDensity
    meanDensity = np.mean(mask)
    print("The mean crust density is: " + str(meanDensity) + " [kg/m3].")

def calculateLowerCrustGravity():
    # calculates the gravity at the crust-mantle boundary
    crustMass = (4 / 3) * np.pi * (R ** 3 - (R - compensationDepth) ** 3) * crustDensity
    residualMass = M - crustMass
    return G * residualMass / (R - compensationDepth) ** 2

def readMoho():
    moho = np.loadtxt(mohoFileName, skiprows=34, usecols=2)
    return np.reshape(moho, (180,360))

if __name__ == "__main__":
    moho = readMoho()
    topo = decreaseTopoResolution(readTopo(), moho.shape) / 1e3
    fig1 = plt.figure()
    plt.imshow(moho, cmap="hot")
    plt.colorbar()
    fig2 = plt.figure()
    plt.imshow(topo, cmap="hot")
    plt.colorbar()
    plt.show()