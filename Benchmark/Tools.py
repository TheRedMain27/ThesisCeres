"""
This file contains some code that is used multiple times, to avoid duplicate code.
"""
import numpy as np
import matplotlib.pyplot as plt

# constants
topoDataType = np.dtype("int16")
topoDataType = topoDataType.newbyteorder(">")

topoFileName = "Data/Earth2014.BED2014.1min.geod.bin" #  (Hirt and Rexer, 2015)

def readTopo():
    # reads topography data from file and returns np array
    data = np.fromfile(topoFileName, dtype=topoDataType, count=-1)
    return np.flip(np.reshape(data, (10800, 21600)), 0)