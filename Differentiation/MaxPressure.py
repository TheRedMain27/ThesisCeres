"""
Several functions which perform calculations which are performed multiple times.
"""
import math
import numpy as np
import scipy.linalg

R = 470 * 10 ** 3  # [m], park 2016
M = 62.62736 * (10 ** 9) / (6.67430 * 10 ** -11)  # [kg], konopliv 2018
I = 0.37 * M * R ** 2  # park, 2016

def determine2LayerDensities(innerRadius):
    A = np.array([[float(innerRadius ** 3), float(R ** 3 - innerRadius ** 3)],
                   [float(innerRadius ** 5), float(R ** 5 - innerRadius ** 5)]])
    b = np.array([M / (4/3 * math.pi), I / (8/15 * math.pi)])
    x = scipy.linalg.solve(A, b)
    return x[0], x[1]

if __name__ == "__main__":
    pass
