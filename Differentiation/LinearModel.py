import math
import numpy as np
import scipy.linalg
import matplotlib
import matplotlib.pyplot as plt

G = 6.67430e-11
R = 470e3  # [m], park 2016
M = 62.62736e9 / G  # [kg], konopliv 2018
I = 0.37 * M * R ** 2  #[kg*m2] park, 2016

def determineCoefficients():
    # solves linear density profile constrained by mass and moment of inertia
    A = np.array([[1/3, -R/4], [1/5, -R/6]])
    b = np.array([M / (4 * math.pi * R ** 3), I / (4 * math.pi * R ** 5)])
    return scipy.linalg.solve(A, b)

if __name__ == "__main__":
    print(determineCoefficients())