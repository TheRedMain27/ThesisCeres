from Tools import *

def determineLinearCoefficients():
    # solves linear density profile constrained by mass and moment of inertia
    # returns coefficients [c0, c1], where rho = c0 - c1 * r
    A = np.array([[1/3, -R/4],
                  [1/5, -R/6]])
    b = np.array([M / (4 * np.pi * R ** 3),
                  I / (4 * np.pi * R ** 5)])
    return scipy.linalg.solve(A, b)
def determineQuadraticCoefficients(c0):
    # solves quadratic density profile constrained by mass and moment of inertia
    # returns coefficients [c1, c2], given c0, where rho = c0 + c1 * r + c2 * r^2
    A = np.array([[R/4, (R**2)/5],
                  [R/6, (R**2)/7]])
    b = np.array([M / (4 * np.pi * R ** 3) - c0 / 3,
                  I / (4 * np.pi * R ** 5) - c0 / 5])
    return scipy.linalg.solve(A, b)

def determineCubicCoefficients(c0):
    # solves cubic density profile constrained by mass and moment of inertia
    # returns coefficients [c2, c3], given c0, where rho = c0 + c2 * r^2 + c3 * r^3
    # assumption: slope of density at r = 0 is zero --> c1 = 0
    A = np.array([[(R**2)/5, (R**3)/6],
                  [(R**2)/7, (R**3)/8]])
    b = np.array([M / (4 * np.pi * R ** 3) - c0 / 3,
                  I / (4 * np.pi * R ** 5) - c0 / 5])
    return scipy.linalg.solve(A, b)

def determine5thOrderCoefficients():
    # solves 5th order polynomial density profile constrained by mass and moment of inertia
    # returns coefficients [c0, c3, c4, c5], where rho = c0 + c3 * r^3 + c4 * r^4 + c5 * r^5
    # assumptions: rho'(0) = rho'(R) = rho"(0) = rho"(R) = 0
    A = np.array([[0, 3, 4*R, 5*R**2],
                  [0, 6, 12*R, 20*R**2],
                  [1/3, (R**3)/6, (R**4)/7, (R**5)/8],
                  [1/5, (R**3)/8, (R**4)/9, (R**5)/10]])
    b = np.array([0,
                  0,
                  M / (4 * np.pi * R ** 3),
                  I / (4 * np.pi * R ** 5)])
    return scipy.linalg.solve(A, b)

def determine6thOrderCoefficients(c0):
    # solves 6th order polynomial density profile constrained by mass and moment of inertia
    # returns coefficients [c3, c4, c5, c6], given c0, where rho = c0 + c3 * r^3 + c4 * r^4 + c5 * r^5 + c6 * r^6
    # assumptions: rho'(0) = rho'(R) = rho"(0) = rho"(R) = 0
    A = np.array([[3, 4*R, 5*R**2, 6*R**3],
                  [6, 12*R, 20*R**2, 30*R**3],
                  [(R**3)/6, (R**4)/7, (R**5)/8, (R**6)/9],
                  [(R**3)/8, (R**4)/9, (R**5)/10, (R**6)/11]])
    b = np.array([0,
                  0,
                  M / (4 * np.pi * R ** 3) - c0/3,
                  I / (4 * np.pi * R ** 5) - c0/5])
    return scipy.linalg.solve(A, b)

def determine7thOrderCoefficients(c0, c8):
    # solves 7th order polynomial density profile constrained by mass and moment of inertia
    # returns coefficients [c3, c4, c5, c6, c7], given c0,
    # where rho = c0 + c3 * r^3 + c4 * r^4 + c5 * r^5 + c6 * r^6 + c7 * r^7
    # assumptions: rho'(0) = rho'(R) = rho"(0) = rho"(R) = 0
    A = np.array([[3, 4*R, 5*R**2, 6*R**3, 7*R**4],
                  [6, 12*R, 20*R**2, 30*R**3, 42*R**4],
                  [R**3, R**4, R**5, R**6, R**7],
                  [(R**3)/6, (R**4)/7, (R**5)/8, (R**6)/9, (R**7)/10],
                  [(R**3)/8, (R**4)/9, (R**5)/10, (R**6)/11, (R**7)/12]])
    b = np.array([0,
                  0,
                  c8 - c0,
                  M / (4 * np.pi * R ** 3) - c0/3,
                  I / (4 * np.pi * R ** 5) - c0/5])
    return scipy.linalg.solve(A, b)

if __name__ == "__main__":
    dr = 1
    rlist = np.arange(0, R + dr, dr)

    # matplotlib.rcParams.update({'font.size': 16})
    # fig, axs = plt.subplots(1, 3, figsize=(15, 6))
    # fig.supylabel("Radius [km]")
    # fig.supxlabel("Density [kg/m$^{3}$]")
    #
    # linearCoefficients = determineLinearCoefficients()
    # c0 = linearCoefficients[0]
    # c1 = linearCoefficients[1]
    # profile = c0 - c1 * rlist
    # axs[0].plot(profile, rlist / 1e3)
    # axs[0].set_title("Linear")
    #
    # c0list = np.arange(1000, 5000, 1000)
    # for c0 in c0list:
    #     coefficients = determineQuadraticCoefficients(c0)
    #     c1 = coefficients[0]
    #     c2 = coefficients[1]
    #
    #     profile = c0 + c1 * rlist + c2 * rlist ** 2
    #
    #     axs[1].plot(profile, rlist/1e3, label = r"$\rho_{core}$ = " + str(c0))
    # axs[1].legend()
    # axs[1].set_title("Quadratic")
    #
    # c0list = np.arange(1000, 5000, 1000)
    # for c0 in c0list:
    #     coefficients = determineCubicCoefficients(c0)
    #     c2 = coefficients[0]
    #     c3 = coefficients[1]
    #
    #     profile = c0 + c2 * rlist ** 2 + c3 * rlist ** 3
    #
    #     axs[2].plot(profile, rlist / 1e3, label=r"$\rho_{core}$ = " + str(c0))
    # axs[2].legend()
    # axs[2].set_title("Cubic")
    #
    # plt.tight_layout()
    # plt.savefig(r"Images/PolynomialModels.pdf")
    # plt.savefig(r"Images/PolynomialModels.png")
    # plt.show()

    coefficients = determine7thOrderCoefficients(7874, 910)
    c3 = coefficients[0]
    c4 = coefficients[1]
    c5 = coefficients[2]
    c6 = coefficients[3]
    c7 = coefficients[4]

    profile = 7874 + c3 * rlist ** 3 + c4 * rlist ** 4 + c5 * rlist ** 5 + c6 * rlist ** 6 + c7 * rlist ** 7

    matplotlib.rcParams.update({'font.size': 18})
    plt.figure(figsize=(6, 6))
    plt.xlabel("Density [kg/m$^{3}$]")
    plt.ylabel("Radius [km]")
    plt.plot(profile, rlist / 1e3)
    plt.tight_layout()
    plt.savefig(r"Images/7thOrderModel.pdf")
    plt.savefig(r"Images/7thOrderModel.png")
    plt.show()