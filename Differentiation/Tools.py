"""
This file contains some code that is used multiple times, to avoid duplicate code.
"""

def profileCalculator(profile, dr):
    # calculates mass, gravity and pressure for a density profile
    # profile is a

    # set up values of r at which variables will be calculated
    rlist = np.arange(0, R, dr)

    # set up rlists for mass calculation
    # mantle gets an additional point to calculate the first point of the crust
    rlistMantle = np.arange(0, innerRadius + dr, dr)
    rlistCrust = np.arange(innerRadius, R, dr)

    # initialize empty lists to hold calculated values
    MlistMantle = np.zeros((rlistMantle.shape[0]))
    MlistCrust = np.zeros((rlistCrust.shape[0]))

    # use Euler integrator to calculate mantle masses
    for i, r in enumerate(rlistMantle[:-1]):
        dMdrMantle = 4 * np.pi * innerDensity * r ** 2
        MlistMantle[i + 1] = MlistMantle[i] + dMdrMantle * dr

    # set first point of crust using mantle
    MlistCrust[0] = MlistMantle[-1]
    # delete last point of mantle (it will be included in crust)
    MlistMantle = MlistMantle[:-1]

    # use Euler integrator to calculate crust masses
    for i, r in enumerate(rlistCrust[:-1]):
        dMdrCrust = 4 * np.pi * outerDensity * r ** 2
        MlistCrust[i + 1] = MlistCrust[i] + dMdrCrust * dr

    # tie together mass lists and scale them
    Mlist = np.concatenate((MlistMantle, MlistCrust))
    Mlist1020 = Mlist / 1e20

    # calculate gravitational acceleration
    g = (G * Mlist) / rlist ** 2
    # correct for division by zero
    g[0] = 0

    # set up rlists for pressure calculation
    # crust gets an additional point to calculate the first point of the mantle
    rlistMantle = np.arange(0, innerRadius, dr)
    rlistCrust = np.arange(innerRadius - dr, R, dr)

    # initialize empty lists to hold calculated values
    plistMantle = np.zeros((rlistMantle.shape[0]))
    plistCrust = np.zeros((rlistCrust.shape[0]))

    # use Euler integrator to calculate crust pressures
    for i, r in enumerate(np.flip(rlistCrust)[:-1]):
        dpdrCrust = G * 4 / 3 * np.pi * outerDensity ** 2 * r
        plistCrust[i + 1] = plistCrust[i] + dpdrCrust * dr

    # set first point of mantle using crust
    plistMantle[0] = plistCrust[-1]
    # delete last point of crust (it will be included in mantle)
    plistCrust = plistCrust[:-1]

    # use Euler integrator to calculate mantle pressures
    for i, r in enumerate(np.flip(rlistMantle)[:-1]):
        dpdrMantle = G * 4 / 3 * np.pi * innerDensity ** 2 * r
        plistMantle[i + 1] = plistMantle[i] + dpdrMantle * dr

    # tie together pressure lists and scale them
    plist = np.concatenate((plistCrust, plistMantle))
    plistMPa = np.flip(plist) / 1e6

    # get radius in km
    rlistkm = rlist / 1e3

    return rlistkm, Mlist1020, g, plistMPa