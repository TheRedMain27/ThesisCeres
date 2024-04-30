from Tools import *

dr = 1
rlist = np.arange(0, R + dr, dr)

outerDensity = 1750
innerDensity = 2495
innerRadius = R - 85.52 * 1e3

profile = np.ones(rlist.shape[0])
profile[:round(innerRadius / dr) + 1] *= innerDensity
profile[round(innerRadius / dr) + 1:] *= outerDensity
MVector, gVector, PVector = profileCalculator(rlist, profile, dr)

matplotlib.rcParams.update({'font.size': 18})
fig, axs = plt.subplots(1, 3, figsize=(15, 6))
fig.supylabel("Radius [km]")
axs[0].plot(profile, rlist / 1e3)
axs[0].set_xlabel("Density [kg/m$^{3}$]")
axs[1].plot(gVector, rlist / 1e3),
axs[1].set_xlabel("Gravity [m/s$^{2}$]")
axs[2].plot(PVector / 1e6, rlist / 1e3)
axs[2].set_xlabel("Pressure [MPa]")
plt.tight_layout()
plt.savefig(r"Images/FinalModel.pdf")
plt.savefig(r"Images/FinalModel.png")
plt.show()