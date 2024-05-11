from Tools import *

def equalMasses():
    topo = readTopo()

    extraDepth = topo * crustDensity / (mantleDensity - crustDensity)
    return compensationDepth + topo + extraDepth

def equalPressures():
    topo = readTopo()

    crustGravity = calculateLowerCrustGravity()

    extraDepth = topo * crustDensity / (mantleDensity - crustDensity) * (g0 / crustGravity)
    return compensationDepth + topo + extraDepth

if __name__ == "__main__":
    Airy = equalMasses()
    print(np.min(Airy))
    print(np.max(Airy))
    plt.imshow(Airy / 1e3, cmap="hot")
    plt.colorbar()
    plt.show()
    Airy.tofile("Results/AiryEqualMasses.txt")

    Airy = equalPressures()
    print(np.min(Airy))
    print(np.max(Airy))
    plt.imshow(Airy / 1e3, cmap="hot")
    plt.colorbar()
    plt.show()
    Airy.tofile("Results/AiryEqualPressures.txt")