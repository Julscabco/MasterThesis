import matplotlib.pyplot as plt
import numpy as np

def monod_growth_law(c, kn, gnmax=0.034):
    return gnmax * c / (kn + c)

gnmax = 0.079
n0 = 1.0e9
kn = n0/5.0

n = np.linspace(0.0,1.0e9,1000)

g = [monod_growth_law(i,kn,gnmax) for i in n]
plt.plot(n,g)
plt.show()