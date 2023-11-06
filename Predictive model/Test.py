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

def nutrient_integral(gmax,B,n0,delta_t,tmax):

    niter = int(tmax/delta_t)
    time = np.zeros(niter)
    n = np.zeros(niter)
    b = np.zeros(niter)
    
    n[0] = n0
    time[0] = 0.0
    b[0] = B

    for ite in range(1,niter):
        n[ite] =  n[ite -1] - b[ite-1]*monod_growth_law(n[ite-1],n0/5.0,gmax)*delta_t
        b[ite] = b[ite-1] + monod_growth_law(n[ite-1],n0/5.0,gmax)*b[ite-1]*delta_t
        time[ite] = time[ite-1] + delta_t
        
    return time,n,b

    
    
t,nutrients,bacteria = nutrient_integral(gnmax,1.0e6,n0,0.01,70.0)

plt.figure()
plt.plot(t,nutrients)

plt.figure()
plt.plot(t,bacteria)