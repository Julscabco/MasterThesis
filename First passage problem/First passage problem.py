import numpy as np
import matplotlib.pyplot as plt

def taul_to_minutes(t, n, k):
    return t * (float(k) / float(n)) * 23.0

# Parameters
k = 0.06
r = 0.0005 #original value 0.0005
n_values = [20,30,40,50]
N = 60
T = 3000 # Total simulation time
dt = 0.01  # Time step

# Number of spatial nodes
num_nodes = N + 1


fig, ax = plt.subplots(2,1,figsize=(10,10))

for n in n_values:
    
    # Initialize probability array
    P = np.zeros((num_nodes, int(T / dt) + 1))

    tt = np.array([0.0])

    # Initial condition
    P[0, 0] = 1.0
    # Forward Euler method
    for t in range(1, P.shape[1]):
        for i in range(0, num_nodes):
        
            new_timer_income = (k * P[i - 1, t - 1] if i > 0 else 0)
            infected_income = (r * P[i + n, t - 1] if N >(n + i) > n else 0)    
            new_timer_outcome = -(k * P[i, t - 1] if i < N else 0) 
            infected_outcome = -(r * P[i, t-1] if N > i > 0  else 0)
            
            if i == 0:
                extrap = r*np.sum(P[0:n,t-1])
            else:
                extrap = 0.0
            P[i, t] = P[i, t - 1] + dt *(new_timer_income+infected_income+new_timer_outcome+infected_outcome+extrap)      
    
        tt = np.append(tt,tt[t-1]+dt)
    

    # Calculate the first-passage time distribution from the CDF
    pdf_first_passage_time = np.diff(P[N,:]) / (np.diff(tt))
    
    time_values = [taul_to_minutes(i, 60, 6)/60.0 for i in tt]
    ax[0].plot(time_values[:-1], P[N,:-1], label='n='+str(n))
    ax[1].plot(time_values[:-1], pdf_first_passage_time, label='n = '+str(n))


ax[0].set_ylabel("Probability Density")
ax[0].set_title("P(N,t)")
ax[0].legend(loc='best')


ax[1].set_xlabel("Time (min)")
ax[1].set_ylabel("Probability Density")
ax[1].set_title("First Passage Time Distribution")
ax[1].legend(loc='best')

