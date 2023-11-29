# -*- coding: utf-8 -*-
"""
Created on Thu Nov  9 13:43:07 2023

@author: Usuario
"""

import os
import csv
import itertools

import time as pytime
import numpy as np

from Utils import taul_to_minutes, minutes_to_taul
from Predictive_model import system_evolution_nlg

# Time units: tau_l = N/k = 23.0 minutes

""" ------- PARAMETERS ------- """

delta_t = 0.01
values_Nb0 = [1000]
values_B = [1.0e8]
values_P = [1.0e9]
N = 60
k = 6.0
nd = 30
nu = 5.0e-10
values_n0 = [1.0e8]
Nimax_values = [90]
Nimax_std_values = [i/3.2 for i in Nimax_values]

# Collapse time in hours
collapse_time_hours = 4.0

# time parameters in minutes
doubling_time_mins = 20.0
time_window_mins_values = [2.0]

# Time parameters in tau_l units
doubling_time = minutes_to_taul(doubling_time_mins, N, k)
time_window_values = [minutes_to_taul(i, N, k) for i in time_window_mins_values]
collapse_time = minutes_to_taul((collapse_time_hours*60.0),N,k)

""" ---------- WE FETCH THE NUMBER OF SIMULATION -------"""

with open('simulation_number.txt', 'r') as file:
    counter = int(file.read())

# Increment the counter
counter += 1

with open('simulation_number.txt', 'w') as file:
    file.write(str(counter))
    
storage_directory = 'Simulation_'+str(counter)

folder_path = os.path.join(os.getcwd(),'RESULTS',storage_directory)
if not os.path.exists(folder_path):
    os.mkdir(folder_path,0o666)


"""------------ NUMBA COMPILATION -------------"""
# We call the function with just one iteration so that numba does the compilation
niter = 1
_, _, _, _, _, _, _, _, _, _, _, _, _ = system_evolution_nlg(delta_t, values_Nb0[0], values_B[0], 
                                                    values_P[0], N, k, nd, nu, values_n0[0], doubling_time, 
                                                    time_window_values[0], Nimax_values[0], Nimax_std_values[0], collapse_time, niter)

print('Numba compilation done')


""" ---------- ACTUAL SIMULATION ------------- """
niter = 20000
start = pytime.time()

param_combinations = list(itertools.product(values_Nb0,values_B, values_P, values_n0, Nimax_values, Nimax_std_values, time_window_values))
ncomb = 0

for combination in param_combinations:
    
    Nb0,B,P,n0,Nimax,Nimax_std, time_window = combination

    time_vector, nphages, nhbacteria, nibacteria, ndbacteria, nnutrients, avg_ninfect, burstmd, burstsizemd, firstinf, lysist, ninfect, collstates = system_evolution_nlg(delta_t, Nb0, 
                                                                                                                                        B, P, N, k, nd, nu, n0, doubling_time, 
                                                                                                                                        time_window, Nimax, Nimax_std, collapse_time, niter, 
                                                                                                                                        constant_nutrients=False, allow_reinfection=True, time_window_track=True)
    total_time = pytime.time()-start
    

    time_hours = [taul_to_minutes(t, N, k)/60.0 for t in time_vector]

    """ --------- STORAGE OF THE DATA IN A FILE --------- """

    file_name = 'Simulation'+'_'+str(counter)+'_Combination_'+str(ncomb)+'.csv'
    complete_path = os.path.join(folder_path,file_name)

    with open(complete_path,'w',newline='') as csvfile:
    
        f = csv.writer(csvfile,delimiter=',',quotechar='"', quoting=csv.QUOTE_MINIMAL)
    

        f.writerow(["PARAMETERS"])
        f.writerow(['Parameter','Name in code', 'Value','Units'])
        f.writerow(['Delta t','delta_t',delta_t, 'tau_l units'])
        f.writerow(['Initial number of bacteria','Nb0',Nb0, '-'])
        f.writerow(['Initial concentration of bacteria','B',B, 'cells/ml'])
        f.writerow(['Initial concentration of phage', 'P',P, 'cells/ml'])
        f.writerow(['Maximum number of timer molecules','N', N, '-'])
        f.writerow(['Production rate of timer molecules','k', k, 'cells/min·ml'])
        f.writerow(['Number of dead timer molecules when superinfected','nd', nd, '-'])
        f.writerow(['Adsorption rate of phage to bacteria','nu', nu, 'cells/min·ml'])
        f.writerow(['Initial concentration of nutrients','n0', n0, 'nutrient/ml'])
        f.writerow(['Maximum infections one bacterium can handle','Nimax',Nimax,'-'])
        f.writerow(['Width of the maximum infection one bacterium can handle','Nimax_std',Nimax_std,'-'])
        f.writerow(['Time window to track maximum infections','time_window', time_window,'tau_l units'])
        
        f.writerow([])
        f.writerow([])

        f.writerow(['ADDITIONAL DATA FROM SIMULATION'])
        f.writerow(['Magnitude', 'Value','Unit'])
        f.writerow(['Duration', total_time, 's' ])
        f.writerow(['Volume of the system', float(Nb0)/float(B), 'ml'])
        f.writerow(['Number of iterations',niter,'-'])

        f.writerow([])
        f.writerow([])

        f.writerow(['DATA FROM SIMULATION'])
        f.writerow(['Time (hours)','Number of phage','Number of healthy bacteria','Number of infected bacteria','Number of dead bacteria', 'Number of nutrients', 'Average MOSI', 'Burst events due to MD'])
        f.writerow(['---'])
    
        for i in range(0,niter):
            f.writerow([time_hours[i],nphages[i],nhbacteria[i],nibacteria[i],ndbacteria[i],nnutrients[i], avg_ninfect[i], burstmd[i], burstsizemd[i]])
        
        f.writerow(['---'])
        f.writerow(['Time of first infection(mins)','Tau (mins)', 'MOSI'])
        f.writerow(['---'])
        for i in range(1, len(lysist)):
            f.writerow([firstinf[i],lysist[i],ninfect[i]])
        
        f.writerow(['---'])
        f.writerow(['Last infections at the time of collapse'])
        f.writerow(['---'])
        for i in range(0,len(collstates)):
            f.writerow([collstates[i]])
            
        

        ncomb += 1
        print('Combination',str(ncomb),'done')
        
print('Code lasted in', total_time, 'seconds')
print('Simulation id:',counter)

















