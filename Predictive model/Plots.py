# -*- coding: utf-8 -*-
"""
Created on Thu Nov  9 13:20:48 2023

@author: Usuario
"""


""" This code has to read the data stored fromt he simulations 
and generate the desired plots """

import os
import csv

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize

from Utils import nutrient_integral


simulation_id = 59
folder_name = 'Simulation'+'_'+str(simulation_id)
storage_directory = os.path.join(os.getcwd(),os.path.join('RESULTS',folder_name))

destiny_directory = os.path.join(os.getcwd(),os.path.join('FIGURES',folder_name))
if not os.path.exists(destiny_directory):
    os.mkdir(destiny_directory,0o666)
   

""" ----- PLOT 1 CONFIG ----- """
# Plot of number of phage, bacteria or nutrients vs time
plot_1 = True
save_plot_1 = False

parameters_in_title_1 = ['Nb0','B','P','n0']
parameters_symbols_1 = [r'${Nb}_{0}$', r'${B}_{0}=$',r'${P}_{0}=$',r'${N}_{0}=$']




""" ------- PLOT 2 CONFIG ------ """
# Plot of bacteria and nutrients vs time compared to theoretical integration
plot_2 =    False
save_plot_2 = False

parameters_in_title_2 = ['B','n0']
parameters_symbols_2 = [r'${B}_{0}$',r'${N}_{0}$']

""" ------ PLOT 3 CONFIG ------- """
# Histogram of lysis times
plot_3 = True
save_plot_3 = False

""" ------ PLOT 4 CONFIG ------- """

# Histogram of superinfections
plot_4 = True
save_plot_4 = False

""" ------ PLOT 5 CONFIG ------ """

plot_5 = True
save_plot_5 = True

parameters_in_title_5 = ['B','P','n0']
parameters_symbols_5 = [r'${B}_{0}=$',r'${P}_{0}=$',r'${N}_{0}=$']


""" ------- FILE READING --------"""

for filename in os.listdir(storage_directory):
    complete_path = os.path.join(storage_directory,filename)
    ncomb = (filename.split('_')[-1]).split('.')[0]

    with open(complete_path,"r") as csvfile:
        csvreader = csv.reader(csvfile)
    
        for row in csvreader:
            if 'PARAMETERS' in row:
                break
    
        for row in csvreader:
            column_names = row
            break
    
        df_parameters= pd.DataFrame(columns=column_names)
        i = 0
        for row in csvreader:
            if len(row)>2:
                df_parameters.loc[i] = row
                i = i + 1
            else:
                break

        for row in csvreader:
            if 'ADDITIONAL DATA FROM SIMULATION' in row:
                break
        
        for row in csvreader:
            column_names = row
            break
        
        df_add_parameters = pd.DataFrame(columns = column_names)
    
        i = 0
        for row in csvreader:
            if len(row)>2:
                df_add_parameters.loc[i] = row
                i = i + 1
            else:
                break
    
    
        time = []
        nphages = []
        nhbacteria = []
        nibacteria = []
        ndbacteria = []
        nnutrients = []
        firstinfect = []
        lysist = []
        ninfect = []
    

        for row in csvreader:
            if '---' in row:
                break
        
        for row in csvreader:
            if '---' in row:
                break
            else:
                t,p,hb,ib,db,n = row
                time.append(float(t))
                nphages.append(float(p))
                nhbacteria.append(float(hb))
                nibacteria.append(float(ib))
                ndbacteria.append(float(db))
                nnutrients.append(float(n))

        
        for row in csvreader:
            if '---' in row:
                break
            
            
        for row in csvreader:
            if '---' in row:
                break
            else:
                fi,l,inf = row
                lysist.append(float(l))
                ninfect.append(int(inf))
                firstinfect.append(float(fi))
        


    time = np.array(time)
    nphages = np.array(nphages)
    nhbacteria = np.array(nhbacteria)
    nibacteria = np.array(nibacteria)
    ndbacteria = np.array(ndbacteria)
    nnutrients = np.array(nnutrients)
    firstinfect = np.array(firstinfect)
    lysist = np.array(lysist)
    ninfect = np.array(ninfect)


"""---------PLOT 1: PHAGE, BACTERIA AND NUTRIENTS VS TIME --------"""

if plot_1 == True  :
    
    fig, ax = plt.subplots(3, 1, figsize=(10, 10))
    
    title = ''
    for i in range(0,len(parameters_in_title_1)):
        title += str(parameters_symbols_1[i])
        value = np.array(df_parameters['Value'][df_parameters['Name in code']==parameters_in_title_1[i]])[0]
        title += "%.2e"%float(value)
        title += ', '
            
    nbacteria = np.add(nhbacteria,nibacteria)

    ax[0].plot(time, nphages, label='Phage')
    ax[0].set_ylabel('Number of phage',fontsize=16)
    ax[0].set_title(title)

    ax[0].set_xlim(0.0,time[-1])

    ax[1].plot(time, nbacteria, label='Bacteria', color='orange')


    ax[1].set_ylabel('Number of bacteria',fontsize=16)

    ax[1].set_xlim(0.0,time[-1])

    ax[2].plot(time, nnutrients)
    ax[2].set_xlabel('Time (h)',fontsize=20)
    ax[2].set_ylabel('Number of nutrients',fontsize=16)
    ax[2].set_xlim(0.0,time[-1])
        
    for axis in ax:
        axis.tick_params(axis='y', labelsize=14)
    
    if save_plot_1 == True:
        figure_name = 'P1_S'+str(simulation_id)+'_comb_'+str(ncomb)+'.png'
        fig.savefig(os.path.join(destiny_directory, figure_name))
        
        print("Plot 1 saved as:",figure_name)
    

    





""" PLOT 2: BACTERIA AND NUTRIENTS VS TIME (SIMULATION VS INTEGRATION COMPARISON) """

if plot_2 == True:
    
    
    B = float(df_parameters.at[(df_parameters['Name in code']=='B').idxmax(),'Value'])
    n0 = float(df_parameters.at[(df_parameters['Name in code']=='n0').idxmax(),'Value'])
    V = float(df_add_parameters.at[(df_add_parameters['Magnitude']=='Volume of the system').idxmax(),'Value'])
    
    # We integrate the differential equations that the bacteria and virus should follow
    t,nutrients,bacteria = nutrient_integral(0.034,float(B),float(n0),(time[1]-time[0])*60.0,time[-1]*60.0)
    
    nbacteria = np.add(nhbacteria,nibacteria)
    t = [j/60.0 for j in t]
    nbacteria_c = [j/V for j in nbacteria]
    nnutrients_c = [j/V for j in nnutrients]
    
    title = ''
    for i in range(0,len(parameters_in_title_2)):
        title += (str(parameters_symbols_2[i]) + '=')
        value = np.array(df_parameters['Value'][df_parameters['Name in code']==parameters_in_title_2[i]])[0]
        title += "%1.1e"%float(value)
        title += ', '
    
    
    fig2, ax2 = plt.subplots()
    
    ax2.plot(t, bacteria, label = 'B(t)',color='black',linestyle='dashed')
    ax2.plot(time, nbacteria_c, label='Simulation',color = 'black')
    ax2.plot(t, nutrients, label = 'n(t)',color = 'red',linestyle='dashed')
    ax2.plot(time, nnutrients_c, label='Simulation', color = 'red')
    
    ax2.set_xlabel('time (mins)')
    ax2.set_ylabel('Concentration (ml-1)')
    ax2.set_title(title)
    ax2.legend(loc='best')
    
    if save_plot_2 == True:
        figure_name = 'P2_S'+str(simulation_id)+'_comb_'+str(ncomb)+'.png'
        fig2.savefig(os.path.join(destiny_directory, figure_name))
        
        print("Plot 2 saved as:",figure_name)
    


""" PLOT 3: HISTOGRAM OF LYSIS TIMES OF ALREADY DEAD CELLS """

if plot_3 == True:
    
    lysis_times_of_dead_cells = lysist[np.where(lysist!=0.0)]
    
    fig_hist1, ax_hist = plt.subplots(figsize=(10, 6))
    
    # Calculation and normamlization of the histogram
    hist, bins = np.histogram(lysis_times_of_dead_cells, bins=50)
    bins = (bins[1:] + bins[:-1])/2.0
    binwidth = bins[2]-bins[1]
    hist = hist/(len(lysist)*binwidth)
    
    
    ax_hist.scatter(bins, hist)
    ax_hist.set_xlim([0.0,200.0])
    ax_hist.set_xlabel('Time from infection to lysis (mins)',fontsize=16)
    ax_hist.set_ylabel('Normalized frequency',fontsize=16)
    plt.tick_params(axis='x', labelsize=14)
    plt.tick_params(axis='y', labelsize=14)
    
    if save_plot_3 == True:
        figure_name = 'P3_S'+str(simulation_id)+'_comb_'+str(ncomb)+'.png'
        fig_hist1.savefig(os.path.join(destiny_directory, figure_name))
    
        print("Plot 3 saved as:",figure_name)
    
    
    
""" PLOT 4: HISTOGRAM OF MOSI """

if plot_4 == True:
    
    
    fig_hist1, ax_hist = plt.subplots(figsize=(10, 6))
    
    # Calculation and normamlization of the histogram
    hist, bins = np.histogram(ninfect, bins=50)
    bins = (bins[1:] + bins[:-1])/2.0
    binwidth = bins[2]-bins[1]
    hist = hist/(len(ninfect)*binwidth)
    
    ax_hist.scatter(bins, hist)
    ax_hist.set_xlim([0.0,20.0])
    ax_hist.set_xlabel('MOSI',fontsize=16)
    ax_hist.set_ylabel('Normalized frequency',fontsize=16)
    plt.tick_params(axis='x', labelsize=14)
    plt.tick_params(axis='y', labelsize=14)
    
    
    
""" PLOT 5: HISTOGRAM OF LYSIS TIMES ACCORDING TO MOSI """  

if plot_5 == True:
    
    dead_cells_indices = np.where(lysist!=0.0)
    
    first_infections = firstinfect[dead_cells_indices]
    norm = Normalize(vmin=np.min(first_infections),vmax=np.max(first_infections))
    
    title = ''
    for i in range(0,len(parameters_in_title_5)):
        title += str(parameters_symbols_5[i])
        value = np.array(df_parameters['Value'][df_parameters['Name in code']==parameters_in_title_5[i]])[0]
        title += "%.2e"%float(value)
        title += ', '
    
    
    fig5, ax5 = plt.subplots()
    scatter = ax5.scatter(lysist[dead_cells_indices],ninfect[dead_cells_indices],
                c=first_infections,cmap='plasma',norm=norm,marker='o',alpha=0.7)
    ax5.set_xlabel('Time from first infection until lysis (min)',fontsize=12)
    ax5.set_ylabel('Multiplicity of superinfection (MOSI)',fontsize=12)
    ax5.set_title(title)
    
    cbar = plt.colorbar(scatter, ax=ax5)
    cbar.set_label('Time of Initial Infection')
    
    
    
    if save_plot_5 == True:
        figure_name = 'P5_S'+str(simulation_id)+'_comb_'+str(ncomb)+'.png'
        fig5.savefig(os.path.join(destiny_directory, figure_name))
    
        print("Plot 5 saved as:",figure_name)
    

    
            
    
    
    
    

    
    
    
                
        


        


    


    



























