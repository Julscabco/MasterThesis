# -*- coding: utf-8 -*-
"""
Created on Tue Sep 26 16:52:55 2023

@author: Usuario
"""

import os

import numpy as np
import matplotlib.pyplot as plt


def get_burst_size(t, tlim, a):
    if t < tlim:
        return 0.0
    elif t >= tlim:
        return a * t - a*tlim

# The plot
fig_hist, ax_hist = plt.subplots(figsize=(10, 6))
file_directory = os.path.join(os.getcwd(),os.path.join('Figure1ForOralPresentation','Prova'))
save_plot = True
burst_histogram = False
p = 0.15
c = 0.1

for filename in os.listdir(file_directory):
    #Event-based Gillespie algorithm
    #n,nd,k,mosi,ti_superinfection,len_superinfection,ndata = filename.split('txt')[0].split('_')
    #output_folder='Corrected exponential'
    
    # Constant time interval
    delta_t,nb,n,k,nd,mosi,ti_superinfection,len_superinfection = filename.split('txt')[0].split('_')
    output_folder='Figure1ForOralPresentation'
    print(n)
                      
    f = open(os.path.join(file_directory,filename),"r")
    tau_vector=[]
    for line in f.readlines():
        tau_vector.append(float(line.split('\n')[0]))


    # Calculation and normamlization of the histogram
    hist, bins = np.histogram(tau_vector, bins=40)
    bins = (bins[1:] + bins[:-1])/2.0
    binwidth = bins[2]-bins[1]
    
    if burst_histogram == True:
        hist = [hist[i]*get_burst_size(bins[i],15.0,10.0) for i in range(0,len(hist))]
    hist = hist/(len(tau_vector)*binwidth)
    
    

    ax_hist.scatter(bins, hist, label= r'${N}_{T}$=' + str(n) + ', ' + r'k='+str(k))

    """
    ax_hist.set_title(r'MOSI='+str(mosi) + ' ,' + r'${\tau}_{0s}=$'+ str(ti_superinfection) +
                      ' ,' + r'${\tau}_{s}=$'+ str(len_superinfection)  + ','+
                      r'${n}_{T}$=' + str(n)
                      + ',' + 'k='+str(k), fontsize=20)
    """
    if burst_histogram ==False:
        ax_hist.set_title(r'MOSI='+str(mosi), fontsize=20)
    elif burst_histogram==True:
        ax_hist.set_title(r'MOSI='+str(mosi), fontsize=20)

    ax_hist.set_xlim(0.0,50.0)
    if burst_histogram==True:
        ax_hist.set_xlabel(r'${\tau}_{l}$ (min)', fontsize=18, labelpad=8)
        ax_hist.set_ylabel('Proportional to burst phages', fontsize=18, labelpad=8)
    elif burst_histogram==False:
        ax_hist.set_xlabel(r'${\tau}_{l}$ (min)',fontsize=18, labelpad=8)
        ax_hist.set_ylabel('Normalized Frequency', fontsize=18, labelpad=8)
    ax_hist.tick_params(axis='both', which='major', labelsize=20)
    ax_hist.legend(loc='best', fontsize=18)
    


if save_plot==True:
    if burst_histogram==True:
        fig_name = 'Phage_'+str(nd)+'_'+str(k)+'_'+ str(mosi)+'_'+str(ti_superinfection)+'_'+str(len_superinfection)+'_'+str(p)+'.png'
    elif burst_histogram==False:
        fig_name = 'Constant_N'+'_'+str(n)+'_'+str(k)+'_'+ str(mosi)+'_'+str(ti_superinfection)+'_'+str(len_superinfection)+'.png'
        
    fig_hist.savefig(os.path.join(os.getcwd(),'parameter_exploring',output_folder,fig_name))
    