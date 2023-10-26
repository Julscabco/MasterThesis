# -*- coding: utf-8 -*-
"""
Created on Tue Sep 26 16:52:55 2023

@author: Usuario
"""

import os

import numpy as np
import matplotlib.pyplot as plt

def get_burst_size(t, c=0.5, p=5.0):
    return c*np.exp(t*p)

# The plot
fig_hist, ax_hist = plt.subplots(figsize=(10, 6))
file_directory = os.path.join(os.getcwd(),os.path.join('ConstantTimeInterval','0.01_60_6_0.4_15_3_10000'))
save_plot = False
burst_histogram = True
p = 0.15
c = 0.1

for filename in os.listdir(file_directory):
    #Event-based Gillespie algorithm
    #n,nd,k,mosi,ti_superinfection,len_superinfection,ndata = filename.split('txt')[0].split('_')
    #output_folder='Corrected exponential'
    
    # Constant time interval
    delta_t,n,nd,k,mosi,ti_superinfection,len_superinfection,ndata = filename.split('txt')[0].split('_')
    output_folder='FIG2_ConstantTimeInterval'
    print(nd)
                      
    f = open(os.path.join(file_directory,filename),"r")
    tau_vector=[]
    for line in f.readlines():
        tau_vector.append(float(line.split('\n')[0]))
       
        
        

    # Calculation and normamlization of the histogram
    hist, bins = np.histogram(tau_vector, bins=60)
    bins = (bins[1:] + bins[:-1])/2.0
    binwidth = bins[2]-bins[1]
    
    if burst_histogram==True:
        hist = [hist[i]*get_burst_size(bins[i],c,p) for i in range(0,len(hist))]
    hist = hist/(len(tau_vector)*binwidth)
    
    
    


    ax_hist.scatter(bins, hist, label=r'${n}_{d}$=' + str(nd))

    """
    ax_hist.set_title(r'MOSI='+str(mosi) + ' ,' + r'${\tau}_{0s}=$'+ str(ti_superinfection) +
                      ' ,' + r'${\tau}_{s}=$'+ str(len_superinfection)  + ','+
                      r'${n}_{T}$=' + str(n)
                      + ',' + 'k='+str(k), fontsize=20)
    """
    if burst_histogram ==False:
        ax_hist.set_title(r'MOSI='+str(mosi) + ' ,' + r'${\tau}_{0s}=$'+ str(ti_superinfection) +
                      ' ,' + r'${\tau}_{s}=$'+ str(len_superinfection)  + ','+
                      r'${n}_{T}$=' + str(n)
                      + ',' + 'k='+str(k) + r'$\Delta t = $' + str(delta_t), fontsize=20)
    elif burst_histogram==True:
        ax_hist.set_title(r'MOSI='+str(mosi) + ' ,' + r'${\tau}_{0s}=$'+ str(ti_superinfection) +
                          ' ,' + r'${\tau}_{s}=$'+ str(len_superinfection)  + ','+
                          r'${n}_{T}$=' + str(n)
                          + ',' + 'k='+str(k) + r'$\Delta t = $' + str(delta_t) + ','
                          + 'p=' + str(p) + ',' + 'c=' + str(c), fontsize=20)

    ax_hist.set(xlim=(0.0, 50.0))
    if burst_histogram==True:
        ax_hist.set_xlabel(r'Number of phages', fontsize=18, labelpad=8)
    elif burst_histogram==False:
        ax_hist.set_xlabel(r'${\tau}_{l}$ (min)',fontsize=18, labelpad=8)
    ax_hist.set_ylabel('Normalized Frequency', fontsize=18, labelpad=8)
    ax_hist.tick_params(axis='both', which='major', labelsize=20)
    ax_hist.legend(loc='best', fontsize=18)
    


if save_plot==True:
    if burst_histogram==True:
        fig_name = 'Phage_'+str(n)+'_'+str(k)+'_'+ str(mosi)+'_'+str(ti_superinfection)+'_'+str(len_superinfection)+'_'+str(p)+'.png'
    elif burst_histogram==False:
        fig_name = str(n)+'_'+str(k)+'_'+ str(mosi)+'_'+str(ti_superinfection)+'_'+str(len_superinfection)+'.png'
        
    fig_hist.savefig(os.path.join(os.getcwd(),'parameter_exploring',output_folder,fig_name))
    