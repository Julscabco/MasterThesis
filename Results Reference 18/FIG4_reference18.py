# -*- coding: utf-8 -*-
"""
Created on Mon Oct  2 10:09:56 2023

@author: Usuario
"""

import os

import numpy as np
import matplotlib.pyplot as plt

# The plot
fig_hist, ax_hist = plt.subplots(figsize=(10, 6))
file_directory = os.path.join(os.getcwd(),os.path.join('CollectiveConstantTimeInterval2_data_noburst_fig4'))
save_plot=True


for filename in os.listdir(file_directory):
    #n,nd,k,mosi,ti,length,ndata = filename.split('txt')[0].split('_')
    
    # For CollectiveConstantTimeInterval2
    deltat,Nb,n,k,nd,mosi,ti,length,burst_size= filename.split('txt')[0].split('_')
    
    f = open(os.path.join(file_directory,filename),"r")
    tau_vector=[]
    for line in f.readlines():
        tau_vector.append(float(line.split('\n')[0]))
        
        
    
    # Calculation and normamlization of the histogram
    hist, bins = np.histogram(tau_vector, bins=60)
    bins = (bins[1:] + bins[:-1])/2.0
    binwidth = bins[2]-bins[1]
    hist = hist/(len(tau_vector)*binwidth)

    ax_hist.scatter(bins, hist, label=r'${\tau}_{0s}$=' + str(ti) +','+ r'${\tau}_{s}=$' + str(length))

    """
    ax_hist.set_title(r'MOSI='+str(mosi) +  ',' +
                      'ndata='+str(20000) + ',' +  
                      r'${n}_{T}$=' + str(n) + ','+
                      r'${n}_{d}=$'+str(nd) + ',' +
                      'k='+str(k), fontsize=20)
    """
    
    # For collective constant time interval 2
    ax_hist.set_title(r'${N}_{b}$='+str(Nb) +','+ r'MOSI='+str(mosi)  + ','+
                  r'${n}_{T}$=' + str(n)
                  + ',' + 'k='+str(k) + ',' + r'${n}_{d}=$'+str(nd) +
                  ','+ r'$\Delta t=$'+str(deltat) +',' + r'$\Delta {N}_{p}=$'+str(burst_size), fontsize=20)    
    

    ax_hist.set(xlim=(30.0, 55.0))
    ax_hist.set_xlabel(r'${\tau}_{l} (min)$', fontsize=18, labelpad=8)
    ax_hist.set_ylabel('Normalized Frequency', fontsize=18, labelpad=8)
    ax_hist.tick_params(axis='both', which='major', labelsize=20)
    ax_hist.legend(loc='best', fontsize=18)
    


if save_plot==True:
    #fig_name = str(n)+'_'+str(nd)+'_'+str(k)+'_'+ str(mosi)+'_'+'.png'
    fig_name = str(n)+'_'+str(nd)+'_'+str(k)+'_'+ str(mosi)+'_'+str(burst_size)+'.png'
    fig_hist.savefig(os.path.join(os.getcwd(),'parameter_exploring','FIG4_CollectiveTimeInterval',fig_name))
    