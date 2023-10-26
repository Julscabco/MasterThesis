# -*- coding: utf-8 -*-
"""
Created on Wed Oct 11 10:19:36 2023

@author: Usuario
"""

import os

import numpy as np
import matplotlib.pyplot as plt

""" This code generates a graph to compare if the results with a constant time 
interval and just one bacteria coincide with the ones for the case with a lot of 
bacteria """
# Save plot?
save_plot = False

# file names
filename_individual = '0.01_60_33_6.0_1.6_15.0_3.0_10000.txt'
filename_several = '0.01_10000_60_6_33_7.6_15_3_10'

# Complete path
folder_individual_bacteria = os.path.join(os.getcwd(),'ConstantTimeInterval','0.01_60_6_1.6_15_3_10000',filename_individual)
folder_several_bacteria = os.path.join(os.getcwd(),'CollectiveConstantInterval2_data',filename_several)

# We extract parameters
deltat,Nb,n,k,nd,mosi,ti,length,burst = filename_several.split('_')
burst = burst.split('.txt')[0]


# We read the data in the files
f1 = open(folder_individual_bacteria,"r")
tau_individual = []
for line in f1.readlines():
    tau_individual.append(float(line.split('\n')[0]))
f1.close()


f2 = open(folder_several_bacteria,"r")
tau_several = []
for line in f2.readlines():
    tau_several.append(float(line.split('\n')[0]))
f2.close()



# We create the figure
fig_hist, ax_hist = plt.subplots(figsize=(10, 6))


# Calculation and normamlization of the histogram
nbins = 60

hist1, bins1 = np.histogram(tau_individual, bins=nbins)
hist2, bins2 = np.histogram(tau_several, bins=nbins)

bins1 = (bins1[1:] + bins1[:-1])/2.0
binwidth1 = bins1[2]-bins1[1]
hist1 = hist1/(len(tau_individual)*binwidth1)

bins2 = (bins2[1:] + bins2[:-1])/2.0
binwidth2 = bins2[2]-bins2[1]
hist2 = hist2/(len(tau_several)*binwidth2)

ax_hist.scatter(bins1, hist1, label=r'${N}_{b}=$'+str(1))
ax_hist.scatter(bins2, hist2, label=r'${N}_{b}=$'+str(Nb) + ', '+ r'${\Delta N}_{p}=$'+str(burst))

ax_hist.set_title(r'MOSI='+str(mosi) + ' ,' + r'${\tau}_{0s}=$'+ str(ti) +
                   ' ,' + r'${\tau}_{s}=$'+ str(length)  + ','+
                   r'${n}_{T}$=' + str(n) + r'${n}_{d}=$' + str(nd)
                   + ',' + 'k='+str(k) + ',' + r'$\Delta t$='+ str(deltat), fontsize=20)

ax_hist.set(xlim=(0.0, 70.0))
ax_hist.set_xlabel(r'${\tau}_{l} (min)$', fontsize=18, labelpad=8)
ax_hist.set_ylabel('Normalized Frequency', fontsize=18, labelpad=8)
ax_hist.tick_params(axis='both', which='major', labelsize=20)
ax_hist.legend(loc='best', fontsize=18)

output_folder = 'Comparison_CTI_ind_sev'
if save_plot==True:
    fig_hist.savefig(os.path.join(os.getcwd(),'parameter_exploring',output_folder,filename_several+'.png'))



