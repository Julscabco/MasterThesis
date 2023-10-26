# -*- coding: utf-8 -*-
"""
Created on Thu Oct  5 15:44:49 2023

@author: Usuario
"""
import os 

import matplotlib.pyplot as plt
import numpy as np

def get_burst_size(t, c=0.5, p=5.0):
    return c*np.exp(t*p)

def get_sigmoidal_burst_size(t, c=0.5, p=0.5, t0 = 30.0, a=1.0):
    tlim = 38.0
    if t<tlim:
        return float(c)/(1.0+np.exp(-(t-t0)*p))
    elif t>=tlim:
        return float(c)/(1.0+np.exp(-(tlim-t0)*p))+ a*(t-tlim)


file_name = '0.01_10000_60_6_33_7.6_15_3_10'
#folder = 'CollectiveGillespie_data'
#folder = 'CollectiveConstantInterval_data'
folder='CollectiveConstantInterval2_data'
file_path = os.path.join(os.getcwd(),os.path.join(folder,file_name))
fig_hist, ax_hist = plt.subplots(figsize=(10, 6))
save_plot=True

# If the histogram needs to be in units of phages
burst_histogram=True
p = 0.4
c = 20.0
t0 = 27.0
a = 0.1


f = open(file_path,"r")
t = []
for line in f.readlines():
    t.append(float(line))

hist, bins = np.histogram(t, bins=50)
bins = (bins[1:] + bins[:-1])/2.0
binwidth = bins[2]-bins[1]

if burst_histogram==True:
    hist = [hist[i]*get_sigmoidal_burst_size(bins[i],c,p,t0,a) for i in range(0,len(hist))]
    hist = hist/(len(t)*binwidth*np.max(hist))
else:
    hist = hist/(len(t)*binwidth)



ax_hist.scatter(bins, hist)

#ax_hist.set(xlim=(0.0, 50.0))
ax_hist.set_xlabel(r'${\tau}_{l} (min)$', fontsize=18, labelpad=8)
ax_hist.set_ylabel('Normalized Frequency', fontsize=18, labelpad=8)
ax_hist.tick_params(axis='both', which='major', labelsize=20)
#ax_hist.legend(loc='best', fontsize=18)

#Nb,n,k,nd,mosi,ti,length,burst_size = file_name.split('_')
#delta_t,Nb,n,k,nd,mosi,ti,length = file_name.split('_')
delta_t,Nb,n,k,nd,mosi,ti,length,burst_size = file_name.split('_')

if burst_histogram==False:
    ax_hist.set_title(r'${N}_{b}$='+str(Nb) +','+ r'MOSI='+str(mosi) + ' ,' + r'${\tau}_{0s}=$'+ str(ti) +
                  ' ,' + r'${\tau}_{s}=$'+ str(length)  + ','+
                  r'${n}_{T}$=' + str(n)
                  + ',' + 'k='+str(k) + ',' + r'${n}_{d}=$'+str(nd) +
                  ','+ r'$\Delta {N}_{p}=$'+str(burst_size), fontsize=20)
elif burst_histogram==True:
    ax_hist.set_title(r'MOSI='+str(mosi) + ',' + 'k='+str(k) + ',' + r'${n}_{d}=$'+str(nd) +
                  ','+ r'$\Delta {N}_{p}=$'+str(burst_size) + ','
                  'p='+str(p)+','+'c='+str(c) +','+ r'${t}_{0}=$' + str(t0)
                  +',' +'a='+str(a), fontsize=20)

"""
ax_hist.set_title(r'${N}_{b}$='+str(Nb) +','+ r'MOSI='+str(mosi) + ' ,' + r'${\tau}_{0s}=$'+ str(ti) +
                  ' ,' + r'${\tau}_{s}=$'+ str(length)  + ','+
                  r'${n}_{T}$=' + str(n)
                  + ',' + 'k='+str(k) + ',' + r'${n}_{d}=$'+str(nd) +
                 ',' + r'$\Delta t$=' + str(delta)
"""

output_folder = 'CollectiveConstantTimeInterval2'
#output_folder = 'CollectiveGillespie'
if save_plot==True:
    
    if burst_histogram==True:
        #figure_name = 'Phage_sigmoidal'+file_name + '.png'
        figure_name = 'Phage_sigmoidal_linear'+file_name+'_'+str(p)+'_'+str(c)+'_'+str(t0)+'_'+str(a)+ '.png'
        
    elif burst_histogram==False:
        figure_name = file_name+'.png'
        
    fig_hist.savefig(os.path.join(os.getcwd(),'parameter_exploring',output_folder,figure_name))
        
