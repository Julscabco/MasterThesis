# -*- coding: utf-8 -*-
"""
Created on Wed Oct 11 13:24:12 2023

@author: Usuario
"""

import os

import numpy as np
import matplotlib.pyplot as plt

def get_burst_size(t, p=1.0):
    return int(t*p)

# Path to the file
folder_path = os.path.join(os.getcwd(),'ConstantTimeInterval','Phage_histogram')
filename = ' '

# We get the parameters
delta_t,n,nd,k,mosi,ti,length,ndata = filename.split('txt')[0].split('_')

# We open the file and read the data
f = open(os.path.join(folder_path,filename),"r")
tau_vector=[]

for line in f.readlines():
    tau_vector.append(float(line.split('\n')[0]))