# -*- coding: utf-8 -*-
"""
Created on Mon Nov 27 11:56:16 2023

@author: Usuario
"""
import numpy as np
import matplotlib.pyplot as plt

def get_collapse_time(df,time,nbacteria):
    initial_bacteria = float(np.array(df['Value'][df['Parameter']=='Initial number of bacteria'])[0])
    collapse_time = time[np.max(np.where(nbacteria > (initial_bacteria*0.0005))[0])]
    return collapse_time

def get_init_collapse_time(df,time,nbacteria):
    initial_bacteria = float(np.array(df['Value'][df['Parameter']=='Initial number of bacteria'])[0])
    collapse_start = time[np.max(np.where(nbacteria > (initial_bacteria*0.8))[0])]
    return collapse_start

def get_collapse_width(df,time,nbacteria):
    collapse_start = get_init_collapse_time(df,time,nbacteria)
    collapse_end = get_collapse_time(df,time,nbacteria)
    return collapse_end - collapse_start

def exists_collapse(time,nbacteria):
    if len(np.where(nbacteria==0.0)[0]) == 0:
        return False
    else:
        return True
    
def generate_colors(num_colors):
    # Create a colormap with the desired number of colors
    cmap = plt.cm.get_cmap('viridis', num_colors)

    # Get the RGB values of the colors in the colormap
    colors_rgb = [cmap(i)[:3] for i in np.linspace(0, 1, num_colors)]

    # Convert RGB values to hex format
    colors_hex = [ "#{:02x}{:02x}{:02x}".format(int(r * 255), int(g * 255), int(b * 255)) for r, g, b in colors_rgb]

    return colors_hex