#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug  9 14:53:35 2021

@author: dileepn

Temperature profile of a geothermal spring
"""
import matplotlib.pyplot as plt
import numpy as np

def temp(x,y):
     # Temp at (x,y)
     return 450./np.sqrt(x**2+y**2+1)+420./np.sqrt((x+10.)**2+(y-20.)**2+1)
     
def plot_2d():
    # Make data.
    X = np.arange(-25, 35, 0.125)
    Y = np.arange(-25, 35, 0.125)
    X, Y = np.meshgrid(X, Y)
    Z = temp(X, Y)
    level_curves = [20, 30, 40, 50, 60, 70, 80, 120, 200]

    # Plot the surface.
    #plt.figure(size=(12,12))
    fig, ax = plt.subplots(figsize=(14,12))
    c = ax.pcolormesh(X, Y, Z, cmap='nipy_spectral', vmin=Z.min(), vmax=Z.max(), shading='auto')
    ax.set_title('Hot Spring')
    # set the limits of the plot to the limits of the data
    ax.axis([X.min(), X.max(), Y.min(), Y.max()])
    fig.colorbar(c, ax=ax)
    CS = ax.contour(X, Y, Z, level_curves, colors='white', linewidths=0.75)
    ax.clabel(CS, inline=True, fmt='%1.0f', fontsize=16)
    plt.show()
    
    
plot_2d()