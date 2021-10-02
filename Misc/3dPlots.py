#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Oct  2 20:05:19 2021

@author: dileepn

3D plots
"""

# Preliminaries
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

x = np.linspace(-2,2,100)
y1 = np.sqrt(4 - x**2)
y2 = -np.sqrt(4 - x**2)
xx, yy1 = np.meshgrid(x,y1)
xx, yy2 = np.meshgrid(x,y2)
zz1 = np.sqrt(xx**2 + yy1**2 - 4)
zz2 = -np.sqrt(xx**2 + yy2**2 - 4)

# Plot surface                                             
fig = plt.figure(figsize=(16,16))
plt.rcParams.update({                                                                             
    "text.usetex": True,
    "font.family": "serif",
    "font.sans-serif": ["Helvetica"]})
ax = Axes3D(fig)
#ax.contourf(xx,yy,zz)
ax.plot_surface(xx, yy1, zz1, cmap='winter', alpha=0.8)
ax.plot_surface(xx, yy2, zz2, cmap='winter', alpha=0.8)
ax.view_init(elev=10, azim=60)
ax.set_xlabel(r'$x$', fontsize=18)
ax.set_ylabel(r'$y$', fontsize=18)
#ax.zaxis.set_rotate_label(False) 
ax.set_zlabel(r'$z$', fontsize=18)
ax.set_title("3D Plot of "+r'$x^2+y^2-z^2=4$', fontsize=20)