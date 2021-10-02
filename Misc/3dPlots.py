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

x = np.linspace(-4,4,100)
y = np.linspace(-4,4,100) 
xx, yy = np.meshgrid(x,y)
zz1 = np.sqrt(xx**2 + yy**2)
zz2 = -np.sqrt(xx**2 + yy**2)

# Plot surface                                             
fig = plt.figure(figsize=(16,16))
plt.rcParams.update({                                                                             
    "text.usetex": True,
    "font.family": "serif",
    "font.sans-serif": ["Helvetica"]})
ax = Axes3D(fig)
#ax.contourf(xx,yy,zz)
ax.plot_surface(xx, yy, zz1, cmap='winter', alpha=0.8)
ax.plot_surface(xx, yy, zz2, cmap='winter', alpha=0.8)
ax.view_init(elev=10, azim=60)
ax.set_xlabel(r'$x$', fontsize=18)
ax.set_ylabel(r'$y$', fontsize=18)
#ax.zaxis.set_rotate_label(False) 
ax.set_zlabel(r'$z$', fontsize=18)
ax.set_title("3D Plot of "+r'$x^2+y^2-z^2=4$', fontsize=20)