#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Oct  2 20:05:19 2021

@author: dileepn

3D plots: x^2 + y^2 - z^2 = 4
"""

# Preliminaries
import numpy as np
import matplotlib.pyplot as plt
#from mpl_toolkits.mplot3d import Axes3D

u = np.linspace(-2,2,200)
v = np.linspace(0,2*np.pi,200)
uu, vv = np.meshgrid(u,v) 

a = 2
b = 2
c = 2

x = a*np.cosh(uu)*np.cos(vv)
y = b*np.cosh(uu)*np.sin(vv)
z = c*np.sinh(uu)

# Plot surface                                             
fig = plt.figure(figsize=(16,16))
plt.rcParams.update({                                                                             
    "text.usetex": True,
    "font.family": "serif",
    "font.sans-serif": ["Helvetica"]})

ax = fig.add_subplot(111, projection='3d')
#ax = Axes3D(fig)
ax.plot_surface(x, y, z, cmap='winter', alpha=0.8)
ax.view_init(elev=10, azim=60)
ax.set_xlabel(r'$x$', fontsize=18)
ax.set_ylabel(r'$y$', fontsize=18)
#ax.zaxis.set_rotate_label(False) 
ax.set_zlabel(r'$z$', fontsize=18)
ax.set_title("3D Plot of "+r'$x^2+y^2-z^2=4$', fontsize=20)