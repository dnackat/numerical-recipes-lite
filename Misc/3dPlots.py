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

xx, yy = np.meshgrid(np.linspace(-4,4,100),np.linspace(-4,4,100))
zz = np.sqrt(xx**2 + yy**2 - 4.)

# Plot surface
fig = plt.figure(figsize=(16,16))
plt.rcParams.update({
    "text.usetex": True,
    "font.family": "serif",
    "font.sans-serif": ["Helvetica"]})
ax = Axes3D(fig)
ax.plot_surface(xx, yy, zz, cmap='winter', alpha=0.2)
ax.view_init(elev=10, azim=60)
ax.set_xlabel(r'$x$', fontsize=12)
ax.set_ylabel(r'$y$', fontsize=12)
#ax.zaxis.set_rotate_label(False) 
ax.set_zlabel(r'$\Phi_3 = x_2^2$', fontsize=12)
ax.set_title("3D Plot of "+r'$x^2+y^2-z^2=4$', fontsize=20)