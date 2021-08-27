#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 26 00:47:31 2021

@author: dileepn

An example of Lagrange multiplier method for constrained optimization

Finding the maximum perimeter rectangle that can be inscribed in an ellipse
given by: x^2 + 4y^2 = 4. 
"""

# Preliminaries
import numpy as np
import matplotlib.pyplot as plt

# Grid
x = np.arange(-7.07, 7.07, 0.125)
y = np.arange(-7.07, 7.07, 0.125)

xx, yy = np.meshgrid(x, y)

# Objective function: perimeter of the inscribed rectangle
f = xx**2 - 6*xx + yy**2 - 6*yy + xx*yy 

# Constraint curve: a circle with radius 5 centered around the origin
g = xx*xx + yy*yy

# Coloring
color = g

# Gradients
fx, fy = np.gradient(f, 0.5, 0.5)
gx, gy = np.gradient(g, 0.5, 0.5)

# Plot

# Gradeint and constraint
cx = np.linspace(-7.07,7.07,100)
plt.figure(figsize=(16,16)) 
plt.plot(cx,np.sqrt(50-cx**2),'r',linewidth=2.5)
plt.plot(cx,-np.sqrt(50-cx**2),'r',linewidth=2.5)
#plt.quiver(xx, yy, fx, fy, color)
#plt.quiver(xx, yy, gx, gy)
plt.xlabel(r'$x$',fontsize=20)
plt.ylabel(r'$y$',fontsize=20)

# Level curves
level_curves = [-10, -6, 0, 10, 30, 40, 50, 60, 75,100]
#c = plt.pcolormesh(xx, yy, f, cmap='nipy_spectral', vmin=f.min(), vmax=f.max(), shading='auto')
plt.axis([xx.min(), xx.max(), yy.min(), yy.max()])
CS = plt.contour(xx, yy, f, level_curves, colors='black', linewidths=0.75)
plt.clabel(CS, inline=True, fmt='%1.0f', fontsize=16)
#plt.colorbar(c)
plt.show()