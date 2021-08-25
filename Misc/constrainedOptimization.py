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
x = np.arange(0, 2.5, 0.25)
y = np.arange(0, 2.5, 0.25)

xx, yy = np.meshgrid(x, y)

# Objective function: perimeter of the inscribed rectangle
f = 4*xx + 4*yy 

# Constraint curve: a circle with radius 2 centered around the origin
g = xx*xx + 4*yy*yy

# Coloring
color = g

# Gradients
fx, fy = np.gradient(f, 0.25, 0.25)
gx, gy = np.gradient(g, 0.25, 0.25)

# Plot
cx = np.linspace(0,2,100)
plt.figure(figsize=(16,16)) 
plt.plot(cx,np.sqrt(1-cx**2/4.),linewidth=2.5)
plt.quiver(xx, yy, fx, fy, color)
plt.quiver(xx, yy, gx, gy)
plt.xlabel(r'$x$',fontsize=20)
plt.ylabel(r'$y$',fontsize=20)