#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 28 19:49:45 2020

@author: dileepn

Vector field plots
"""

import numpy as np
import matplotlib.pyplot as plt

#%% Vector plots
x,y = np.meshgrid(np.linspace(-2,2,10),np.linspace(-2,2,10))

u = -y/np.sqrt(x**2 + y**2)
v = x/np.sqrt(x**2 + y**2)

#z = 1/(2*np.pi)*np.exp(-(x**2 + y**2)/2)

plt.quiver(x,y,u,v,color=['red','green','blue'])
plt.grid(axis='both')

#%% Harmonic functions
x,y = np.meshgrid(np.linspace(-5,5,100),np.linspace(-5,5,100))

a = -1.25
b = 2.75
k = 2.5

z = a*(x**2 - y**2) + b*x*y
z1 = np.exp(k*x)*np.sin(k*y)  

plt.contourf(x,y,z,cmap='winter')

#%% Plane curve where every point bisects the part of the tangent in the first quadrant passing through that point
c = 2

x = np.linspace(0.1,5,100)

# Curve that satisfies the condition
y = c/x

# Point for checking lengths
x1 = 1.25
y1 = c/x1

# Points of intersection of tnagnet line with x and y-axes
a = 2*x1
b = 2*y1

# Calculate lengths
len1 = np.sqrt(x1**2 + (b-y1)**2)
len2 = np.sqrt((a-x1)**2 + y1**2)

if abs(len1-len2) < 1e-3:
     print("Close enough. P is a bisector of the tangent line.")
else:
     print("Nope. P is not the bisector of the tangent line.")


#%% Item response curve

a = 5
b = 0
c = 0.33

theta = np.linspace(-5,5,100)

P = c + (1-c)/(1+np.exp(-a*(theta-b)))

