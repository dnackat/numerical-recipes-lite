#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 13 01:52:40 2021

@author: dileepn

Plot of a hyperboloid: isobar in an ideal vortex
"""

# Preliminaries
import numpy as np
import matplotlib.pyplot as plt

# Radius of vessel
R = 1.0
r1 = np.linspace(0.1,R,100)
r2 = np.linspace(-R,-0.1,100)

# Free surface (air-water interface height)
H0 = 2.

# v(theta component) = C/r
C = 2.

# Acceleration due to gravity
g = 9.805

# Isobar curve
z1 = H0 + (C**2/(2*g))*(1/R**2 - 1/r1**2)
z2 = H0 + (C**2/(2*g))*(1/R**2 - 1/r2**2)

# Plot
plt.figure(figsize=(12,12))
plt.plot(r1,z1,'b',linewidth=2)
plt.plot(r2,z2,'b',linewidth=2)
plt.grid(axis="both")
plt.show()
