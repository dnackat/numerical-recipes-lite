#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 23 23:13:54 2021

@author: dileepn

Streamlines passing through origin of an oscillating garden hose at t = 0 and
t = pi/(2*omega), where omega is the speed of oscillation 
"""

# Preliminaries
import numpy as np 
import matplotlib.pyplot as plt

# Constants
u_0 = 1. # m/s
v_0 = 2. # m/s
omega = 2. # oscillations per second


# Create a 1D grid
y = np.linspace(0,10,100)

# Streamline function for t = 0 and (x,y) = (0,0)
x1 = (u_0/omega)*(np.cos(omega*y/v_0)-1.)

# Streamline function for t = pi/(2*omega) and (x,y) = (0,0)
x2 = (u_0/omega)*(np.sin(omega*y/v_0))

# Create plot
plt.figure(figsize=(16,16))
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.plot(x1,y,'r',linewidth=2,label=r'$t = 0$')
plt.plot(x2,y,'b',linewidth=2,label=r'$t = \frac{\pi}{2\omega}$')
plt.xlabel(r'$x_s$',fontsize=20)
plt.ylabel(r'$y_s$',fontsize=20)
plt.legend(fontsize=20)
plt.grid()
plt.show()