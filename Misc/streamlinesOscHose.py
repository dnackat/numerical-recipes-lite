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
y = np.linspace(-5,5,100)

# Streamline function for t = 0 and (x,y) = (0,0)
x1 = (u_0/omega)*(np.cos(omega*y/v_0)-1.)

# Streamline function for t = pi/(2*omega) and (x,y) = (0,0)
x2 = (u_0/omega)*(np.sin(omega*y/v_0))

# Create plot
plt.figure()
#plt.rc('text', usetex=True)
#plt.rc('font', family='serif')
plt.plot(y,x1,'r',linewidth=2,label="r'$t = 0$'")
plt.plot(y,x2,'b',linewidth=2,label="r't = \frac{$\pi$}{$2*\omega$}'")
plt.xlabel("r'$x_s$'")
plt.ylabel("r'$y_s$'")
plt.legend(loc="best")
plt.grid()
plt.show()