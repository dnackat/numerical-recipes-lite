#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 30 20:43:21 2021

@author: dileepn

Simplified planar case of Hagen-Poiseuille flow of oil in a lubricated pipe. 
"""

# Preliminaries
import numpy as np
import matplotlib.pyplot as plt

#### Non-dimensionalized velocity profile in oil ####
y_by_H = np.linspace(0.,1.,100) # Non-dimensional vertical distance
H1_by_H = 0.8 # Ratio of oil height to pipe height

# Viscosity ratio of oil to water in the pipe
visc_ratio0 = 1
visc_ratio1 = 2
visc_ratio2 = 3
visc_ratio3 = 10
visc_ratio4 = 100

V0 = (1 - H1_by_H**2) + visc_ratio0*(H1_by_H**2 - y_by_H**2)
V1 = (1 - H1_by_H**2) + visc_ratio1*(H1_by_H**2 - y_by_H**2)
V2 = (1 - H1_by_H**2) + visc_ratio2*(H1_by_H**2 - y_by_H**2)
V3 = (1 - H1_by_H**2) + visc_ratio3*(H1_by_H**2 - y_by_H**2)
V4 = (1 - H1_by_H**2) + visc_ratio4*(H1_by_H**2 - y_by_H**2)

## Plot the velocity profile ##
plt.figure(figsize=(16,16))
plt.rcParams.update({
    "text.usetex": True,
    "font.family": "sans-serif",
    "font.sans-serif": ["Helvetica"]})
plt.plot(y_by_H,V0,'b',label=r'$\frac{\mu_{o}}{\mu_{w}} = 1$',linewidth=2)
plt.plot(y_by_H,V1,'g',label=r'$\frac{\mu_{o}}{\mu_{w}} = 2$',linewidth=2)
plt.plot(y_by_H,V2,'o',label=r'$\frac{\mu_{o}}{\mu_{w}} = 3$',linewidth=2)
plt.plot(y_by_H,V3,'r',label=r'$\frac{\mu_{o}}{\mu_{w}} = 10$',linewidth=2)
plt.plot(y_by_H,V4,'b',label=r'$\frac{\mu_{o}}{\mu_{w}} = 100$',linewidth=2)
plt.xlabel(r'$\frac{V_{x}(y)}{(-\frac{dP}{dx})(\frac{H^2}{2\mu_{w}})}$',fontsize=16)
plt.ylabel(r'$\frac{y}{H}$',fontsize=16)
plt.yticks(


