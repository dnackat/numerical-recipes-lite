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

u = x/np.sqrt(x**2 + y**2)
v = y/np.sqrt(x**2 + y**2)

#z = 1/(2*np.pi)*np.exp(-(x**2 + y**2)/2)

plt.quiver(x,y,u,v,color='rgb')
plt.grid(axis='both')
