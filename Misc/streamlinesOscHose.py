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
omega = 2 # oscillations per second


# Create a 1D grid
y = np.linspace(-5,5,100)

# Streamline function for t = 0 and (x,y) = (0,0)
x = 