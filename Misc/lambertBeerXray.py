#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 14 19:43:21 2021

@author: dileepn

Lambert-Beer model for X-ray attenuation
"""
# Preliminaries
import numpy as np
import matplotlib.pyplot as plt

# Intial intensity 
I_0 = 1.0

# Attenuation coefficient
mu = 1.0

# length of passage
dx = 1.0


# Combinations to try
data = [(0.5,1),(0.5,2),(1,2),(2,0.5),(3,3)]

for entry in data:
     # Unpack tuple
     mu, dx = entry
     
     # Intensity calculation
     I = I_0*np.exp(-mu*dx)
     
     # Print result     
     print("Intensity for mu = {:.1f} and dx = {} is {:.2E}".format(mu,dx,I))
     