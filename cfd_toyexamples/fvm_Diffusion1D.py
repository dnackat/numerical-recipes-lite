#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Aug 15 20:07:57 2021

@author: dileepn

This is a toy example of using the Finite-Volume method to solve the 1-D 
diffusion equation. 
"""

# Preliminaries
import numpy as np
import matplotlib.pyplot as plt

# Parameters
Gamma = 1.0    # Diffusivity


# Create the grid
fc = np.arange(0.,6.,1)  # Face centroids
cc = (fc + (fc + 1.))/2.      # Cell centroids

# Vector of unknowns
phi = np.ones((len(cc),1))
      

