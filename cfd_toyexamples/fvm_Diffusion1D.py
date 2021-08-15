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

# Cell width
dx = fc[1] - fc[0]

# Source terms
S = 6.*cc*dx      # Linear source term, S*d(vol) = 6*x*(dx*1.0)

# Vector of unknowns
phi = np.ones((len(cc),1)) # Length should be cell centroids, 
                              # boundary values lumped with sources
      
# Boundary conditions
phi[0] = 10.
phi[len(phi)-1] = 135.

# Gauss-Seidel method
def gauss(A, b, x, n):

    L = np.tril(A)
    U = A - L
    for i in range(n):
        xprev = x
        x = np.dot(np.linalg.inv(L), b - np.dot(U, x))
        if np.linalg.norm(x - xprev) < 1.e-6:
            break
    return x

