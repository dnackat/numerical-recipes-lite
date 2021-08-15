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
numcells = 5   # Number of cells in the grid
dx = 1.   # Cell width
fc = np.arange(0.,numcells+1,dx)  # Face centroids
cc = (fc + (fc + 1.))/2.      # Cell centroids

# Source terms
S = 6.*cc*dx      # Linear source term, S*d(vol) = 6*x_centroid*(dx*1.0)

# Vector of unknowns
phi = np.ones((len(cc),1)) # Length should be equal to no. of cell centroids; 
                              # boundary values lumped with constants in b
                              
# Boundary conditions for phi (included with constants in b vector)


# Vector of constants
b = np.zeros((len(cc),1)) # Length = no. of cell centroids

# Populate the vector of constants
b = S 
b[0] += Gamma/(cc[0]-fc[0])   # Add the boundary value for first face
b[-1] += Gamma/(fc[-1]-cc[-1])     #  Add the boundary value for last face
      
# Matrix of coefficients (size = size(b)*size(phi))
A = np.zeros((len(b),len(phi)))

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

