#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Aug 15 20:07:57 2021

@author: dileepn

This is a toy example of using the Finite-Volume method to solve the 1-D 
diffusion equation of the form:
     d/dx(Gamma*d(phi)/dx) + source term = 0
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
cc = (fc[:-1] + (fc[:-1] + 1.))/2.      # Cell centroids

# Source terms
S = 6.*cc*dx      # Linear source term, S*d(vol) = 6*x_centroid*(dx*1.0)

# Vector of unknowns
phi = np.ones((len(cc),1)) # Length should be equal to no. of cell centroids; 
                           # boundary values lumped with constants in b vector
                              
# Boundary conditions for phi (included with constants in b vector)
phi_first = 10.     # phi on first boundary face
phi_last = 135.     # phi on last boundary face

# Vector of constants
b = np.zeros((len(cc),1)) # Length = no. of cell centroids

# Populate the vector of constants
b = S     # Source term contribution
b[0] += (Gamma/(cc[0]-fc[0]))*phi_first    # Add the boundary value for first face
b[-1] += (Gamma/(fc[-1]-cc[-1]))*phi_last  #  Add the boundary value for last face
      
# Matrix of coefficients (size = size(b)*size(phi))
A = np.zeros((len(b),len(phi)))

# Populate A. It is sparse and diagonally dominant, so a single loop should suffice
for i in range(len(b)):
     # First and last rows correspond to boundary cells and have just one element  
     if i == 0: # First boundary
          A[i,i] = Gamma*(1./(cc[i]-fc[i]) + 1./(cc[i+1]-cc[i]))
     elif i == len(b)-1: # Last boundary
          A[i,i] = Gamma*(1./(fc[i]-cc[i]) + 1./(cc[i]-cc[i-1]))
     else: # Interior cells
          A[i,i] = Gamma/()
     

# Gauss-Seidel iterative method
def gauss(A, b, x, n):
     """ This is a function that uses the Gauss-Seidel iterative scheme from 
         numerical linear algebra to solve a linear system of the form 
         Ax = b. The scheme updates x using the following algorithm:
         x(k+1) = inv(L)*(b - U*x(k)), where U and L are upper and lower
         triangular matrices obtained by performing LU decmposition of A."""

     L = np.tril(A)
     U = A - L
     for i in range(n):
        xprev = x
        x = np.dot(np.linalg.inv(L), b - np.dot(U, x))
        if np.linalg.norm(x - xprev) < 1.e-6:
            break
     return x

