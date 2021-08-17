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
dx = 1   # Cell width
fc = np.arange(0.,numcells+dx,dx)  # Face centroids
cc = (fc[:-1] + (fc[:-1] + dx))/2.      # Cell centroids

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
b = S.reshape((len(S),1))     # Source term contribution
b[0] += (Gamma/(cc[0]-fc[0]))*phi_first    # Add the boundary value for first face
b[-1] += (Gamma/(fc[-1]-cc[-1]))*phi_last  #  Add the boundary value for last face
      
# Matrix of coefficients (size = size(b)*size(phi))
A = np.zeros((len(b),len(phi)))

# Populate A. It is sparse and diagonally dominant, so a single loop should suffice
for i in range(len(b)): 
     if i == 0: # First boundary
          A[i,i] = Gamma*(1./(cc[i]-fc[i]) + 1./(cc[i+1]-cc[i]))
          A[i,i+1] = -Gamma/(cc[i+1]-cc[i]) # Right of diagonal
     elif i == len(b)-1: # Last boundary
          A[i,i] = Gamma*(1./(fc[i+1]-cc[i]) + 1./(cc[i]-cc[i-1]))
          A[i,i-1] = -Gamma/(cc[i]-cc[i-1]) # Left of diagonal
     else: # Interior cells
          A[i,i] = Gamma*(1./(cc[i+1]-cc[i]) + 1./(cc[i]-cc[i-1])) # Diagonal element
          A[i,i-1] = -Gamma/(cc[i]-cc[i-1]) # Left of diagonal
          A[i,i+1] = -Gamma/(cc[i+1]-cc[i]) # Right of diagonal
     

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
            print("\n")
            print("Converged in",i+1,"iterations.")
            break
     return x

# Solve the system 
phi_fvm = gauss(A, b, phi, 100)

# Analytical (exact) solution
x = np.linspace(0.,5.,100)
phi_exact = 10 + 50*x - x**3

# Print the result
print("\n")
print("---------- Solution ----------")
print("\n")
print("PHI values from FVM method:\n")
print(phi_fvm.round(2))

# Plot the FVM and exact solution
cc_plot = np.zeros(len(cc)+2)   # Add boundary points just for plotting
cc_plot[-1] = fc[len(cc)]
cc_plot[1:-1] = cc
phi_plot = np.zeros(len(cc)+2)  # Add boundary points just for plotting
phi_plot[0] = phi_first
phi_plot[-1] = phi_last
phi_plot[1:-1] = phi_fvm.reshape(len(phi_fvm))

plt.figure(figsize=(16,16))
plt.plot(x,phi_exact,'k',linewidth=1,label="Exact solution") # Plot of the exact solution
plt.scatter(cc_plot,phi_plot,c="red",marker="o",s=200,label="FVM soluton") # Scatter plot of FVM resluts
for c in cc:
     plt.axvline(x=c,color='b',linestyle='--',linewidth=0.5)    # Vertical lines from centroid for clarity
plt.grid(axis="both")
plt.xlim((0.,5.))
plt.xlabel(r'$x$',fontsize=16)
plt.ylabel(r'$\phi$',fontsize=16)
plt.legend(loc="best",fontsize=16)
plt.title("Comparison between FVM and exact solution for the 1D Diffusion Equation"\
          ,fontsize=20) 
