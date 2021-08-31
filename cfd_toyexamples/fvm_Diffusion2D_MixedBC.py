#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 31 21:03:25 2021

@author: dileepn

This is a toy example of using the Finite-Volume method to solve the 2-D 
heat diffusion problem with mixed boundary conditons (Dirichlet, von Neumann,
and a combination of the two). 
The physics is governed by the equation:
     div(rho*u*T) = div(k*grad(T)) + source term
"""

# Preliminaries
import numpy as np
import matplotlib.pyplot as plt

#%%
# Parameters
k = 1.0    # Constant thermal conductivity in W/m-K
h = 10.0   # Constant convective heat transfer coefficient in W/m2-K
domain_length = 1.0 # In both x and y directions in m

# Create the grid
numcells = 5   # Number of cells in the x and y directions
fc_x = np.linspace(0.,domain_length,numcells+1)  # x face centroids
fc_y = np.linspace(0.,domain_length,numcells+1)  # y face centroids
fxx = np.tile(fc_x,numcells+1)       # Long vector of x face centroids
fyy = np.repeat(fc_y,numcells+1)       # Long vector of y face centroids
cc_x = (fc_x[:-1] + (fc_x[1:]))/2.      # x cell centroids
cc_y = (fc_y[:-1] + (fc_y[1:]))/2.      # y cell centroids
cxx = np.tile(cc_x,numcells)       # Long vector of x cell centroids
cyy = np.repeat(cc_y,numcells)       # Long vector of y cell centroids

# Calculate all the cell distances for convenience as they remain constant in this problem
delta_x = fc_x[2] - fc_x[1] # Distance between faces in the x direction
delta_y = fc_y[2] - fc_y[1] # Distance between faces in the y direction
del_x = cc_x[2] - cc_x[1]   # Distance between cell centroids in the x direction 
del_x_b = cc_x[0] - fc_x[0] # Distance between cell centroid and boundary face in x (right, left)
del_y = cc_y[2] - cc_y[1]   # Distance between cell centroids in the y direction 
del_y_b = cc_y[0] - fc_y[0] # Distance between cell centroid and boundary face in y (top, bottom)

# Cell areas and volumes
thickness = 1.0     # 2D case
area_x = delta_y*thickness # Areas of faces with normal vector aligned with the x-axis
area_y = delta_x*thickness # Areas of faces with normal vector aligned with the y-axis
cell_volume = delta_x*delta_y*thickness

#%%
# Volumetric source. Depends on cell centroid values
S = 100.*(1. + 5*cxx + 5*cyy)*cell_volume      # Linear source term, S*d(vol) = 6*x_centroid*(dx*1.0)

 #%%
# Vector of unknown temperatures
T = np.ones((len(cc_x),1)) # Length should be equal to no. of cell centroids; 
                           # boundary values lumped with constants in b vector
                              
# Boundary conditions for T (included with constants in b vector)
T_right = 500.    # T in K on right boundary face
T_left = 500.     # T in K on left boundary face
T_top = 500.      # T in K on top boundary face
T_inf = 300.      # T in K of flow across bottom face

# Vector of constants.
b = np.zeros((len(cc_x),1)) # Length = no. of cell centroids

# Populate the vector of constants
b = S.reshape((len(S),1))     # Source term contribution
b[0] += (k/(cc[0]-fc[0]))*T_first    # Add the boundary value for first face
b[-1] += (k/(fc[-1]-cc[-1]))*T_last  #  Add the boundary value for last face
      
# Matrix of coefficients (size = size(b)*size(T))
A = np.zeros((len(b),len(T)))

# Populate A. It is sparse and diagonally dominant, so a single loop should suffice
for i in range(len(b)): 
     if i == 0: # First boundary
          A[i,i] = k*(1./(cc[i]-fc[i]) + 1./(cc[i+1]-cc[i]))
          A[i,i+1] = -k/(cc[i+1]-cc[i]) # Right of diagonal
     elif i == len(b)-1: # Last boundary
          A[i,i] = k*(1./(fc[i+1]-cc[i]) + 1./(cc[i]-cc[i-1]))
          A[i,i-1] = -k/(cc[i]-cc[i-1]) # Left of diagonal
     else: # Interior cells
          A[i,i] = k*(1./(cc[i+1]-cc[i]) + 1./(cc[i]-cc[i-1])) # Diagonal element
          A[i,i-1] = -k/(cc[i]-cc[i-1]) # Left of diagonal
          A[i,i+1] = -k/(cc[i+1]-cc[i]) # Right of diagonal
     

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
T_fvm = gauss(A, b, T, 100)

# Analytical (exact) solution
x = np.linspace(0.,5.,100)
T_exact = 10 + 50*x - x**3

# Print the result
print("\n")
print("---------- Solution ----------")
print("\n")
print("T values from FVM method:\n")
print(T_fvm.round(2))

# Plot the FVM and exact solution
cc_plot = np.zeros(len(cc)+2)   # Add boundary points just for plotting
cc_plot[-1] = fc[len(cc)]
cc_plot[1:-1] = cc
T_plot = np.zeros(len(cc)+2)  # Add boundary points just for plotting
T_plot[0] = T_first
T_plot[-1] = T_last
T_plot[1:-1] = T_fvm.reshape(len(T_fvm))

plt.figure(figsize=(16,16))
plt.plot(x,T_exact,'k',linewidth=1,label="Exact solution") # Plot of the exact solution
plt.scatter(cc_plot,T_plot,c="red",marker="o",s=200,label="FVM soluton") # Scatter plot of FVM resluts
for c in cc:
     plt.axvline(x=c,color='b',linestyle='--',linewidth=0.5)    # Vertical lines from centroid for clarity
plt.grid(axis="both")
plt.xlim((0.,5.))
plt.xlabel(r'$x$',fontsize=16)
plt.ylabel(r'$\T$',fontsize=16)
plt.legend(loc="best",fontsize=16)
plt.title("Comparison between FVM and exact solution for the 1D Diffusion Equation"\
          ,fontsize=20) 
