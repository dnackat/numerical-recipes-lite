#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct  6 22:25:40 2021

@author: dileepn

Solution 2D Convection-Diffusion equation with Dirichlet boundary conditions 
using FVM
"""
# Preliminaries
import numpy as np
import matplotlib.pyplot as plt


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

# Volumetric source. Depends on cell centroid values
S = 100.*(1. + 5*cxx + 5*cyy)*cell_volume


# Vector of unknown temperatures
T = np.ones((len(cxx),1)) # Length should be equal to no. of cell centroids; 
                           # boundary values lumped with constants in b vector

                              
# Boundary conditions for T (included with constants in b vector)
T_right = 500.    # T in K on right boundary face
T_left = 500.     # T in K on left boundary face
T_top = 500.      # T in K on top boundary face
T_inf = 300.      # T in K of flow across bottom face

mixed_bc_coeff = (h*k/del_y_b)/(h + k/del_y_b) # Coeff for bottom face with mixed BC

# Vector of constants
b = np.zeros((len(cxx),1)) # Length = no. of cell centroids

# Populate the vector of constants
b = np.copy(S.reshape((len(S),1)))     # Source term contribution at every cell centroid

# Boundary face indices
left_face_indices = np.array([i for i in range(0,len(b),numcells)])
right_face_indices = np.array([i for i in range(numcells-1,len(b),numcells)])
top_face_indices = np.array([i for i in range(len(b)-numcells,len(b))])
bottom_face_indices = np.array([i for i in range(0,numcells)])

# Contributions to b from boundaries 
b[bottom_face_indices] += mixed_bc_coeff*area_y*T_inf # Mixed BC contribution on bottom face
b[top_face_indices] += k*area_y*T_top/del_y_b # Dirichlet boundary on top face
b[left_face_indices] += k*area_x*T_left/del_x_b # Dirichlet boundary on left faces
b[right_face_indices] += k*area_x*T_right/del_x_b # Dirichlet boundary on right faces
  

# Matrix of coefficients (size = size(b)*size(T))
A = np.zeros((len(b),len(b))) # Each row corresonds to a cell centroid

# Populate A. It is sparse and diagonally dominant, so a single loop should suffice.
for i in range(len(b)): # Fill from bottom to top
          # Bottom boundary cells including corners
          if i in bottom_face_indices: 
               if i == 0: # Bottom left cell
                    A[i,i+1] = -k*area_x/del_x # East cell
                    A[i,i+numcells] = -k*area_y/del_y # North cell
                    A[i,i] = abs(A[i,i+1]) + abs(A[i,i+numcells]) + mixed_bc_coeff*area_y \
                         + k*area_x/del_x_b
               elif i == numcells-1: # Bottom right cell
                    A[i,i-1] = -k*area_x/del_x # West cell
                    A[i,i+numcells] = -k*area_y/del_y # North cell
                    A[i,i] = abs(A[i,i-1]) + abs(A[i,i+numcells]) + mixed_bc_coeff*area_y \
                         + k*area_x/del_x_b
               else: # Remaining bottom cells
                    A[i,i+1] = -k*area_x/del_x # East cell
                    A[i,i-1] = -k*area_x/del_x # West cell
                    A[i,i+numcells] = -k*area_y/del_y # North cell
                    A[i,i] = abs(A[i,i+1]) + abs(A[i,i-1]) + abs(A[i,i+numcells]) + \
                         mixed_bc_coeff*area_y
          # Top boundary cells including corners
          elif i in top_face_indices: 
               if i == len(b)-numcells: # Top left cell
                    A[i,i+1] = -k*area_x/del_x # East cell
                    A[i,i-numcells] = -k*area_y/del_y # South cell
                    A[i,i] = abs(A[i,i+1]) + abs(A[i,i-numcells]) + k*area_y/del_y_b \
                         + k*area_x/del_x_b
               if i == len(b)-1: # Top right cell
                    A[i,i-1] = -k*area_x/del_x # West cell
                    A[i,i-numcells] = -k*area_y/del_y # South cell
                    A[i,i] = abs(A[i,i-1]) + abs(A[i,i-numcells]) + k*area_y/del_y_b \
                         + k*area_x/del_x_b
               else: # Remaining top faces
                    A[i,i+1] = -k*area_x/del_x # East cell
                    A[i,i-1] = -k*area_x/del_x # West cell
                    A[i,i-numcells] = -k*area_y/del_y # South cell
                    A[i,i] = abs(A[i,i-1]) + abs(A[i,i+1]) + abs(A[i,i-numcells]) \
                         + k*area_x/del_x_b
          # Right boundary cells excluding corners
          elif i in right_face_indices:
               if i != numcells-1 and i != len(b)-1:
                    A[i,i-1] = -k*area_x/del_x # West cell
                    A[i,i+numcells] = -k*area_y/del_y # North cell
                    A[i,i-numcells] = -k*area_y/del_y # South cell
                    A[i,i] = abs(A[i,i-1]) + abs(A[i,i+numcells]) + abs(A[i,i-numcells]) \
                         + k*area_x/del_x_b
          # Left boundary cells excluding corners
          elif i in left_face_indices: 
               if i != 0 and i != len(b)-numcells: 
                    A[i,i+1] = -k*area_x/del_x # East cell
                    A[i,i+numcells] = -k*area_y/del_y # North cell
                    A[i,i-numcells] = -k*area_y/del_y # South cell
                    A[i,i] = abs(A[i,i+1]) + abs(A[i,i+numcells]) + abs(A[i,i-numcells]) \
                         + k*area_x/del_x_b
          # Interior cells
          else: 
                    A[i,i-1] = -k*area_x/del_x # West cell
                    A[i,i+1] = -k*area_x/del_x # East cell
                    A[i,i+numcells] = -k*area_y/del_y # North cell
                    A[i,i-numcells] = -k*area_y/del_y # South cell
                    A[i,i] = abs(A[i,i-1]) + abs(A[i,i+1]) + abs(A[i,i+numcells]) + \
                         abs(A[i,i-numcells]) # Diagonal element


######## Gauss-Seidel iterative method ########
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

########################

# Solve the system 
T_fvm = gauss(A, b, T, 1000)
Z = T_fvm.reshape((numcells,numcells))

# Print the result
print("\n")
print("Vector of constants is:\n")
print(b.reshape((numcells,numcells)))
print("\n")
print("Coefficient matrix is:\n")
print(A)
print("\n")
print("---------- Solution ----------")
print("\n")
print("T values from FVM method:\n")
print(Z.round(2))

# Plot temperature along central x and y axes
x_plot = np.zeros(len(cc_x)+2) # To add boundary points for plotting
x_plot[1:-1] = cc_x
x_plot[len(x_plot)-1] = domain_length

y_plot = np.zeros(len(cc_y)+2) # To add boundary points for plotting
y_plot[1:-1] = cc_y
y_plot[len(y_plot)-1] = domain_length

T_plot_x = np.hstack((500.,Z[int(Z.shape[0]/2),:],500.))

plt.figure(num=1,figsize=(12,12))
plt.title(r'Temperature profile across the refractory brick',fontsize=16)
plt.plot(x_plot,T_plot_x,'r',label='Central x-axis',linewidth=2)
plt.plot(cc_x,Z[:,int(Z.shape[1]/2)],'b',label='Central y-axis',linewidth=2)
plt.grid(axis='both')
plt.xlabel(r'$x, y$ [m]',fontsize=16)
plt.ylabel(r'$T [K]$',fontsize=16)
plt.legend(loc='best',fontsize=16)

# Contour plot
level_curves=[i for i in range(400,700,20)]
plt.figure(num=2,figsize=(16,16))
plt.xlabel(r'$x [m]$',fontsize=16)
plt.ylabel(r'$y [m]$',fontsize=16)
plt.title("Temperature Profile on the refractory brick")
c = plt.pcolormesh(cc_x, cc_y, Z, cmap='nipy_spectral', vmin=Z.min(), vmax=Z.max(), \
                   shading='auto')
CS = plt.contour(cc_x,cc_y,Z,level_curves,colors='white')
plt.clabel(CS, inline=True, fmt='%1.0f',fontsize=16)
plt.colorbar(c)