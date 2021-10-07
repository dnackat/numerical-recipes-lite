#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct  6 22:25:40 2021

@author: dileepn

Solution 2D Convection-Diffusion equation with Dirichlet boundary conditions 
using FVM on a uniform, Cartesian mesh.

Governing PDE: div(Gamma*grad(phi)) - div(rho*vel_vector*phi) + Sources = 0
"""
# Preliminaries
import numpy as np
#import matplotlib.pyplot as plt


# Parameters
k = 0.0    # Constant thermal conductivity in W/m-K
rho = 2.0  # Density of the fluid in kg/m3
h = 0.0   # Constant convective heat transfer coefficient in W/m2-K
domain_length = 3.0 # In both x and y directions in m

# Create the grid
numcells = 3   # Number of cells in the x and y directions
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

# Velocity field at face centroids
u = 1. + fxx**2
v = 1. + fyy**2

# Volumetric source. Depends on cell centroid values
S = 0.5*4*(cxx + cyy)*cell_volume

# Fluxes
F_i = rho*area_x*u
F_j = rho*area_y*v

# Vector of unknown temperatures
phi = np.ones((len(cxx),1)) # Length should be equal to no. of cell centroids; 
                           # boundary values lumped with constants in b vector

                              
# Boundary conditions for phi (included with constants in b vector)
phi_right = 0.    # phi on right boundary face
phi_left = 0.     # phi on left boundary face
phi_top = 0.      # phi on top boundary face
phi_bottom = 1.      # phi of flow across bottom face

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
b[bottom_face_indices] += F_j[0]*phi_bottom # Dirichlet BC contribution on bottom face
b[top_face_indices] += k*area_y*phi_top/del_y_b # Dirichlet boundary on top face
b[left_face_indices] += F_i[0]*phi_left + k*area_x*phi_left/del_x_b # Dirichlet boundary on left faces
b[right_face_indices] += k*area_x*phi_right/del_x_b # Dirichlet boundary on right faces
  

# Matrix of coefficients (size = size(b)*size(phi))
A = np.zeros((len(b),len(b))) # Each row corresonds to a cell centroid

#%%
# Populate A. It is sparse and diagonally dominant, so a single loop should suffice.
for i in range(len(b)): # Fill from bottom to top
          # Bottom boundary cells including corners
          if i in bottom_face_indices: 
               if i == 0: # Bottom left cell
                    A[i,i+1] = -max(-F_i[i+1],0.) - k*area_x/del_x # East cell
                    A[i,i+numcells] = -max(-F_j[i+numcells+1],0.) - k*area_y/del_y # North cell
                    A[i,i] = abs(A[i,i+1]) + abs(A[i,i+numcells]) \
                         + k*area_x/del_x_b + \
                         (F_i[i+1] + F_j[i+numcells+1])  
               elif i == numcells-1: # Bottom right cell
                    A[i,i-1] = -max(F_i[i],0.) - k*area_x/del_x # West cell
                    A[i,i+numcells] = -max(-F_j[i+numcells+1],0.) - k*area_y/del_y # North cell
                    A[i,i] = abs(A[i,i-1]) + abs(A[i,i+numcells]) \
                         + k*area_x/del_x_b + \
                         (F_i[i+1] - F_i[i] + F_j[i+numcells+1])
               else: # Remaining bottom cells
                    A[i,i+1] = -max(-F_i[i+1],0.) - k*area_x/del_x # East cell
                    A[i,i-1] = -max(F_i[i],0.) - k*area_x/del_x # West cell
                    A[i,i+numcells] = -max(-F_j[i+numcells+1],0.) - k*area_y/del_y # North cell
                    A[i,i] = abs(A[i,i+1]) + abs(A[i,i-1]) + abs(A[i,i+numcells]) + \
                    (F_i[i+1] - F_i[i] + F_j[i+numcells+1])
          # Top boundary cells including corners
          elif i in top_face_indices: 
               if i == len(b)-numcells: # Top left cell
                    A[i,i+1] = -max(-F_i[i+numcells],0.) - k*area_x/del_x # East cell
                    A[i,i-numcells] = -max(F_j[i+numcells-1],0.) -k*area_y/del_y # South cell
                    A[i,i] = abs(A[i,i+1]) + abs(A[i,i-numcells]) + k*area_y/del_y_b \
                         + k*area_x/del_x_b + \
                         (F_i[i+numcells-1] + F_j[i+numcells+2] - F_j[i+2])
               if i == len(b)-1: # Top right cell
                    A[i,i-1] = -max(F_i[i+numcells-1],0.) - k*area_x/del_x # West cell
                    A[i,i-numcells] = -max(F_j[i+numcells-1],0.) - k*area_y/del_y # South cell
                    A[i,i] = abs(A[i,i-1]) + abs(A[i,i-numcells]) + k*area_y/del_y_b \
                         + k*area_x/del_x_b + \
                         (F_i[i+numcells] - F_i[i+numcells-1] + F_j[i+numcells+3] - F_j[i+2])
               else: # Remaining top faces
                    A[i,i+1] = -max(-F_i[i+numcells],0.) - k*area_x/del_x # East cell
                    A[i,i-1] = -max(F_i[i+numcells-1],0.) - k*area_x/del_x # West cell
                    A[i,i-numcells] = -max(F_j[i+numcells-1],0.) - k*area_y/del_y # South cell
                    A[i,i] = abs(A[i,i-1]) + abs(A[i,i+1]) + abs(A[i,i-numcells]) \
                         + k*area_x/del_x_b + \
                         (F_i[i+numcells] - F_i[i+numcells-1] + F_j[i+numcells+3] - F_j[i+2])
          # Right boundary cells excluding corners
          elif i in right_face_indices:
               if i != numcells-1 and i != len(b)-1:
                    A[i,i-1] = -max(F_i[i+1],0.) - k*area_x/del_x # West cell
                    A[i,i+numcells] = -max(-F_j[i+numcells+2],0.) - k*area_y/del_y # North cell
                    A[i,i-numcells] = -max(F_j[i+numcells-1],0.) - k*area_y/del_y # South cell
                    A[i,i] = abs(A[i,i-1]) + abs(A[i,i+numcells]) + abs(A[i,i-numcells]) \
                         + k*area_x/del_x_b + \
                         (F_i[i+2] - F_i[i+1] + F_j[i+numcells+2] - F_j[i+1])
          # Left boundary cells excluding corners
          elif i in left_face_indices: 
               if i != 0 and i != len(b)-numcells: 
                    A[i,i+1] = -max(-F_i[i+2],0.) - k*area_x/del_x # East cell
                    A[i,i+numcells] = -max(-F_j[i+numcells+2],0.) - k*area_y/del_y # North cell
                    A[i,i-numcells] = -max(F_j[i+numcells-1],0.) - k*area_y/del_y # South cell
                    A[i,i] = abs(A[i,i+1]) + abs(A[i,i+numcells]) + abs(A[i,i-numcells]) \
                         + k*area_x/del_x_b + \
                         (F_i[i+2] - F_i[i+1] + F_j[i+numcells+2] - F_j[i+1])
          # Interior cells
          else: 
                    A[i,i-1] = -max(F_i[i+1],0.) - k*area_x/del_x # West cell
                    A[i,i+1] = -max(-F_i[i+2],0.) - k*area_x/del_x # East cell
                    A[i,i+numcells] = -max(-F_j[i+numcells+2],0.) - k*area_y/del_y # North cell
                    A[i,i-numcells] = -max(F_j[i+numcells-1],0.) - k*area_y/del_y # South cell
                    A[i,i] = abs(A[i,i-1]) + abs(A[i,i+1]) + abs(A[i,i+numcells]) + \
                         abs(A[i,i-numcells]) + \
                         (F_i[i+2] - F_i[i+1] + F_j[i+numcells+2] - F_j[i+1])  # Diagonal element


######## Gauss-Seidel iterative method ########
def gauss(A, b, x, n):
     """ phihis is a function that uses the Gauss-Seidel iterative scheme from 
         numerical linear algebra to solve a linear system of the form 
         Ax = b. phihe scheme updates x using the following algorithm:
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

################################################

#Solve the system 
phi_fvm = gauss(A, b, phi, 1000)
Z = phi_fvm.reshape((numcells,numcells))

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
print("phi values from FVM method:\n")
print(Z.round(2))

