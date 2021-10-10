#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct  6 22:25:40 2021

@author: dileepn

Solution 2D Convection-Diffusion equation with Dirichlet boundary conditions 
using FVM on a uniform, Cartesian mesh. 

Governing PDE: div(Gamma*grad(phi)) - div(rho*vel_vector*phi) + Sources = 0

Assumptions: 1. Linear profile for the gradient
             2. First-order Upwind Difference Scheme (UDS) for the convection term  
"""
# Preliminaries
import numpy as np
import matplotlib.pyplot as plt


# Parameters
k = 3.  # Constant thermal conductivity in W/m-K
rho = 2.0  # Density of the fluid in kg/m3
h = 0.0   # Constant convective heat transfer coefficient in W/m2-K
domain_length = 3.0 # In both x and y directions in m

# Create the grid
numcells = 5   # Number of cells in the x and y directions
numfaces = numcells + 1 # Number of faces in the x and y directions
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
u = 1. + fc_x**2    # x-component of velocity 
v = 1. + fc_y**2    # y-component of velocity

# Volumetric source. Depends on cell centroid values
S = 0.5*4*(cxx + cyy)*cell_volume

# Fluxes through faces
F_i = rho*area_x*u
F_j = rho*area_y*v

# Vector of unknown temperatures
phi = np.ones((len(cxx),1)) # Length should be equal to no. of cell centroids; 
                           # boundary values lumped with constants in b vector

                              
# Boundary conditions for phi (included with constants in b vector)
phi_right = 0.    # phi on right boundary face
phi_left = 0.     # phi on left boundary face
phi_top = 0.      # phi on top boundary face
phi_bottom = 1.   # phi of flow across bottom face

# Vector of constants
b = np.zeros((len(cxx),1)) # Length = no. of cell centroids

# Populate the vector of constants
b = np.copy(S.reshape((len(S),1)))     # Source term contribution at every cell centroid

# Boundary cell indices
left_cell_indices = np.array([i for i in range(0,len(b),numcells)])
right_cell_indices = np.array([i for i in range(numcells-1,len(b),numcells)])
top_cell_indices = np.array([i for i in range(len(b)-numcells,len(b))])
bottom_cell_indices = np.array([i for i in range(0,numcells)])

# Contributions to b from boundaries 
b[bottom_cell_indices] += F_j[0]*phi_bottom + k*area_y*phi_bottom/del_y_b # Dirichlet BC contribution on bottom face
b[top_cell_indices] += k*area_y*phi_top/del_y_b # Dirichlet boundary on top face
b[left_cell_indices] += F_i[0]*phi_left + k*area_x*phi_left/del_x_b # Dirichlet boundary on left faces
b[right_cell_indices] += k*area_x*phi_right/del_x_b # Dirichlet boundary on right faces
  

# Matrix of coefficients (size = size(b)*size(phi))
A = np.zeros((len(b),len(b))) # Each row corresonds to a cell centroid

##### Function to determine face indices based on cell centroids #####
def faceid(a, name="west"):
     
     # What row am I on?
     rowid = int(np.floor(a/numcells))
     
     # What column am I on?
     colid = int(np.mod(a, numcells))

     if name == "west":
          return colid
     elif name == "east": 
          return colid + 1
     elif name == "south":
          return rowid
     elif name == "north":
          return rowid + 1
     else:
          print("Invalid face index! \n")
          return None
     
     return None

#######################################################################

# Populate A. It is sparse and diagonally dominant, so a single loop should suffice.
for i in range(len(b)): # Fill from bottom to top
     # Bottom boundary cells including corners
     if i in bottom_cell_indices: 
           if i == 0: # Bottom left cell
               A[i,i+1] = -max(-F_i[faceid(i,"east")],0.) - k*area_x/del_x # East cell
               A[i,i+numcells] = -max(-F_j[faceid(i,"north")],0.) - k*area_y/del_y # North cell
               A[i,i] = abs(A[i,i+1]) + abs(A[i,i+numcells]) \
                     + k*area_x/del_x_b + k*area_y/del_y_b + \
                     (F_i[faceid(i,"east")] + F_j[faceid(i,"north")])  
           elif i == numcells-1: # Bottom right cell
               A[i,i-1] = -max(F_i[faceid(i,"west")],0.) - k*area_x/del_x # West cell
               A[i,i+numcells] = -max(-F_j[faceid(i,"north")],0.) - k*area_y/del_y # North cell
               A[i,i] = abs(A[i,i-1]) + abs(A[i,i+numcells]) \
                     + k*area_x/del_x_b + k*area_y/del_y_b + \
                     (F_i[faceid(i,"east")] - F_i[faceid(i,"west")] + F_j[faceid(i,"north")])
           else: # Remaining bottom cells
               A[i,i+1] = -max(-F_i[faceid(i,"east")],0.) - k*area_x/del_x # East cell
               A[i,i-1] = -max(F_i[faceid(i,"west")],0.) - k*area_x/del_x # West cell
               A[i,i+numcells] = -max(-F_j[faceid(i,"north")],0.) - k*area_y/del_y # North cell
               A[i,i] = abs(A[i,i+1]) + abs(A[i,i-1]) + abs(A[i,i+numcells]) \
                     + k*area_y/del_y_b + \
               (F_i[faceid(i,"east")] - F_i[faceid(i,"west")] + F_j[faceid(i,"north")])
     # Top boundary cells including corners
     elif i in top_cell_indices: 
           if i == len(b)-numcells: # Top left cell
               A[i,i+1] = -max(-F_i[faceid(i,"east")],0.) - k*area_x/del_x # East cell
               A[i,i-numcells] = -max(F_j[faceid(i,"south")],0.) - k*area_y/del_y # South cell
               A[i,i] = abs(A[i,i+1]) + abs(A[i,i-numcells]) + k*area_y/del_y_b \
                     + k*area_x/del_x_b + \
                     (F_i[faceid(i,"east")] + F_j[faceid(i,"north")] - F_j[faceid(i,"south")])
           if i == len(b)-1: # Top right cell
               A[i,i-1] = -max(F_i[faceid(i,"west")],0.) - k*area_x/del_x # West cell
               A[i,i-numcells] = -max(F_j[faceid(i,"south")],0.) - k*area_y/del_y # South cell
               A[i,i] = abs(A[i,i-1]) + abs(A[i,i-numcells]) + k*area_y/del_y_b \
                     + k*area_x/del_x_b + \
                     (F_i[faceid(i,"east")] - F_i[faceid(i,"west")] + F_j[faceid(i,"north")] - F_j[faceid(i,"south")])
           else: # Remaining top faces
               A[i,i+1] = -max(-F_i[faceid(i,"east")],0.) - k*area_x/del_x # East cell
               A[i,i-1] = -max(F_i[faceid(i,"west")],0.) - k*area_x/del_x # West cell
               A[i,i-numcells] = -max(F_j[faceid(i,"south")],0.) - k*area_y/del_y # South cell
               A[i,i] = abs(A[i,i-1]) + abs(A[i,i+1]) + abs(A[i,i-numcells]) \
                     + k*area_y/del_y_b + \
                     (F_i[faceid(i,"east")] - F_i[faceid(i,"west")] + F_j[faceid(i,"north")] - F_j[faceid(i,"south")])
     # Right boundary cells excluding corners
     elif i in right_cell_indices:
           if i != numcells-1 and i != len(b)-1:
               A[i,i-1] = -max(F_i[faceid(i,"west")],0.) - k*area_x/del_x # West cell
               A[i,i+numcells] = -max(-F_j[faceid(i,"north")],0.) - k*area_y/del_y # North cell
               A[i,i-numcells] = -max(F_j[faceid(i,"south")],0.) - k*area_y/del_y # South cell
               A[i,i] = abs(A[i,i-1]) + abs(A[i,i+numcells]) + abs(A[i,i-numcells]) \
                     + k*area_x/del_x_b + \
                     (F_i[faceid(i,"east")] - F_i[faceid(i,"west")] + F_j[faceid(i,"north")] - F_j[faceid(i,"south")])
     # Left boundary cells excluding corners
     elif i in left_cell_indices: 
           if i != 0 and i != len(b)-numcells: 
               A[i,i+1] = -max(-F_i[faceid(i,"east")],0.) - k*area_x/del_x # East cell
               A[i,i+numcells] = -max(-F_j[faceid(i,"north")],0.) - k*area_y/del_y # North cell
               A[i,i-numcells] = -max(F_j[faceid(i,"south")],0.) - k*area_y/del_y # South cell
               A[i,i] = abs(A[i,i+1]) + abs(A[i,i+numcells]) + abs(A[i,i-numcells]) \
                     + k*area_x/del_x_b + \
                     (F_i[faceid(i,"east")] + F_j[faceid(i,"north")] - F_j[faceid(i,"south")])
     # Interior cells
     else: 
               A[i,i-1] = -max(F_i[faceid(i,"west")],0.) - k*area_x/del_x # West cell
               A[i,i+1] = -max(-F_i[faceid(i,"east")],0.) - k*area_x/del_x # East cell
               A[i,i+numcells] = -max(-F_j[faceid(i,"north")],0.) - k*area_y/del_y # North cell
               A[i,i-numcells] = -max(F_j[faceid(i,"south")],0.) - k*area_y/del_y # South cell
               A[i,i] = abs(A[i,i-1]) + abs(A[i,i+1]) + abs(A[i,i+numcells]) + \
                     abs(A[i,i-numcells]) + \
                     (F_i[faceid(i,"east")] - F_i[faceid(i,"west")] + F_j[faceid(i,"north")] - F_j[faceid(i,"south")])  # Diagonal element


######## Gauss-Seidel iterative method ########
def gauss(A, b, x, itr=1000, tol=1.e-6):
     """ phihis is a function that uses the Gauss-Seidel iterative scheme from 
         numerical linear algebra to solve a linear system of the form 
         Ax = b. phihe scheme updates x using the following algorithm:
         x(k+1) = inv(L)*(b - U*x(k)), where U and L are upper and lower
         triangular matrices obtained by performing LU decmposition of A."""

     L = np.tril(A)
     U = A - L
     for i in range(itr):
        xprev = x
        x = np.dot(np.linalg.inv(L), b - np.dot(U, x))
        if np.linalg.norm(x - xprev) < tol:
            print("\n")
            print("Converged in",i+1,"iterations.")
            break
     return x

################################################

# Solve the system 
phi_fvm = gauss(A, b, phi, 1000).reshape((numcells, numcells))

# Invert phi so that it corresponds with the given geometry
Z = np.flip(phi_fvm, axis=0)

# Print the result
print("\n")
print("Vector of constants, b is:\n")
print(b)
print("\n")
print("Coefficient matrix, A is:\n")
print(A)
print("\n")
print("------------ Solution of the linear system A*phi = b ------------")
print("\n")
print("phi values from FVM method:\n")
print(Z.round(3))

##### Plot the results (without interpolation to interior nodes) #####
cells_x, cells_y = np.meshgrid(cc_x, cc_y)
level_curves = [i for i in np.arange(0.1,1.1,0.1)]

fig, ax = plt.subplots(figsize=(16,16))
#c = ax.pcolormesh(cells_x, cells_y, phi_fvm, cmap='nipy_spectral', \
                  #vmin=phi_fvm.min(), vmax=phi_fvm.max(), shading='auto')
ax.axis([cells_x.min(), cells_x.max(), cells_y.min(), cells_y.max()])
#fig.colorbar(c, ax=ax)
CS = ax.contour(cells_x, cells_y, phi_fvm, level_curves, linewidths=2)
ax.clabel(CS, inline=True, fmt='%1.3f', fontsize=16)
plt.show()


