#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 12 20:40:19 2021

@author: dileepn

Solution of 1D, steady, incompressible flow through a porous medium on a staggered 
grid along with SIMPLE scheme for pressure-velocity coupling. 
 
"""
# Preliminaries
import numpy as np
#import matplotlib.pyplot as plt

# Given parameters
uLB = 1.0 # Velocity at left boundary (inflow) in m/s
uRB = 1.0 # Velocity at right boundary (outflow) in m/s
rho = 1000 # Density of fluid in kg/m3
C = 0.1 # Constant of porosity in 1/m
L = 0.3 # Length of the domain in m

# Initial guess for SIMPLE scheme
uA = 0.5
uB = 0.5
p1 = 0.0
p2 = 0.0
p3 = 0.0

# Geometry of flow domain
numcells = 3
dx = L/numcells

# Under-relaxation and tolerance for convergence
tolerance = 1.e-6 
alphaU = 1     # URF for x-velocity (no relaxation if set equal to 1)
alphaP = 1     # URF for pressure (no relaxation if set equal to 1)

# Maximum iterations permitted
maxiter = 100

# Begin the main loop 
for i in range(maxiter):
     
     # Calculate/update momentum coefficients
     aA = C*rho*uA*dx
     dA = 1./aA
     aB = C*rho*uB*dx
     dB = 1./aB
     
     # Constants, b
     bA = 0.5*rho*C*uA**2*dx
     bB = 0.5*rho*C*uB**2*dx
     
     # Calculate residual for the momentum equations
     u_residual = abs(uA*aA/alphaU - ((p1 - p2) + bA + uA*aA*(1.0 - alphaU)/alphaU)) + \
 		     abs(uB*aB/alphaU - ((p2 - p3) + bB + uB*aB*(1.0 - alphaU)/alphaU))
     u_residual = u_residual/(abs(uA*aA/alphaU + uB*aB/alphaU)) # Normalize to take care of roundoff errors
     
     # Solve momentum equations (under-relaxed)
     uA = ((p1 - p2) + bA + uA*aA*(1.0 - alphaU)/alphaU)*alphaU/aA
     uB = ((p2 - p3) + bB + uB*aB*(1.0 - alphaU)/alphaU)*alphaU/aB
     
     # Calculate residual for the continuity equations
     c_residual = (abs(uA - uLB) + abs(uB - uA) + abs(uRB - uB))
     c_residual = c_residual/(0.5*((abs(uA) + abs(uB - uA) + abs(uB)))) # Normalize
     
     # Check for convergence
     if (u_residual + c_residual < tolerance):
          print("Converged solution is: p1 = {:.2f}, p2 = {:.2f}, p3 = {:.2f}, uA = {:.2f}, uB = {:.2f}\n".format(p1, p2, p3, uA, uB))
          break
                
     
     # Solve the linear system for pprimes 
     coeffMatrix = np.array([[dA, -dA], [-dA, dA + dB]])
     
     # Gauss-Seidel to iteratively the linear system to get pprimes
     itr = 100
     tol = 1.e-4
     x = np.zeros((numcells-1,1))
     for i in range(itr):
          xprev = x
          a = uLB - uA
          b = uA - uB + dB*p3
          bVector = np.array([[a],[b]])
          L = np.tril(coeffMatrix)
          U = coeffMatrix - L
          x = np.dot(np.linalg.inv(L), bVector - np.dot(U, x))
          p3 = x[1] - (uRB - uB)/dB
          if np.sqrt(np.sum((x - xprev)**2)) < tol:
               print("\n")
               print("GS converged in",i,"iterations.")
               print("pprimes for iteration {} are: {:.2f}, {:.2f}, and {:.2f}".format(i, x[0], x[1], p3))
               break
     
     p1prime = x[0]
     p2prime = x[1]
     p3prime = p2prime - (uRB - uB)/dB
     
     # Correct the velocities and pressures
     uAprime = dA*(p1prime - p2prime)
     uBprime = dB*(p2prime - p3prime)
     uA = uA + uAprime
     uB = uB + uBprime
     p1 = p1 + alphaP*p1prime
     p2 = p2 + alphaP*p2prime
     p3 = p3 + alphaP*p3prime
     
     # Check if continuity is satisfied 
     cont1 = uA - uLB
     cont2 = uB - uA
     cont3 = uRB - uB
     
     # Output 
     if i == 0:
          print("Iter \t uA \t uB \t p1 \t p2 \t p3 \t u_res \t cont_res")
          
     print("{0:2d} {0:1.3e} {0:1.3e} {0:1.3e} {0:1.3e} {0:1.3e} {0:1.3e} {0:1.3e}".format(i, uA, uB, p1, p2, p3, u_residual, c_residual))  