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
     c_residual = (abs(uA - 1) + abs(uB - uA) + abs(1 - uB))
     c_residual = c_residual/(0.5*((abs(uA) + abs(uB - uA) + abs(uB)))) # Normalize
     
     # Check for convergence
     if (u_residual + u_residual < tolerance):
          print("Converged solution is: p1 = {:.2f} \n \
                p2 = {:.2f} \n p3 = {:.2f} \n uA = {:.2f} \n \
                uB = {:.2f}\n".format(p1, p2, p3, uA, uB))
                
     
     # Solve the linear system for pprimes 
     coeffMatrix = np.array([[dA, -dA, 0], [dA, -(dA + dB), dB], [0, dB, -dB]])
     bVector = np.array([[1 - uA],[uB - uA],[1 - uB]])
     
     x, y, z = np.linalg.solve(coeffMatrix, bVector)
     
     p1prime = x.item()
     p2prime = y.item()
     p3prime = z.item()
     
     # Correct the velocities and pressures
     uAprime = dA*(p1prime - p2prime)
     uBprime = dB*(p2prime - p3prime)
     uA = uA + uAprime
     uB = uB + uBprime
     p1 = p1 + alphaP*p1prime
     p2 = p2 + alphaP*p2prime
     p3 = p3 + alphaP*p3prime
     
     # Check if continuity is satisfied 
     cont1 = uA - 1.
     cont2 = uB - uA
     cont3 = 1. - uB
     
     # Output 
     if i == 0:
          print("Iter \t uA \t uB \t p1 \t p2 \t p3 \t u_res \t cont_res \t cont1 \t cont2 \t cont3")
          
     print("{0:2d} {0:1.3e} {0:1.3e} {0:1.3e} {0:1.3e} {0:1.3e} {0:1.3e} {0:1.3e} {0:1.3e} {0:1.3e} {0:1.3e}".format(i, uA, uB, p1, p2, p3, u_residual, c_residual, cont1, cont2, cont3))  