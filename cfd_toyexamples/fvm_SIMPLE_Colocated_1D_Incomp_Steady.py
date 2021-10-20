#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 20 22:39:07 2021

@author: dileepn
"""
import numpy as np
#import matplotlib.pyplot as plt

# Given parameters
u1 = 1.0 # Velocity at left boundary (inflow) in m/s
p4 = 10.0 # Pressure at right boundary (outflow) in N/m2
rho = 1000 # Density of fluid in kg/m3
C = 0.1 # Constant of porosity in 1/m
L = 0.3 # Length of the domain in m

# Initial guess for SIMPLE scheme
uA = 0.5
u2 = 0.5
uB = 0.5
u3 = 0.5
uC = 0.5
u4 = 0.5
p1 = 10.0
pA = 10.0
p2 = 10.0
pB = 10.0
p3 = 10.0
pC = 10.0

# Only velocity BCs are given, so set one of the pressue corrections to zero
p3prime = 0.0

# Geometry of flow domain
numcells = 3
dx = L/numcells

# Under-relaxation and tolerance for convergence
tolerance = 1.e-6 
alphaU = 1.0    # URF for x-velocity (no relaxation if set equal to 1)
alphaP = 1.0     # URF for pressure (no relaxation if set equal to 1)

# Maximum iterations permitted
maxiter = 100

# Begin the main loop 
for i in range(maxiter):
     
     # Calculate/update momentum coefficients
     aA = C*rho*uA*dx/alphaU
     bA = 0.5*rho*C*uA**2*dx + (1.0 - alphaU)*aA*uA
     dA = 1./aA
     
     aB = C*rho*uB*dx/alphaU
     bB = 0.5*rho*C*uB**2*dx + (1.0 - alphaU)*aB*uB
     dB = 1./aB
     
     aC = C*rho*uC*dx/alphaU
     bC = 0.5*rho*C*uC**2*dx + (1.0 - alphaU)*aC*uC
     dC = 1./aC
     
     a1 = C*rho*u1*dx/(2.0*alphaU)
     b1 = 0.5*rho*C*u1**2*dx + (1.0 - alphaU)*a1*u1
     d1 = 1./a1
     
     a4 = C*rho*u4*dx/(2.0*alphaU)
     b4 = 0.5*rho*C*u4**2*dx + (1.0 - alphaU)*a4*u4
     d4 = 1./a4
     
     # Calculate residual for the momentum equations
     u_residual = abs(uA*aA - (p1 - p2) - bA) + \
      		      abs(uB*aB - (p2 - p3) - bB) + \
                  abs(uC*aC - (p3 - p4) - bC) 
     u_residual = u_residual/(abs(uA*aA + uB*aB + uC*aC)) # Normalize to take care of roundoff errors
     
     # Solve momentum equations (under-relaxed)
     uA = ((p1 - p2) + bA)/aA
     uB = ((p2 - p3) + bB)/aB
     uC = ((p3 - p4) + bC)/aC
     
     # Calculate hat velocities
     uAhat = bA/aA
     uBhat = bB/aB
     uChat = bC/aC
     
     u2hat = (uAhat + uBhat)/2.0
     u3hat = (uBhat + uChat)/2.0
     
     u1hat = b1/a1
     u4hat = b4/a4
     
     d2 = (dA + dB)/2.0
     d3 = (dB + dC)/2.0
     
     # Using Rhie-Chow interpolation, calculate face star velocities
     
     # Calculate residual for the continuity equations
     c_residual = (abs(uA - uLB) + abs(uB - uA) + abs(uRB - uB))
     c_residual = c_residual/(0.5*((abs(uA) + abs(uB - uA) + abs(uB)))) # Normalize
     
     # Check for convergence
     if (u_residual + c_residual < tolerance):
          print("\n")
          print("=============================================================================")
          print("Converged solution is: p1 = {:.2f}, p2 = {:.2f}, p3 = {:.2f}, uA = {:.2f}, uB = {:.2f}\n".format(p1, p2, p3, uA, uB))
          print("=============================================================================")
          print("\n")
          break
                
     
     # Solve the linear system for pprimes
     coeffMatrix = np.array([[dA, -dA], [-dA, dA + dB]]) 
     bVector = np.array([[uLB - uA],[uA - uB + dB*p3]])
     
     p1prime, p2prime = np.linalg.solve(coeffMatrix, bVector)
     p1prime = p1prime.item()
     p2prime = p2prime.item()
     
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
          print("----------------------------------------------------------------------")
          print("It |  uA \t uB \t    p1 \t   p2 \t  p3  \t    u_res       cont_res")
          print("----------------------------------------------------------------------")
          
     print("{:2d} | {:.3f} |  {:.3f} |  {:.3f} |{:.3f} |{:.3f} |  {:1.3e} |  {:1.3e}".format(i, uA, uB, p1, p2, p3, u_residual, c_residual))  