#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 20 22:39:07 2021

@author: dileepn

Solution of 1D, steady, incompressible flow through a porous medium on a colocated  
grid along with SIMPLE scheme for pressure-velocity coupling along with Rhie-Chow 
momentum interpolation.
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

# Geometry of flow domain
numcells = 3
dx = L/numcells

# Under-relaxation and tolerance for convergence
tolerance = 1.e-6 
alphaU = 0.8    # URF for x-velocity (no relaxation if set equal to 1)
alphaP = 0.7     # URF for pressure (no relaxation if set equal to 1)

# Maximum iterations permitted
maxiter = 100
u_residual = 10.0
c_residual = 10.0

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
     b1 = 0.5*rho*C*u1**2*(dx/2.0) + (1.0 - alphaU)*a1*u1
     d1 = 1./a1
     
     a4 = C*rho*u4*dx/(2.0*alphaU)
     b4 = 0.5*rho*C*u4**2*(dx/2.0) + (1.0 - alphaU)*a4*u4
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
     u2 = u2hat + d2*(pA - pB)
     u3 = u3hat + d3*(pB - pC)
     u4 = u4hat + d4*(pC - p4)
     
     # Check for convergence
     if (u_residual + c_residual < tolerance):
          print("\n")
          print("=============================================================================")
          print("The solution is converged in",i,"iterations.\n")
          print("Converged solution is:\n uA = {:.2f}, uB = {:.2f}, uC = {:.2f}, \n u1 = {:.2f}, u2 = {:.2f}, u3 = {:.2f}, u4 = {:.2f}, \n p1 = {:.2f}, p4 = {:.2f}, \n pA = {:.2f}, pB = {:.2f}, pC = {:.2f}, \n".format(uA, uB, uC, u1, u2, u3, u4, p1, p4, pA, pB, pC))
          print("=============================================================================")
          print("\n")
          break

     # Solve the linear system for pprimes
     coeffMatrix = np.array([[d2, -d2, 0.],[-d2, d2 + d3, -d3],[0., -d3, d3 + d4]]) 
     bVector = np.array([[u1 - u2],[u2 - u3],[u3 - u4]])
     
     pAprime, pBprime, pCprime = np.linalg.solve(coeffMatrix, bVector)
     pAprime = pAprime.item()
     pBprime = pBprime.item()
     pCprime = pCprime.item()
	 
	 # Update boundary pressure corrections
     p1prime = pAprime
     p2prime = (pAprime + pBprime)/2.0
     p3prime = (pBprime + pCprime)/2.0
	 
	 # Correct face velocities
     u2prime = d2*(pAprime - pBprime)
     u2 = u2 + u2prime
     
     u3prime = d3*(pBprime - pCprime)
     u3 = u3 + u3prime
     
     u4prime = d4*pCprime # p4 is given, so p4prime is zero
     u4 = u4 + u4prime
	
	 # Correct cell pressures
     pA = pA + alphaP*pAprime
     pB = pB + alphaP*pBprime
     pC = pC + alphaP*pCprime
	
	 # Correct cell velocities to improve convergence
     uA = uA + dA*(p1prime - p2prime)
     uB = uB + dB*(p2prime - p3prime)
     uC = uC + dC*p3prime # p4 is given, so p4prime is zero
    	
     # Update face pressures from boundary-cell momentum equations
     p1 = pA + (u1 - u1hat)/d1
     p2 = (pA + pB)/2.0
     p3 = (pB + pC)/2.0
     
     # Calculate continuity residual
     c_residual = abs(u1 - u2) + abs(u2 - u3) + abs(u3 - u4)
     c_residual = c_residual/(0.5*(u1 + u2 + u3))
     
     # Output 
     if i == 0:
          print("----------------------------------------------------------------------")
          print("It |  uA \t uB \t uC \t  pA \t   pB \t  pC  \t    u_res       cont_res")
          print("----------------------------------------------------------------------")
          
     print("{:2d} | {:.3f} |  {:.3f} |  {:.3f} | {:.3f} |{:.3f} |{:.3f} |  {:1.3e} |  {:1.3e}".format(i, uA, uB, uC, pA, pB, pC, u_residual, c_residual))  