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

# Under-relaxation and tolerance for convergence
tolerance = 1.e-6 
alphaU = 1     # URF for x-velocity
alphaP = 1     # URF for pressure

# Maximum iterations permitted
maxiter = 100 

# Begin the main loop 

 