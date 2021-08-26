#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 26 13:00:54 2021

@author: dileepn

Best fit parabola using least-squares approximation
"""
# Preliminaries
import numpy as np
import matplotlib.pyplot as plt

# Data
data = [(1,4.4),(2,39.41),(3,54.83),(4,50.31)]

def lsParabola(data):
     """ This function fits a parabola of the form ax^2 + bx + c to data of
         the form (xi,yi) for i from 1 to N using the least squares approximation. """
     
     # Get all the coefficients for the linear 3x3 system in a, b, and c
     a1 = 0.; b1 = 0.; c1 = 0.; a2 = 0.; b2 = 0.; c2 = 0.; a3 = 0.; b3 = 0.; c3 = 0.
     d1 = 0.; d2 = 0.; d3 = 0.;
     
     for point in data:
          # First row
          a1 += point[0]**4
          b1 += point[0]**3
          c1 += point[0]**2
          d1 += point[1]*point[0]**2
          # Second row
          a2 += point[0]**3
          b2 += point[0]**2
          c2 += point[0]
          d2 += point[1]*point[0]
          # Third row
          a3 += point[0]**2
          b3 += point[0]
          c3 += 1
          d3 += point[1]
          
     # Construct the system
     A = np.array([[a1,b1,c1],[a2,b2,c2],[a3,b3,c3]])
     b = np.array([[d1],[d2],[d3]])
     
     # Solve the system
     x = np.dot(np.linalg.inv(A),b)
     
     # Plot
     xp = np.linspace(np.min([data[i][0] for i in range(len(data))]), \
                      np.max([data[i][0] for i in range(len(data))]))
     yp = x[0]*xp**2 + x[1]*xp + x[2]
     plt.figure(figsize=(16,16)) 
     plt.xlabel(r'$x$',fontsize=16)
     plt.ylabel(r'$y$',fontsize=16)
     plt.xlim((np.min(xp),np.max(xp)))
     plt.scatter([data[i][0] for i in range(len(data))],[data[i][1] for i in range(len(data))],s=150,color='r')
     plt.grid(axis="both")
     plt.plot(xp,yp,linewidth=2)
     
     # Return the solution
     return x

# Solution
a, b, c = lsParabola(data)