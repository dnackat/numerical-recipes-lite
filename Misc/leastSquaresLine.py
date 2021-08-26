#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 26 15:38:59 2021

@author: dileepn

Best fit line using least-squares approximation
"""
# Preliminaries
import numpy as np
import matplotlib.pyplot as plt

# Data
data = [(0,4.4),(1,39.41),(2,54.83),(3,50.31)]

def lsLine(data):
     """ This function fits a line of the form ax + b to data of
         the form (xi,yi) for i from 1 to N using the least squares approximation. """
     
     # Get all the coefficients for the linear 3x3 system in a, b, and c
     a1 = 0.; b1 = 0.; a2 = 0.; b2 = 0.; d1 = 0.; d2 = 0.
     
     for point in data:
          # First row
          a1 += point[0]**2
          b1 += point[0]
          d1 += point[1]*point[0]
          # Second row
          a2 += point[0]
          b2 += 1
          d2 += point[1]
          
     # Construct the system
     A = np.array([[a1,b1],[a2,b2]])
     b = np.array([[d1],[d2]])
     
     # Solve the system
     x = np.dot(np.linalg.inv(A),b)
     
     # Plot
     xp = np.linspace(np.min([data[i][0] for i in range(len(data))]), \
                      np.max([data[i][0] for i in range(len(data))]))
     yp = x[0]*xp + x[1]
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
a, b = lsLine(data)