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
X = \
"65.2, 21.9,-14.1,29.6,2.1,64.5,5.5,36.7,43.0,6.2,21.6,-1.7,-12.1,22.8,41.7,68.3,64.2,53.0,62.8,-14.5"
Y = \
"-48.2,53.9,86.1,15.1,68.3,-41.7,80.8,2.6,1.8,73.0,53.0,49.3,77.7,27.7,9.5,-37.5,-30.5,-37.1,-28.8,67.2"

dataX = [float(x) for x in X.split(",")]
dataY = [float(y) for y in Y.split(",")]

data = []

for i in range(len(dataX)):
     data.append((dataX[i],dataY[i])) 

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