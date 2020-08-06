#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug  6 11:35:52 2020

@author: dileepn

This toy script creates a Pascal's triangle matrix
"""

import numpy as np

def pascal(max_size=100):
     """
     Input:
          n: The desired size of the Pascal's triangle matrix (if input is n, size is n by n)
          max_size: The maximum size of the triangle permitted. 
          
     Output:
          Pascal's triangle matrix of size n by n as a numpy array
          
     Description:
          This function generates a Pascal's triangle matrix of size n by n. An element 
          (i,j) of a Pascal's triangle matrix is "i choose j" (number of ways of picking)
          j elements out of i things. These elements can also be thought as 
          coefficients of a binomial expansion.
     """
     
     # Get input for size
     n = input("Enter the desired size of the Pascal's triangle (positive integer):\n")
     
     # Check if size entered is acceptable
     try:
          n = int(n)
     except ValueError:
          print("Invalid entry for the size, n. Please try again.")
          return
     
     # Generate an array to store the Pascal matrix
     P = np.zeros((n,n))
     
     # Start the loop to fill in entries
     for i in range(n):
          for j in range(i+1):
               
               # First row has just one element filled in
               if i == 0 and j == 0:
                    P[i,j] = 1
               # Use rule for all other elements
               else:
                    P[i,j] = P[i-1,j-1] + P[i-1,j]
                    
     # Output the Pascal matrix
     print("Pascal's triangle matrix of size",n,"is:\n",P)
     
     return
                    