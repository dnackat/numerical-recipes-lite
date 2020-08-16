#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Aug 16 12:40:54 2020

@author: dileepn

This script creates permutation matrices of specified sizes for row or column exchanges
"""

import numpy as np

def permMatrix(dim=3,ex='r',order=[0]):
     """
     Inputs:
          dim (integer): The size of the permutation matrix (square, so 3 means a 3x3 matrix)
          ex (r or c): r means row exchanges; c means column exchanges
          order (list of integers): The order in which to shuffle rows or columns.
                                   e.g. 4321 means 4th row is switched with 1st row,
                                        3rd row is switched with the 2nd row, etc.
                                   Note that the length of this list must be equal to 'dim'.
          
          
     Output:
          The desired permutation matrix
          
     Description:
          This function creates a permutation matrix of the desired size for row
          or column exchanges. This function is not optimized for resource utilization
          as an identity matrix is created each time an exchange is to be done.
          This is done to keep the script remains short and simple.
     """
     
     # Check if dimension okay
     try:
          dim = int(dim)
     except ValueError:
          print("\nInvalid entry for dimension. Please enter an integer.\n")
          return
     
     # Check order list
     for i in range(len(order)):
          try:
               int(order[i])
          except ValueError:
               print("Invalid 'order' list! It must contain integers. Try again.")
     
     # Check if 'ex' is okay
     if ex == 'r':
          print("\nRow exchange option selected.\n")
     elif ex == 'c':
          print("\nColumn exchange option selected.\n")
     else:
          print("\nInvalid entry for exchange paramter. Check function help and try again.\n")
          return

     
     # Create an array to hold the permutation matrix
     P = np.identity(dim)
     
     # Build the desired permutation matrix
     if len(order) == 1 and order[0] == 0:
          # Do nothing. Just return the identity matrix
          print("\nNo exchanges performed.\n")
          return P
     elif len(order) > 1:
          # Iterate through 'order' list and only populate rows/cols that need to be exchanged
          for i in range(len(order)):
               if not (order[i] - 1 == i):
                    if ex == 'r':
                    # Row exchanges
                         P[i,:] = np.identity(dim)[order[i]-1,:]
                    else:
                    # Column exchange
                         P[:,i] = np.identity(dim)[:,order[i]-1]
                         
     # If all is well, return the desired perm matrix
     return P