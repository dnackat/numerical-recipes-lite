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
                                   Note that the size of the list must be equal to 'dim'.
          
          
     Output:
          The desired permutation matrix
          
     Description:
          This function creates a permutation matrix of the desired size for row
          or column exchanges.
     """
     
     # Check if dimension okay
     try:
          dim = int(dim)
     except ValueError:
          print("\nInvalid entry for dimension. Please enter an integer.\n")
          return
     
     # Check if 'ex' is okay
     if ex == 'r':
          print("\nDoing row exchanges.\n")
     elif ex == 'c':
          print("\nDoing column exchanges.\n")
     else:
          print("\nInvalid entry for exchange paramter. Check function help and try again.\n")
          return
     
     # Create an array to hold the permutation matrix
     P = np.identity(dim)
     
     # Build the desired permutation matrix
     if len(order) == 1 and order[0] == 0:
          # Do nothing. Just return the identity matrix
          return P
     elif len(order) > 1:
          # Iterate through 'order' list and only populate rows/cols that need to be exchanged
          for i in range(len(order)):
               if not (order[i] - 1 == i):
                    if ex == 'r':
                    # Row exchanges
                         P[i,:] = np.identity(dim)[order[i]-1,:]