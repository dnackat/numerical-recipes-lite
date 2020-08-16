#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Aug 16 12:40:54 2020

@author: dileepn

This script creates permutation matrices of specified sizes for row or column exchanges
"""

import numpy as np

def permMatrix(dim=3,ex='r',order=(1,2,3,4)):
     """
     Inputs:
          dim (integer): The size of the permutation matrix (square, so 3 means a 3x3 matrix)
          ex (r or c): r means row exchanges; c means column exchanges
          
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
     
     # Create 