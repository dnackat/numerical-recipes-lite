#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Aug 16 12:40:54 2020

@author: dileepn

This script creates permutation matrices of specified sizes for row or column exchanges
"""

import numpy as np

def permMatrix(dim=3,ex='r'):
     """
     Inputs:
          dim (integer): The size of the permutation matrix (square, so 3 means a 3x3 matrix)
          ex (r or c): r means row exchanges; c means column exchanges
          
     Output:
          The desired permutation matrix
          
     Description:
          This function
     """