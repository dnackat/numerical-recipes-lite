#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug  6 11:35:52 2020

@author: dileepn

This toy script creates a Pascal's triangle 
"""

import numpy as np

def pascal(max_size=100):
     """
     Input:
          n: The desired size of the Pascal's triangle (if input is n, size is n by n)
          max_size: The maximum size of the triangle permitted. 
          
     Output:
          Pascal's triangle of size n by n as a numpy array
          
     Description:
          This function generates a Pascal's triangle of size n by n. An element 
          (i,j) of a Pascal's triangle is "i choose j" (number of ways of picking)
          j elements out of i things. These elements can also be thought as 
          coefficients of a binomial expansion.
     """
     n = input("Enter the desired size of the Pascal's triangle")
     try:
          n = int(n)
