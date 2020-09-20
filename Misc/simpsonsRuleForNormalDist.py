#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Sep 20 20:19:00 2020

@author: dileepn
"""
import numpy as np

# Ends
a = 0
b = 96

# Number of intervals
n = 8

# dx
dx = (b-a)/n

# Function definition
def f(x):
     return (1./(2.8*np.sqrt(2*np.pi)))*np.exp(-(x-69)**2/5.6)

#### Approximate the integral using the Simpson's rule ####

# Variable to hold sum
simp_sum = 0.0

# x_i values
x = np.arange(a,b+1,n)

for i in range(n):
     # Set weight of first and last terms equal to unity
     if i == 0 or i == (n-1):
         simp_sum += f(x(i))
     else:
          if (i+1) % 2 