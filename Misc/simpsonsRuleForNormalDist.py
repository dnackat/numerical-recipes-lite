#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Sep 20 20:19:00 2020

@author: dileepn
"""
import numpy as np

# Ends
a = 96
b = 100

# dx: width of the interval
dx = 2

# Number of intervals
n = int((b-a)/dx)



# Function definition
def f(x):
     return (1./(2.8*np.sqrt(2*np.pi)))*np.exp(-((x-69.)**2)/5.6)

#### Approximate the integral using the Simpson's rule ####

# Variable to hold sum
simp_sum = 0.0

# x_i values
x = np.arange(a,b+1,dx)

for i in range(n+1):
     # Set weight of first and last terms equal to unity
     if i == 0 or i == n:
          simp_sum += f(x[i])
     else:
          # Alternating coeffs. of 2 and 4
          if (i+1) % 2 == 0:
               # Even index
               simp_sum += 4.*f(x[i])
          else:
               # Odd index
               simp_sum += 2.*f(x[i])
               
simp_sum *= dx/3.0

print("\nApproximation for the integral using Simpson's rule is: %10.3E\n" %(simp_sum))