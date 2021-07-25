#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jul 25 19:50:15 2021

@author: dileepn

Contour plots for multivariable calculus course
"""

# Preliminaries
import numpy as np
import matplotlib.pyplot as plt

# Create a grid
x = np.linspace(-3,3,60)
y = np.linspace(-3,3,60)

X, Y = np.meshgrid(x,y)

# Function definitions
def f1(x,y):
     return 2*x + y

def f2(x,y):
     return x**2 + y**2 

def f3(x,y):
     return x**2 - y**2

# Calculate function values
Z1 = f1(X,Y)
Z2 = f2(X,Y)
Z3 = f3(X,Y)

# Create plots
