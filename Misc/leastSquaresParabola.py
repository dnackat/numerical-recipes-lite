#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 26 13:00:54 2021

@author: dileepn

Best fit parabola using least-squares approximation
"""
# Preliminaries
import numpy as np
import matplotlib.pyplot as plt

# Data
data = [(0,1),(2,1),(4,3)]

# Get all the coefficients for the linear 3x3 system in a, b, and c
a1 = 0.; b1 = 0.; c1 = 0.; a2 = 0.; b2 = 0.; c2 = 0.; a3 = 0.; b3 = 0.; c3 = 0.
d1 = 0.; d2 = 0.; d3 = 0.;

for point in data:
     # First row
     a1 += point[0]**4
     b1 += point[0]**3
     c1 += point[0]**2
     d1 += point[1]*point[0]**2
     # Second row
     a2 += point[0]**3
     b2 += point[0]**2
     c2 += point[0]
     d2 += point[1]*point[0]
     # Third row
     a3 += point[0]**2
     b3 += point[0]
     c3 += 1
     d3 += point[1]