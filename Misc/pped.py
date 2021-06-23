#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 23 12:45:50 2021

@author: dileepn

Applied Calculus: Point price elasticity

"""
import numpy as np
import matplotlib.pyplot as plt

p = np.linspace(0,5.5,100)

E = -29*p/(-29*p + 166)

plt.figure()
plt.plot(p,E,'b',linewidth=2)
plt.grid()

#%% Price point elasticity using a nonlinear function

# Price
p = np.linspace(2,5,100)

# Demand function
q = 2000/p**3

# Derivative of demand function
dq = -6000/p**2

# Price point elasticity of demand
E = dq*p/q

# Revenue function
R = p*q

# Derivative of revenue function
dR = q + p*dq

plt.figure()
plt.plot(p,R,'b',linewidth=2)
plt.grid()
