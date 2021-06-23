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

p = np.linspace(2,5,100)

q = 2000/p**3

dq = -6000/p**2

E = dq*p/q

plt.figure()
plt.plot(p,E,'b',linewidth=2)
plt.grid()
