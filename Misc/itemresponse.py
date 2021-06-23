#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 23 12:47:19 2021

@author: dileepn\

Applied Calculus: Item repsonse curve
"""
import numpy as np
import matplotlib.pyplot as plt

a = 1
b = -1
c = 0

theta = np.linspace(-5,5,100)
y1 = np.ones(len(theta))

P = c + (1-c)/(1+np.exp(-a*(theta-b)))

plt.figure()
plt.plot(theta,P,'r',linewidth=2)
plt.plot(theta,y1,'--',linewidth=1)
plt.plot(theta,np.repeat(np.min(P),len(theta)),'--',linewidth=1)
plt.grid()