# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

#%%
v = np.linspace(0,100,1000)
c = 3*1e5 
T1 = 1/np.sqrt(1-v**2/c**2)
T2 = 1 + (1/2)*v**2/c**2 


plt.figure()
plt.grid(axis='both')
plt.plot(v,T1,'r',label='T1')
plt.plot(v,T2,'k',label='T2')

plt.legend(loc='best')

#%% 3D Plot attempt
#
#x = np.linspace(-2,2,100)
#y = np.linspace(-2,2,100)
#
#xx, yy = np.meshgrid(x,y)
#
#zz = np.exp(-(xx**2 + yy**2))
#
#fig = plt.figure()
#ax = Axes3D(fig)
#ax.plot_surface(xx,yy,zz,cmap='winter')