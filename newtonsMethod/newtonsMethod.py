#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Aug  1 14:25:47 2020

@author: dileepn

Newton's method: A toy example
"""
import numpy as np

f = np.exp

def fn(f,x):
     """
     Inputs:
          - f: function definition (e.g. f = exp)
          - x: initial guess or subsequent predicitons
     
     Output:
          - Function value at x i.e. f(x)
     """
     
     return float(map(f,x))

def d_fn(x,h=1e-6):
     """
     Inputs:
          - x: initial guess or subsequent predicitons
          
     Outputs:
          - Derivative of f at x i.e. df(x)/dx
     
     Note: Derivative is computed numerically using the forward difference method
           df(x)/dx = (f(x+h)-f(x))/h
           To change the interval, h, pass it as the second argument to the function
           (e.g. d_fn(x,1e-4))
     """
     if h < 1e-8:
          print("Interval, h, is too low. Resetting to default value.")
          h = 1e-6
     
     return (fn(x+h) - f(x))/h

def newt(f,x,tol=1e-5):
     """
     Inputs:
          - f: function definition (in f(x) = 0 form e.g. if you wish to solve x^2 = 5, f would 
                                    be x^5 - 5)
          - x: initial guess or subsequent predicitons
          
     Output:
          - Best approximation to the solution of f(x) = 0
     
     Note: 
     """

     
     
