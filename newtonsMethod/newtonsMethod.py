#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Aug  1 14:25:47 2020

@author: dileepn

Newton's method: A toy example
"""
import numpy as np

f = lambda x: np.exp(x) + x

def fn(f,x):
     """
     Inputs:
          - f: function definition (e.g. f = exp)
          - x: initial guess or subsequent predicitons
     
     Output:
          - Function value at x i.e. f(x)
     """
     
     return f(x)

def d_fn(f,x,h=1e-6):
     """
     Inputs:
          - f: function definition (e.g. f = x**2 - 5)
          - x: initial guess or subsequent predicitons
          
     Outputs:
          - Derivative of f at x i.e. df(x)/dx
     
     Note: 
          Derivative is computed numerically using the forward difference method
          df(x)/dx = (f(x+h)-f(x))/h
          To change the interval, h, pass it as the second argument to the function
          (e.g. d_fn(x,1e-4))
     """
     
     if h < 1e-8:
          print("Interval, h, is too low. Setting to 1e-8.\n")
          h = 1e-8
     
     return (fn(f,x+h) - fn(f,x))/h

def newt(f,x0,tol=1e-5,iters=1000):
     """
     This function solves the input equation iteratively using Newton's method.
     x_i+1 = x_i - f(x_i)/f'(x_i) 
     
     Inputs:
          - f: function definition (in f(x) = 0 form e.g. if you wish to solve x^2 = 5, f would 
                                    be x^5 - 5). While entering function definition,
               use lambda method e.g. for e^x + x = 0, f = lambda x: np.exp(x) + x 
          - x0: initial guess or subsequent predicitons
          
     Output:
          - Best approximation to the solution of f(x) = 0
     
     Notes:
          
          - Default tolerance used to check convergence is 1e-5. Change it by passing
            it as the third argument to the function (e.g. newt(f,x0,1e-6,iters,span)).
          - By default this function uses at most 100 iterations to check if 
            there is convergence. If you wish to change the number of iterations,
            pass it as the fourth argument to the function (e.g. newt(f,x0,tol,500,span)) 
     """
     if tol < 1e-5:
          print("Interval, h, is too low. Resetting to default vallue of 1e-5.\n")
          tol = 1e-5
          
     try:
          iters = int(iters)
     except ValueError:
          print("Invalid value for ietrations. Try again.\n")
          return

     # Array to hold initial guess and successive predictions          
     x = np.zeros(iters)
     
     # Start loop
     for i in range(iters):
          if i == 0:
               x[i] = x0
          else:
               x[i] = x[i-1] - fn(f,x[i-1])/d_fn(f,x[i-1])
               
          # Check for convergence
          if i > 0 and abs(x[i] - x[i-1]) < tol:
               print("\nThe solution is: {:.4f}\n".format(x[i]))
               return
          
     print("\nSorry, I couldn't find a solution in",iters,"iterations.")
     
     return
          
     
     
     
     
