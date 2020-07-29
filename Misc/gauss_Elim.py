# -*- coding: utf-8 -*-
"""
Spyder Editor

Gaussian Elimination: A Toy Example (3x3 system)
"""
import numpy as np
import sys

def gaussElim(dim):
          """ 
          Inputs:
          - dim (integer): Dimension of the matrix (if you enter 3, it means a 3x3 matrix)
          - A (real): Matrix of coefficients (has to be a square matrix and 3x3 or lower)
                         Enter A as a series of numbers separated by single spaces 
                         starting from A[1,1] and proceeding row-wise
                         (e.g. The 2x2 matrix [1 2 
                                               4 5] would be entered
                          as 1 2 4 5)
          - b (real): Coefficients on the right hand side in Ax = b
                         Enter b as a series of numbers starting from b[1,1] and 
                         proceeding down the column 
                         (e.g. The 3x1 column [1 
                                               2 
                                               3] 
                          would be entered as 1 2 3)
          
          Output: 
          - The final form of the augmented matrix, [U|b]
          - Solution for the system, if it exists
          
          Description:
          This function performs Gaussian elimination on square matrices of 
          dimension 3 or lower and calculates the soluion for the system provided 
          the matrix is non-singular. It has the capability to do row exchanges 
          whenever they are necessary. The tolerance used to check if a quanitity 
          is zero is 1x10^-8 i.e. if the absolute value of an entry is less than 1x10^-8, 
          it is considered to be equal to zero. Also, row exchanges will only 
          be done if the element below is greater than 1x10^-5. 
          """
          ### Input handling ###
          
          # prompt for input. We want A and b
          A = input("Enter matrix A ({}x{}) as a list: ".format(dim,dim))
          b = input("Enter column vector b ({}x1) as a list: ".format(dim))
          
          ### Convert matrix and column vector to numpy arrays ###
          
          # Matrix of coefficients
          A = np.array(A).reshape((dim, dim))
          
          # Right hand side coefficients (form: Ax = b)
          b = np.array(b).reshape((dim,1))
          
          # Augmented matrix (form: [AUG|b])
          AUG = np.concatenate((A,b),axis=1)
          print("Augmented Matrix is:\n", AUG)
          
          
          ### Start looping row-wise ###
          # Remember the current pivot and have a temp row for exchanges
          curr_pivot = 0
          temp_row = np.zeros(AUG.shape[1])
          
          # Start loop
          for i in range(AUG.shape[0]):
                    # Check if pivot is non-zero
                    if i < AUG.shape[0]-1 and abs(AUG[i,i]) < 1e-8:
                         
                         # Try exchanging rows
                         if i == 0:
                              if abs(AUG[i+1,i]) > 1e-5:
                                   temp_row = AUG[i,:].copy()
                                   AUG[i,:] = AUG[i+1,:]
                                   AUG[i+1,:] = temp_row
                              elif abs(AUG[i+2,i]) > 1e-5:
                                   temp_row = AUG[i,:].copy()                        
                                   AUG[i,:] = AUG[i+2,:]
                                   AUG[i+2,:] = temp_row
                              else:
                                   print("I couldn't find a pivot. Possible singularity.")
                                   sys.exit(0)
                         elif i == 1:
                              if abs(AUG[i+1,i]) > 1e-5:
                                   temp_row = AUG[i,:].copy()
                                   AUG[i,:] = AUG[i+1,:]
                                   AUG[i+1,:] = temp_row
                              else:
                                   print("I couldn't find a pivot. Possible singularity.")
                                   sys.exit(0)
          
                    
                    # The current pivot
                    curr_pivot = AUG[i,i]
                    print("Pivot for row",i+1,"is:",curr_pivot)
                    
                    # Do row operations
                    zero_pivot = (abs(curr_pivot) < 1e-8)
                    if zero_pivot:
                         print("At least one pivot is zero. The matrix is singular.")
                         sys.exit(0)
                    elif i < 1:
                         if abs(AUG[i+1,i]) > 1e-5: 
                              AUG[i+1,:] = AUG[i+1,:] - (AUG[i+1,i]/curr_pivot)*AUG[i,:]
                         if abs(AUG[i+2,i]) > 1e-5: 
                              AUG[i+2,:] = AUG[i+2,:] - (AUG[i+2,i]/curr_pivot)*AUG[i,:]
                    elif i < 2:
                         if abs(AUG[i+1,i]) > 1e-5: 
                              AUG[i+1,:] = AUG[i+1,:] - (AUG[i+1,i]/curr_pivot)*AUG[i,:]
          
          ### Output ###
          
          # Matrix output               
          print("Upper triangular (augmented) matrix, [U|b]:\n", AUG)
          
          # Calculate solution
          z = AUG[AUG.shape[0]-1,AUG.shape[1]-1]/AUG[AUG.shape[0]-1,AUG.shape[1]-2]
          y = (AUG[AUG.shape[0]-2,AUG.shape[1]-1] - \
               z*AUG[AUG.shape[0]-2,AUG.shape[1]-2])/AUG[AUG.shape[0]-2,AUG.shape[1]-3]
          x = (AUG[AUG.shape[0]-3,AUG.shape[1]-1] - z*AUG[AUG.shape[0]-3,AUG.shape[1]-2] \
               - y*AUG[AUG.shape[0]-3,AUG.shape[1]-3])/AUG[AUG.shape[0]-3,AUG.shape[1]-4]
               
          # Print solution
          print("Solution is: x = {:.2f}, y = {:.2f}, z = {:.2f}".format(x,y,z))
               
          return
          