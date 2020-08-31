# -*- coding: utf-8 -*-
"""
Spyder Editor

Gaussian Elimination: A Toy Example (for square systems of dimension 4 or lower)
"""
import numpy as np
import string
import sys

def gaussElim(A, b, tol=1e-6):
          """ 
          Inputs:
          - A (real, numpy array): Matrix of coefficients (has to be a square matrix of dimension 4 or lower)
                         Enter A as a numpy array
                         (e.g. The 2x2 matrix [1 2 
                                               4 5] would be entered
                          as np.array([[1 2], [4 5]])
          - b (real, numpy array): Coefficients on the right hand side in Ax = b
                         Enter b as a numpy array (column vector) 
                         (e.g. The 3x1 column [1 
                                               2 
                                               3] 
                          would be entered as np.array([[1],[2], [3]])
          
          Output: 
          - The final form of the augmented matrix, [U|b]
          - Solution for the system, if it exists
          
          Description:
          This function performs Gaussian elimination on square matrices of 
          dimension 3 or lower and calculates the soluion for the system Ax = b provided 
          the matrix is non-singular. It has the capability to do row exchanges 
          whenever they are necessary. The default tolerance used to check if a quanitity 
          is zero is 1x10^-6 i.e. if the absolute value of an entry is less than 1x10^-6, 
          it is considered to be equal to zero. Also, row exchanges will only 
          be done if the element below is greater than the tolerance. To change 
          the tolerance, pass it as the second argument to the function 
          (e.g. gaussElim(3,1e-4)).
          """
          ################## Inputs ##################
          
          ### Input error handling ###
          while True:
               
               # Check if A is okay
               try:
                    # Get dimension
                    dim = A.shape[0]
                    # Check if dimension is okay
                    if dim < 2 or dim > 4:
                         print("Dimension of A too low or too high. Check function help and try again.")
                         return
                    # Check if matrix entries are okay
                    for i in range(dim):
                         for j in range(dim):
                              A[i,j] = float(A[i,j])
               except:
                    print("Invalid input for A. Check help for accepted input types and size.")
                    inp = input("Do you want to try again? y or n: \n")
                    if inp == "y":
                         continue
                    if inp == "n":
                         print("Goodbye!")
                         return 
                    else:
                         print("I didn't understand that. Goodbye!")
                         return
               
               # Check if b is okay
               try:
                    # Chack b's shape
                    if not (b.shape[1] == 1 or b.shape[0] == dim):
                         print("Dimension of b too low or too high. Check function help and try again.")
                    for i in range(dim):
                         b[i] = float(b[i])
                    break
               except:
                    print("Invalid input for b. Check help for accepted input types.")
                    inp = input("Do you want to try again? y or n: \n")
                    if inp == "y":
                         continue
                    if inp == "n":
                         print("Goodbye!")
                         return
                    else:
                         print("I did not understand that. Goodbye!")
                         return
                    
          ################## Calculations ##################
          
          # Augmented matrix (form: [AUG|b])
          AUG = np.concatenate((A,b),axis=1)
          print("\nAugmented Matrix is:\n", AUG)
          print("\n")
          
          ### Start looping row-wise ###
          
          # Remember the current pivot
          curr_pivot = 0
          
          # Start loop
          for i in range(AUG.shape[0] - 1):
               
                    # Check if pivot is non-zero
                    if abs(AUG[i,i]) < tol: 
                         
                         # If it is zero, try exchanging rows
                         if i == 0:
                              # Pivot in row 1 is zero
                              if abs(AUG[i+1,i]) > tol:
                                   AUG[[i,i+1],:] = AUG[[i+1,i],:]
                              elif abs(AUG[i+2,i]) > tol:                         
                                   AUG[[i,i+2],:] = AUG[[i+2,i],:]
                              elif dim > 3 and abs(AUG[i+3,i]) > tol:                        
                                   AUG[[i,i+3],:] = AUG[[i+3,i],:]
                              else:
                                   print("Reduced matrix:\n",AUG)
                                   print("\nI couldn't find a pivot. Possible singularity.\n")
                                   return
                         elif i == 1:
                              # Pivot in row 2 is zero
                              if abs(AUG[i+1,i]) > tol:
                                   AUG[[i,i+1],:] = AUG[[i+1,i],:]
                              elif dim > 3 and abs(AUG[i+2,i]) > tol:
                                   AUG[[i,i+2],:] = AUG[[i+2,i],:]
                              else:
                                   print("Reduced matrix:\n",AUG)
                                   print("\nI couldn't find a pivot. Possible singularity.\n")
                                   return
                         elif i == 2:
                              # Pivot in row 3 is zero
                              if abs(AUG[i+1,i]) > tol:
                                   AUG[[i,i+1],:] = AUG[[i+1,i],:]
                              else:
                                   print("Reduced matrix:\n",AUG)
                                   print("\nI couldn't find a pivot. Possible singularity.\n")
                                   return
          
                    
                    # The current pivot
                    curr_pivot = AUG[i,i]
                    print("Pivot for row {} is {:.2f}\n".format(i+1, curr_pivot))
                    
                    # Do row operations
                    zero_pivot = (abs(curr_pivot) < tol)
                    if zero_pivot:
                         print("At least one pivot is zero. The matrix is singular.\n")
                         print("Current form of the matrix:\n",AUG)
                         return
                    elif i < 1:
                         if abs(AUG[i+1,i]) > tol: 
                              AUG[i+1,:] = AUG[i+1,:] - (AUG[i+1,i]/curr_pivot)*AUG[i,:]
                         if dim > 2 and abs(AUG[i+2,i]) > tol: 
                              AUG[i+2,:] = AUG[i+2,:] - (AUG[i+2,i]/curr_pivot)*AUG[i,:]
                         if dim > 3 and abs(AUG[i+3,i]) > tol: 
                              AUG[i+3,:] = AUG[i+3,:] - (AUG[i+3,i]/curr_pivot)*AUG[i,:]
                    elif i < 2:
                         if abs(AUG[i+1,i]) > tol: 
                              AUG[i+1,:] = AUG[i+1,:] - (AUG[i+1,i]/curr_pivot)*AUG[i,:]
                         if dim > 3 and abs(AUG[i+2,i]) > tol: 
                              AUG[i+2,:] = AUG[i+2,:] - (AUG[i+2,i]/curr_pivot)*AUG[i,:]
                    elif i < 3 and dim > 3:
                         if abs(AUG[i+1,i]) > tol: 
                              AUG[i+1,:] = AUG[i+1,:] - (AUG[i+1,i]/curr_pivot)*AUG[i,:]                         
                         
          
          ################## Output ##################
          
          # Matrix output
          print("Pivot for row {} is {:.2f}\n".format(dim-1, AUG[dim-1,dim-1]))
          print("################ OUTPUT ###############\n")               
          print("Upper triangular (augmented) matrix, [U|b]:\n", AUG)
          
          ### Calculate solution by back substitution ###
          
          x = np.zeros(dim)
          # Iterate in reverse for back substitution
          for i in range(dim-1,-1,-1):
               if i == dim-1:
                    x[i] = AUG[i,-1]/AUG[i,i]
               else:
                    x[i] = (AUG[i,-1] - np.dot(AUG[i,i+1:-1],x[i+1:]))/AUG[i,i]
               
          # Print solution
          if dim == 2:
               print("\nSolution is: x = {:.2f}, y = {:.2f}".format(x[0],x[1]))
          elif dim == 3: 
               print("\nSolution is: x = {:.2f}, y = {:.2f}, z = {:.2f}".format(x[0],x[1],x[2]))
          else:     
               print("\nSolution is: x = {:.2f}, y = {:.2f}, z = {:.2f}, t = {:.2f}".format(x[0],x[1],x[2],x[3]))
               
          print("\n################ END OF OUTPUT ###############")   
          return
          